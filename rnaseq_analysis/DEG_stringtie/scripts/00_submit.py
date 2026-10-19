### 脚本名称：00_submit.py
### 功能：冻结代码与配置，为每材料生成独立Slurm脚本并记录作业ID；支持不提交预览。
### 输入文件：--input 材料信息表；输出文件：--outdir 下代码快照、配置、作业与日志。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python scripts/00_submit.py --input config/materials/sample_info.txt --outdir runs/批次 --no-submit。
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape
import pandas as pd
import importlib
io = importlib.import_module('00_common')

def main():
    p = io.parser('生成并提交逐材料StringTie差异分析作业')
    p.add_argument('--no-submit',action='store_true')
    args = p.parse_args()
    io.require(not Path(args.outdir).exists(), '批次目录已存在；请用新批次名，避免重复提交或覆盖')
    out, log = io.start(args); out = out.resolve()
    root = Path(__file__).resolve().parent.parent
    code = out/'code'; code.mkdir()
    for item in ['scripts','config','assets']:
        shutil.copytree(root/item,code/item,ignore=shutil.ignore_patterns('__pycache__','materials'))
    for item in ['main.sh','submit_stringtie.sh','README.md']:
        shutil.copy2(root/item,code/item)
    io.require((code/'assets/Arial.ttf').is_file(), '缺少Arial字体')
    (code/'assets/fonts.conf').write_text('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig><dir>'+escape(str(code/'assets'))+'</dir><cachedir>'+escape(str(out/'fontcache'))+'</cachedir></fontconfig>')
    hashes = [{'path':str(x.relative_to(code)),'sha256':io.digest(x)} for x in sorted(code.rglob('*')) if x.is_file()]
    io.write(pd.DataFrame(hashes),out/'00_code_sha256.tsv')
    records = io.read(args.input)
    io.require(set(records.columns)=={'material','config_dir'} and not records.material.duplicated().any(), '材料表错误')
    generated = []
    for row in records.itertuples():
        io.require(re.fullmatch('[A-Za-z0-9_.-]+',row.material), '材料名不合法')
        config = out/'configs'/row.material; shutil.copytree(row.config_dir,config)
        work = out/row.material; work.mkdir()
        slurm = work/f'{row.material}.slurm'
        q = shlex.quote
        lines = ['#!/bin/bash',f'#SBATCH -J stDE_{row.material}',f'#SBATCH -p {os.environ.get("PARTITION","x86_64_4cpu")}',
                 '#SBATCH -N 1',f'#SBATCH -n {int(os.environ.get("THREADS","2"))}',
                 f'#SBATCH -o {work}/slurm.%j.out',f'#SBATCH -e {work}/slurm.%j.err',
                 '### 功能：单材料StringTie全基因差异分析；环境与参数来自冻结配置。','set -euo pipefail']
        for key in ['PIPE_PYTHON','PIPE_RSCRIPT','MIN_COUNT','MIN_SAMPLES','FDR','ABS_LOG2FC']:
            if key in os.environ:
                lines.append(f'export {key}={q(os.environ[key])}')
        lines += [f'cd {q(str(code))}',f'bash {q(str(code/"main.sh"))} {q(str(config))} {q(str(work/"results"))} all']
        slurm.write_text('\n'.join(lines)+'\n')
        subprocess.run(['bash','-n',str(slurm)],check=True)
        generated.append({'material':row.material,'script':str(slurm),'job_id':''})
    io.write(pd.DataFrame(generated),out/'00_submission.tsv')
    ### 全部脚本通过语法检查后才提交；每次提交立即落盘，以便中途失败时恢复。
    for record in generated:
        log('计划：sbatch '+record['script'])
        if not args.no_submit:
            result = subprocess.run(['sbatch','--parsable',record['script']],check=True,text=True,capture_output=True)
            record['job_id'] = result.stdout.strip().split(';')[0]
            io.require(record['job_id'].isdigit(), f'无法识别作业ID：{result.stdout}')
            io.write(pd.DataFrame(generated),out/'00_submission.tsv')
            log(f'{record["material"]} 已提交：{record["job_id"]}')
    return args

if __name__ == '__main__':
    io.run(main)
