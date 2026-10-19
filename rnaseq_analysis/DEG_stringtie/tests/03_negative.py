### 脚本名称：03_negative.py
### 功能：验证counts与TPM冲突、ID碰撞及缺少上游时会明确拒绝。
### 输入文件：--input 模拟测试工作目录；输出文件：--outdir 异常用例日志和通过记录。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python tests/03_negative.py --input tests/work --outdir tests/negative。
import sys
from pathlib import Path
import shutil
import subprocess
import importlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
io = importlib.import_module('00_common')

def main():
    args = io.parser('验证错误输入拒绝').parse_args()
    out, log = io.start(args); out = out.resolve()
    root = Path(__file__).resolve().parents[1]
    work = Path(args.input).resolve()
    for case in ['discordance','collision']:
        audit = out/case/'audit'; shutil.copytree(work/'results/01_audit',audit)
        sources = io.read(audit/'01_sources.tsv')
        original = sources.tpm.iloc[0]; tab = io.read(original)
        if case=='discordance':
            tab.loc[tab.gene_id=='all_zero',tab.columns[1]] = '100'
        else:
            tab.loc[1,'gene_id'] = 'gene:'+tab.loc[0,'gene_id']
        tpm = out/case/'invalid_tpm.tsv'; io.write(tab,tpm)
        sources.loc[0,'tpm'] = str(tpm); io.write(sources,audit/'01_sources.tsv')
        hashes = io.read(audit/'01_input_sha256.tsv')
        hashes.loc[hashes.path==original,['path','sha256']] = [str(tpm),io.digest(tpm)]
        io.write(hashes,audit/'01_input_sha256.tsv')
        res = subprocess.run([sys.executable,str(root/'scripts/02_matrices.py'),'--input',str(audit),'--outdir',str(out/case/'result')],capture_output=True,text=True)
        (out/f'03_{case}.log').write_text(res.stdout+res.stderr)
        io.require(res.returncode!=0,f'错误输入未被拒绝：{case}')
    res = subprocess.run(['bash',str(root/'main.sh'),str(work/'config'),str(out/'missing_upstream'),'03'],capture_output=True,text=True)
    (out/'03_missing_upstream.log').write_text(res.stdout+res.stderr)
    io.require(res.returncode!=0 and '缺少步骤02' in res.stderr,'上游检查失败')
    (out/'03_NEGATIVE_PASS').write_text('三项错误输入正确拒绝\n')
    log('错误输入拒绝测试通过')
    return args

if __name__=='__main__':
    io.run(main)
