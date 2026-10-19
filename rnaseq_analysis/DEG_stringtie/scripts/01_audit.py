### 脚本名称：01_audit.py
### 功能：核验StringTie表达估计模式、BAM和参考注释，冻结输入来源及读长信息。
### 输入文件：--input 材料配置目录；输出文件：01配置副本、来源校验和与定量命令表。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python scripts/01_audit.py --input config/materials/B73 --outdir results/B73/01_audit。
import json
import shlex
from pathlib import Path
import pandas as pd
import importlib
io = importlib.import_module('00_common')

def main():
    args = io.parser('核验StringTie输入与来源').parse_args()
    out, log = io.start(args)
    root = Path(args.input)
    samples = io.read(root/'00_samples.tsv')
    sources = io.read(root/'00_sources.tsv')
    contrasts = io.read(root/'00_contrasts.tsv')
    io.require(len(sources)==1 and sources.method.iloc[0]=='STRINGTIE_PREPDE_FORMULA', '定量来源必须为StringTie')
    io.require(not samples.sample_id.duplicated().any(), '样本重复')
    io.require(samples.annotation_gtf.nunique()==1, '材料内部注释不一致')
    io.require(samples.groupby('condition').size().min() >= 2, '每条件至少两个独立重复')
    io.require(set(contrasts.numerator)|set(contrasts.denominator) == set(samples.condition), '比较与条件不匹配')
    hashes = []
    paths = list(root.glob('00_*.tsv')) + [Path(sources.counts.iloc[0]), Path(sources.tpm.iloc[0]), Path(samples.annotation_gtf.iloc[0])]
    quant = Path(sources.quant_root.iloc[0])
    paths += list((quant/'scripts').glob('05_*counts.py')) + list((quant/'scripts').glob('05_run_prepde.py'))
    commands, total_bases, total_reads = [], 0, 0
    for row in samples.itertuples():
        sample_dir = Path(row.bam).parent.parent
        gtf = sample_dir/'03_stringtie'/f'{row.sample_id}.stringtie.gtf'
        abundance = sample_dir/'03_stringtie'/f'{row.sample_id}.gene_abund.tsv'
        qc = sample_dir/'01_fastp'/f'{row.sample_id}.fastp.json'
        for path in [Path(row.bam),gtf,abundance,qc]:
            io.require(path.is_file() and path.stat().st_size>0, f'输入缺失：{path}')
        with gtf.open() as handle:
            command, version = handle.readline().strip(), handle.readline().strip()
        tokens = shlex.split(command.lstrip('# '))
        io.require('-e' in tokens and '-G' in tokens, f'不是固定注释表达估计模式：{gtf}')
        io.require(tokens[tokens.index('-G')+1] == row.annotation_gtf, f'注释路径不符：{gtf}')
        io.require(row.bam in tokens, f'BAM来源不符：{gtf}')
        after = json.loads(qc.read_text())['summary']['after_filtering']
        total_bases += int(after['total_bases']); total_reads += int(after['total_reads'])
        commands.append({'sample_id':row.sample_id,'gtf':str(gtf),'command':command,'version':version,
                         'gtf_bytes':gtf.stat().st_size,'gtf_mtime_ns':gtf.stat().st_mtime_ns,
                         'average_clean_read_length':int(after['total_bases'])/int(after['total_reads'])})
        paths.append(qc)
    for path in paths:
        io.require(path.is_file(), f'输入缺失：{path}')
        hashes.append({'path':str(path),'sha256':io.digest(path)})
    io.write(samples, out/'01_samples.tsv'); io.write(contrasts, out/'01_contrasts.tsv')
    io.write(sources, out/'01_sources.tsv'); io.write(pd.DataFrame(hashes), out/'01_input_sha256.tsv')
    io.write(pd.DataFrame(commands), out/'01_quantification_provenance.tsv')
    io.write(pd.DataFrame([{'average_clean_read_length':max(1,round(total_bases/total_reads)),
                           'count_method':'ceil(cov * exon_length / read_length), summed by gene',
                           'note':'复用上游公式实现的矩阵；非本流程重新运行官方prepDE'}]), out/'01_count_method.tsv')
    (out/'01_SUCCESS').write_text('输入来源通过\n')
    log(f'核验完成：{len(samples)}个StringTie样本')
    return args

if __name__ == '__main__':
    io.run(main)
