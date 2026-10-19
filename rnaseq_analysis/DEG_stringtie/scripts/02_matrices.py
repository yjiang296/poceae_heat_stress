### 脚本名称：02_matrices.py
### 功能：按真实样本ID截取同源StringTie全基因counts和TPM，拒绝缺失、重复及非法数值。
### 输入文件：--input 为01输出目录；输出文件：02全基因矩阵、元数据、QC与校验和。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python scripts/02_matrices.py --input results/B73/01_audit --outdir results/B73/02_matrices。
from pathlib import Path
import numpy as np
import pandas as pd
import importlib
io = importlib.import_module('00_common')

def main():
    args = io.parser('构建统一StringTie输入矩阵').parse_args()
    out, log = io.start(args)
    root = Path(args.input)
    io.require((root/'01_SUCCESS').exists(), '缺少01结果，请先运行01或all')
    samples, sources = io.read(root/'01_samples.tsv'), io.read(root/'01_sources.tsv').iloc[0]
    hashes = io.read(root/'01_input_sha256.tsv').set_index('path').sha256
    frames, mappings = {}, []
    for metric in ['counts','tpm']:
        path = sources[metric]
        io.require(io.digest(path)==hashes[path], f'源矩阵在审核后发生变化：{path}')
        tab = io.read(path)
        ### 仅去除已知GFF风格gene:前缀；保留原ID映射，碰撞时立即报错。
        original_ids = tab.gene_id.copy()
        tab['gene_id'] = tab.gene_id.str.removeprefix('gene:')
        mappings.extend({'source':metric,'original_gene_id':a,'gene_id':b} for a,b in zip(original_ids,tab.gene_id))
        io.require(tab.columns[0]=='gene_id' and not tab.gene_id.duplicated().any() and (tab.gene_id!='').all(), f'{metric}基因ID错误')
        io.require(set(samples.sample_id)<=set(tab.columns), f'{metric}缺少样本列')
        tab = tab.set_index('gene_id')[list(samples.sample_id)].apply(pd.to_numeric)
        values = tab.to_numpy()
        io.require(np.isfinite(values).all() and (values>=0).all(), f'{metric}有负数或非有限数')
        if metric=='counts':
            io.require((values==np.floor(values)).all() and values.max()<=2147483647, 'counts必须为合法非负整数')
            io.require((values.sum(axis=0)>0).all(), '有样本counts总和为0')
            tab = tab.astype(np.int64)
        frames[metric] = tab
    counts, tpm = frames['counts'], frames['tpm']
    io.require(set(counts.index)==set(tpm.index), 'counts与TPM基因集合不同；禁止隐式补零')
    tpm = tpm.reindex(counts.index)
    qc = []
    for sample in samples.sample_id:
        discordance = (counts[sample]==0)&(tpm[sample]>=1)
        qc.append({'sample_id':sample,'estimated_count_sum':int(counts[sample].sum()),
                   'genes_count_ge10':int((counts[sample]>=10).sum()),'zero_count_tpm_ge1':int(discordance.sum())})
    io.write(pd.DataFrame(qc), out/'02_quantification_qc.tsv')
    io.write(pd.DataFrame(mappings), out/'02_gene_id_mapping.tsv')
    io.require(sum(x['zero_count_tpm_ge1'] for x in qc)==0, '同源counts与TPM仍有count=0且TPM>=1，检查上游批次')
    io.write(counts.reset_index(), out/'02_gene_counts.tsv')
    io.write(tpm.reset_index(), out/'02_gene_tpm.tsv')
    samples['material'] = sources.material
    io.write(samples, out/'02_samples.tsv')
    io.write(io.read(root/'01_contrasts.tsv'), out/'02_contrasts.tsv')
    (out/'02_SUCCESS').write_text('矩阵检查通过\n')
    log(f'矩阵完成：{len(counts)}个基因、{len(samples)}个样本；counts来自StringTie覆盖度估计')
    return args

if __name__ == '__main__':
    io.run(main)
