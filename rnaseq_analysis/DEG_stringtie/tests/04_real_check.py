### 脚本名称：04_real_check.py
### 功能：核验正式批次的完整基因结果、响应阈值、图件及三个问题基因的StringTie计数。
### 输入文件：--input 正式批次目录；输出文件：--outdir 材料QC、关键基因counts与统计。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python tests/04_real_check.py --input runs/formal_20261014 --outdir runs/formal_20261014/validation。
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import importlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
io = importlib.import_module('00_common')

def main():
    args = io.parser('验收正式StringTie差异分析').parse_args()
    out, log = io.start(args)
    root = Path(args.input)
    manifest = io.read(root/'00_submission.tsv')
    qc, target_counts, target_stats = [], [], []
    targets = {'B73':['Zm00001eb056980','Zm00001eb384760'],'Bd21':['BRADI_3g44700v3']}
    for item in manifest.itertuples():
        base = root/item.material/'results'
        for folder in ['01_audit','02_matrices','03_deseq2','04_plots']:
            io.require((base/folder/f'{folder[:2]}_SUCCESS').exists(),f'{item.material}：步骤未完成：{folder}')
        counts = io.read(base/'02_matrices/02_gene_counts.tsv').set_index('gene_id')
        quant = io.read(base/'02_matrices/02_quantification_qc.tsv')
        io.require(pd.to_numeric(quant.zero_count_tpm_ge1).sum()==0,'仍有定量冲突')
        contrasts = io.read(base/'03_deseq2/03_contrasts.tsv')
        for cid in contrasts.contrast_id:
            tab = io.read(base/'03_deseq2'/f'03_{item.material}_{cid}_all_genes.tsv')
            io.require(len(tab)==len(counts) and set(tab.gene_id)==set(counts.index),'全基因结果不完整')
            p = pd.to_numeric(tab.padj,errors='coerce'); fc = pd.to_numeric(tab.log2FoldChange,errors='coerce')
            expected = np.where(p.notna()&fc.notna(),'NO_DETECTED_RESPONSE','NOT_EVALUABLE')
            expected[(p<.05)&(fc>=1)] = 'UP'; expected[(p<.05)&(fc<=-1)] = 'DOWN'
            io.require((expected==tab.response).all(),'响应标签与阈值不一致')
            selected = tab[tab.gene_id.isin(targets.get(item.material,[]))].copy()
            if len(selected):
                selected.insert(0,'material',item.material); selected.insert(1,'contrast_id',cid)
                target_stats.extend(selected.to_dict('records'))
        pdfs = list((base/'04_plots').glob('*.pdf'))
        io.require(len(pdfs)==len(contrasts)+2,'图件数量不对')
        for pdf in pdfs:
            for ext in ['.pdf','.png','.tsv','.ggplot.rds']:
                path = pdf.with_suffix(ext)
                io.require(path.exists() and path.stat().st_size>0,f'缺少图件：{path}')
        for gene in targets.get(item.material,[]):
            row = pd.to_numeric(counts.loc[gene])
            io.require(row.max()>0,f'问题基因仍全零：{gene}')
            target_counts.extend({'material':item.material,'gene_id':gene,'sample_id':sample,'stringtie_count':int(value)} for sample,value in row.items())
            if gene=='Zm00001eb056980':
                io.require(row['B73heat0_5h_1']==10893 and row['B73heat0_5h_2']==12434,'ZmHSF17与用户提供的StringTie计数不同')
        qc.append({'material':item.material,'job_id':item.job_id,'genes':len(counts),'samples':len(counts.columns),
                   'contrasts':len(contrasts),'plot_sets':len(pdfs),'zero_count_tpm_ge1':0,'status':'PASS'})
    io.write(pd.DataFrame(qc),out/'04_material_validation.tsv')
    io.write(pd.DataFrame(target_counts),out/'04_target_gene_counts.tsv')
    io.write(pd.DataFrame(target_stats),out/'04_target_gene_statistics.tsv')
    (out/'04_REAL_PASS').write_text('完整性、阈值、图件和关键基因验收通过\n')
    log(f'正式批次验收通过：{len(qc)}材料')
    return args

if __name__=='__main__':
    io.run(main)
