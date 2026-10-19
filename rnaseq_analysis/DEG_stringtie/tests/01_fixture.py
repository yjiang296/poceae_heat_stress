### 脚本名称：01_fixture.py
### 功能：生成带已知响应、零基线和全零基因的模拟输入，以验证完整流程。
### 输入文件：--input 流程根目录；输出文件：--outdir 下模拟定量、配置和样本源文件。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python tests/01_fixture.py --input . --outdir tests/work。
import sys
from pathlib import Path
import json
import numpy as np
import pandas as pd
import importlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
io = importlib.import_module('00_common')

def main():
    args = io.parser('生成模拟测试数据').parse_args()
    out, log = io.start(args); out = out.resolve()
    rng = np.random.default_rng(42)
    quant = out/'quant'; quant.mkdir()
    annotation = out/'fixture.gtf'; annotation.write_text('模拟注释\n')
    conditions = ['CK','heat0_5h','heat1h','heat2h','heat4h']
    samples = [f'{c}_{r}' for c in conditions for r in [1,2]]
    baseline = rng.uniform(30,800,500)
    counts = np.column_stack([rng.negative_binomial(30,30/(30+baseline)) for _ in samples])
    extra = np.array([[100,110,3000,3100,3200,3000,3300,3100,3200,3300],
                      [3000,3200,100,110,90,105,100,95,110,100],
                      [0,0,5000,5100,4000,4100,3000,3100,2000,2100],[0]*10,[1]*10])
    counts = np.vstack([counts,extra]); genes = [f'g{i}' for i in range(500)]+['known_up','known_down','zero_baseline','all_zero','low_count']
    count_file = out/'fixture_counts.tsv'; tpm_file = out/'fixture_tpm.tsv'
    frame = pd.DataFrame(counts,columns=samples); frame.insert(0,'gene_id',genes); io.write(frame,count_file)
    tpm = pd.DataFrame(counts/counts.sum(axis=0)*1e6,columns=samples); tpm.insert(0,'gene_id',genes); io.write(tpm,tpm_file)
    rows = []
    for sample in samples:
        sample_dir = quant/'results/samples'/sample
        (sample_dir/'02_hisat2').mkdir(parents=True); (sample_dir/'03_stringtie').mkdir(); (sample_dir/'01_fastp').mkdir()
        bam = sample_dir/'02_hisat2'/f'{sample}.sorted.bam'; bam.write_text('模拟占位文件\n')
        (sample_dir/'03_stringtie'/f'{sample}.stringtie.gtf').write_text(f'# stringtie {bam} -G {annotation} -e\n# StringTie version fixture\n')
        (sample_dir/'03_stringtie'/f'{sample}.gene_abund.tsv').write_text('模拟占位文件\n')
        (sample_dir/'01_fastp'/f'{sample}.fastp.json').write_text(json.dumps({'summary':{'after_filtering':{'total_bases':15000,'total_reads':100}}}))
        rows.append({'sample_id':sample,'condition':sample.rsplit('_',1)[0],'bam':str(bam),'annotation_gtf':str(annotation)})
    config = out/'config'; config.mkdir()
    io.write(pd.DataFrame(rows),config/'00_samples.tsv')
    io.write(pd.DataFrame([{'contrast_id':c,'numerator':c,'denominator':'CK'} for c in conditions[1:]]),config/'00_contrasts.tsv')
    io.write(pd.DataFrame([{'material':'fixture','counts':str(count_file),'tpm':str(tpm_file),'quant_root':str(quant),'method':'STRINGTIE_PREPDE_FORMULA'}]),config/'00_sources.tsv')
    io.write(pd.DataFrame(columns=['path','sha256']),config/'00_original_comparisons.tsv')
    log('模拟数据生成完成')
    return args

if __name__=='__main__':
    io.run(main)
