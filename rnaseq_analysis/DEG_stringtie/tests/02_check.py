### 脚本名称：02_check.py
### 功能：验证模拟差异方向、零基线检验、低表达状态与图件配套源数据。
### 输入文件：--input 模拟结果目录；输出文件：--outdir 验证记录。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python tests/02_check.py --input tests/work/results --outdir tests/work/validation。
import sys
from pathlib import Path
import importlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
io = importlib.import_module('00_common')

def main():
    args = io.parser('验证模拟全流程').parse_args()
    out, log = io.start(args)
    root = Path(args.input)
    files = list((root/'03_deseq2').glob('*_all_genes.tsv'))
    io.require(len(files)==4,'比较数量错误')
    for path in files:
        tab = io.read(path).set_index('gene_id')
        io.require(len(tab)==505,'丢失低表达基因')
        for gene,status in [('known_up','UP'),('known_down','DOWN'),('zero_baseline','UP')]:
            io.require(tab.loc[gene,'response']==status,f'方向不符：{path}:{gene}')
        io.require(tab.loc['all_zero','test_status']=='LOW_COUNT_FILTERED','全零状态错误')
        io.require(tab.loc['low_count','test_status']=='LOW_COUNT_FILTERED','低表达状态错误')
    pdfs = list((root/'04_plots').glob('*.pdf'))
    io.require(len(pdfs)==6,'图件数量错误')
    for pdf in pdfs:
        for ext in ['.pdf','.png','.tsv','.ggplot.rds']:
            path = pdf.with_suffix(ext)
            io.require(path.exists() and path.stat().st_size>0,f'图件配套文件缺失：{path}')
    (out/'02_TEST_PASS').write_text('四比较方向、零基线、低表达和六套图件全部通过\n')
    log('模拟端到端验证通过')
    return args

if __name__=='__main__':
    io.run(main)
