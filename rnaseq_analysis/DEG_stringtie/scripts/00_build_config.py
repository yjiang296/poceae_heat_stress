### 脚本名称：00_build_config.py
### 功能：从旧比较配置和定量样本表生成独立的逐材料配置；按真实BAM映射样本别名。
### 输入文件：--input 旧比较目录；--quant-root 原StringTie流程目录。
### 输出文件：每材料样本、比较和定量来源表，以及 sample_info.txt。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：python scripts/00_build_config.py --input 旧目录 --outdir config/materials --quant-root 定量目录。
import re
from pathlib import Path
import pandas as pd
import importlib
io = importlib.import_module('00_common')

def main():
    p = io.parser('根据已有比较生成StringTie流程配置')
    p.add_argument('--quant-root', required=True)
    args = p.parse_args()
    out, log = io.start(args)
    quant = Path(args.quant_root).resolve()
    catalog = {}
    for path in sorted(quant.glob('sample_info*.txt')):
        rows = [line.split() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
        for values in rows[1:]:
            item = dict(zip(rows[0], values))
            sample = item['sample']
            key = str(quant/'results/samples'/sample/'02_hisat2'/f'{sample}.sorted.bam')
            io.require(key not in catalog or catalog[key] == item, f'原样本配置冲突：{sample}')
            catalog[key] = item
    materials, comparisons, sources, controls = {}, {}, {}, {}
    paths = sorted(Path(args.input).glob('comparison.*.tsv'))
    paths = [x for x in paths if x.name != 'comparison.example.tsv']
    io.require(paths, '没有正式比较配置')
    for path in paths:
        tab = io.read(path)
        io.require(set(tab.columns) == {'sample_id','condition','bam','annotation_gtf','reference_name'}, f'表头错误：{path}')
        io.require(not tab.sample_id.duplicated().any(), f'比较内样本重复：{path}')
        io.require(set(tab.condition) == {'control','treatment'}, f'必须有对照与处理：{path}')
        io.require(tab.groupby('condition').size().min() >= 2, f'生物学重复少于2：{path}')
        refs = set()
        cid = path.name.removeprefix('comparison.').removesuffix('.tsv')
        io.require(re.fullmatch(r'[A-Za-z0-9_.-]+', cid), '比较名含不安全字符')
        for row in tab.itertuples(index=False):
            io.require(row.bam in catalog, f'BAM无法对应StringTie样本：{row.bam}')
            original = catalog[row.bam]
            ref, sample = original['reference_name'], original['sample']
            refs.add(ref)
            io.require(original['annotation_gtf'] == row.annotation_gtf, f'注释不一致：{sample}')
            condition = 'CK' if row.condition == 'control' else cid
            item = {'sample_id':sample, 'condition':condition, 'bam':row.bam, 'annotation_gtf':row.annotation_gtf}
            materials.setdefault(ref,{})
            io.require(sample not in materials[ref] or materials[ref][sample] == item, f'样本分组冲突：{sample}')
            materials[ref][sample] = item
        io.require(len(refs) == 1, f'跨材料比较不支持：{path}')
        ref = next(iter(refs))
        control_bams = set(tab.loc[tab.condition=='control','bam'])
        io.require(ref not in controls or controls[ref]==control_bams, f'共享CK不是同一批样本：{ref}')
        controls[ref] = control_bams
        comparisons.setdefault(ref, []).append({'contrast_id':cid,'numerator':cid,'denominator':'CK'})
        sources.setdefault(ref, []).append({'path':str(path.resolve()),'sha256':io.digest(path)})
    manifest = []
    for ref, samples in sorted(materials.items()):
        dest = out/ref
        io.require(not dest.exists(), f'拒绝覆盖现有配置：{dest}')
        dest.mkdir()
        rows = pd.DataFrame(samples.values()).sort_values(['condition','sample_id'])
        io.write(rows, dest/'00_samples.tsv')
        io.write(pd.DataFrame(comparisons[ref]), dest/'00_contrasts.tsv')
        base = quant/'results/summary'
        source = {'material':ref,'counts':str(base/'05_counts'/ref/f'{ref}.gene_count_matrix.tsv'),
                  'tpm':str(base/'06_abundance'/ref/f'{ref}.gene_tpm_matrix.tsv'),
                  'quant_root':str(quant),'method':'STRINGTIE_PREPDE_FORMULA'}
        io.write(pd.DataFrame([source]), dest/'00_sources.tsv')
        io.write(pd.DataFrame(sources[ref]), dest/'00_original_comparisons.tsv')
        manifest.append({'material':ref,'config_dir':str(dest.resolve())})
        log(f'{ref}：{len(rows)}个样本，{len(comparisons[ref])}个比较')
    io.write(pd.DataFrame(manifest), out/'sample_info.txt')
    return args

if __name__ == '__main__':
    io.run(main)
