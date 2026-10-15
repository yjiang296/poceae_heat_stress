### 脚本名称：validators.py；功能：直接核对注释输出，保留无命中与缺失的区别；公共模块非入口。
### 输入：蛋白与工具输出路径；输出：经校验字典；环境：conda activate base；使用：步骤 02 至 04 导入。
import json
from common import fasta
from pathlib import Path

def ipr(path,proteins):
 expected={h.split()[0]:s for h,s in fasta(proteins)};observed={}
 for r in json.loads(path.read_text()).get('results',[]):
  for ref in r.get('xref',[]):
   q=ref['id']
   if q in observed or q not in expected or r['sequence'].upper()!=expected[q]:raise ValueError('InterPro 序列/ID 不一致：'+q)
   observed[q]=r.get('matches',[])
 if set(observed)!=set(expected):raise ValueError('InterPro 未完整处理输入：'+str(path))
 return observed

def egg_table(path,key):
 fields=None;done=False;rows={}
 for line in path.read_text().splitlines():
  if line.startswith('#'+key+'\t'):fields=line[1:].split('\t')
  elif line.startswith('## Total time (seconds):'):done=True
  elif line and not line.startswith('#'):
   cols=line.split('\t')
   if not fields or len(cols)!=len(fields):raise ValueError('eggNOG 列数/标题异常：'+str(path))
   row=dict(zip(fields,cols));q=row[key]
   if q in rows:raise ValueError('eggNOG 重复 ID：'+q)
   rows[q]=row
 if not fields or not done:raise ValueError('eggNOG 缺少标题或完成尾注：'+str(path))
 return rows

def egg(prefix,proteins):
 ann=egg_table(Path(str(prefix)+'.emapper.annotations'),'query')
 seeds=egg_table(Path(str(prefix)+'.emapper.seed_orthologs'),'qseqid')
 expected={h.split()[0] for h,s in fasta(proteins)}
 if not set(ann)<=set(seeds)<=expected:raise ValueError('eggNOG 出现未知 ID 或注释无对应种子')
 for q,r in ann.items():
  if r['seed_ortholog']!=seeds[q]['sseqid']:raise ValueError('eggNOG 种子与注释不一致')
 return ann,seeds
