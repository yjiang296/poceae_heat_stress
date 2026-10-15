### 脚本名称：01_prepare.py；功能：检查蛋白、建立原 ID 对照、规范末端终止符并分块。
### 输入：--input 蛋白 FASTA 或 gz；输出：--outdir 下输入表、蛋白及分块；预期路径：任意。
### 环境：conda activate base；示例：python 01_prepare.py --input a.fa --outdir x/01 --config config.json --sample A
### 参数：分块大小和异常字符策略来自配置；不选最长转录本，不删除相同序列，不修改原文件。
from common import *

def run(a):
 if not a.input.is_file():raise ValueError('蛋白 FASTA 不存在')
 if a.cfg['invalid_policy'] not in ('fail','quarantine'):raise ValueError('invalid_policy 必须为 fail 或 quarantine')
 size=int(a.cfg['chunk_size'])
 if size<1:raise ValueError('chunk_size 必须为正整数')
 sig=signature(a,[a.input]);marker=a.outdir/'01_complete.json'
 if cached(marker,sig):log(a,'01 已完成且校验一致，跳过');return
 rows=[];ready=[];seen=set();letters=set()
 for i,(header,seq) in enumerate(fasta(a.input),1):
  original=header.split()[0]
  if original in seen:raise ValueError('重复原始蛋白 ID：'+original)
  seen.add(original);letters.update(seq);q=a.sample+'__p'+str(i).zfill(9)
  norm=seq[:-1] if seq.endswith('*') else seq
  invalid=''.join(sorted(set(norm)-set('ACDEFGHIKLMNPQRSTVWYX')))
  status='ready' if norm and not invalid else 'quarantined'
  rows.append(dict(sample=a.sample,query_id=q,original_id=original,original_header=header,length=len(norm),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),normalization='removed_terminal_stop' if seq!=norm else 'none',input_status=status,invalid_characters=invalid or ('empty_sequence' if not norm else '')))
  if status=='ready':ready.append((q,norm))
 if not rows:raise ValueError('空 FASTA')
 if letters<=set('ACGTN*'):raise ValueError('输入全部仅含核酸字母，请确认不是基因组/CDS FASTA')
 mapping=a.outdir/('01_'+a.sample+'_id_map.tsv');tsv(mapping,rows)
 bad=len(rows)-len(ready)
 if bad and a.cfg['invalid_policy']=='fail':raise ValueError(f'{bad} 条蛋白含非法字符/为空；见 {mapping}。修复输入或显式使用 quarantine')
 if not ready:raise ValueError('没有可注释蛋白')
 aa=a.outdir/('01_'+a.sample+'_proteins.faa');chunks=a.outdir/'chunks';chunks.mkdir(exist_ok=True)
 with aa.open('w') as f:
  for q,s in ready:f.write('>'+q+'\n'+s+'\n')
 paths=[]
 for n,start in enumerate(range(0,len(ready),size),1):
  p=chunks/f'01_{a.sample}_{n:05d}.faa'
  with p.open('w') as f:
   for q,s in ready[start:start+size]:f.write('>'+q+'\n'+s+'\n')
  paths.append(p)
 finish(marker,sig,[mapping,aa]+paths,input_path=str(a.input),input_sha256=sha(a.input),members=len(rows),ready=len(ready),quarantined=bad,mapping=str(mapping),proteins=str(aa),chunks=[str(p) for p in paths])
 log(a,f'01 完成：总数 {len(rows)}，进入注释 {len(ready)}，隔离 {bad}')

if __name__=='__main__':entry(run)
