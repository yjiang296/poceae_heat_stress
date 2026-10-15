### 脚本名称：04_merge.py；功能：恢复原蛋白 ID 并整合两套注释，保留无命中及隔离记录。
### 输入：--input 为单材料根目录，含 01/02/03；输出：--outdir 下成员表、结构域表和计数；预期路径：任意。
### 环境：conda activate base；示例：python 04_merge.py --input x --outdir x/04 --config config.json --sample A
### 参数：双向覆盖度 80% 仅作描述性筛选，不是功能确认阈值；无随机过程。
from common import *
from validators import ipr,egg

def run(a):
 prep=require(a.input/'01','01');ip=require(a.input/'02','02');eg=require(a.input/'03','03')
 markers=[a.input/s/(s+'_complete.json') for s in ('01','02','03')]
 sig=signature(a,markers);marker=a.outdir/'04_complete.json'
 if cached(marker,sig):log(a,'04 已完成且校验一致，跳过');return
 observed={}
 for chunk in ip['chunks']:
  data=ipr(Path(chunk['json']),Path(chunk['fasta']))
  if set(data)&set(observed):raise ValueError('InterPro 跨块重复 ID')
  observed.update(data)
 annotations,seeds=egg(Path(eg['prefix']),Path(prep['proteins']))
 rows=table(prep['mapping']);expected={r['query_id'] for r in rows if r['input_status']=='ready'}
 if set(observed)!=expected:raise ValueError('InterPro 总集合不匹配')
 domains=[];merged=[]
 for r in rows:
  q=r['query_id'];hits=[];gos=set();iprs=set()
  for match in observed.get(q,[]):
   s=match['signature'];lib=s.get('signatureLibraryRelease',{});entry=s.get('entry') or {}
   if entry.get('accession'):iprs.add(entry['accession'])
   for g in entry.get('goXRefs',[]):gos.add(g['id'])
   for loc in match.get('locations',[]):
    h=dict(sample=a.sample,original_id=r['original_id'],query_id=q,database=lib.get('library',''),version=lib.get('version',''),signature=s['accession'],description=s.get('description') or s.get('name') or '',interpro=entry.get('accession',''),entry_type=entry.get('type',''),start=loc['start'],end=loc['end'],score=loc.get('score',match.get('score','')),evalue=loc.get('evalue',match.get('evalue','')))
    hits.append(h);domains.append(h)
  row=dict(r);pf=sorted([h for h in hits if h['database'].upper()=='PFAM'],key=lambda h:(h['start'],h['end'],h['signature']))
  row.update(ipr_status='not_run_quarantined' if q not in observed else ('hit' if observed[q] else 'no_match'),pfam_architecture='|'.join(h['signature'] for h in pf),interpro_accessions=';'.join(sorted(iprs)),interpro_go=';'.join(sorted(gos)))
  row['eggnog_status']='not_run_quarantined' if r['input_status']!='ready' else ('annotation_row' if q in annotations else ('seed_without_annotation' if q in seeds else 'no_reported_seed'))
  for col in ('sseqid','evalue','bitscore','pident','qcov','scov'):row['eggnog_seed_'+col]=seeds.get(q,{}).get(col,'')
  for col in ('Description','Preferred_name','COG_category','max_annot_lvl','eggNOG_OGs','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module','KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs'):row['eggnog_'+col]=annotations.get(q,{}).get(col,'')
  row['seed_both_coverage_ge80']=int(q in seeds and min(float(seeds[q]['qcov']),float(seeds[q]['scov']))>=80);merged.append(row)
 mp=a.outdir/('04_'+a.sample+'_members.tsv');dp=a.outdir/('04_'+a.sample+'_domains.tsv');sp=a.outdir/('04_'+a.sample+'_summary.json')
 tsv(mp,merged);tsv(dp,domains,['sample','original_id','query_id','database','version','signature','description','interpro','entry_type','start','end','score','evalue'])
 summary=dict(sample=a.sample,input_members=len(rows),annotation_ready=prep['ready'],quarantined=prep['quarantined'],interpro_hit=sum(r['ipr_status']=='hit' for r in merged),eggnog_seed=len(seeds),eggnog_annotation=len(annotations),neither_ipr_nor_seed=sum(r['ipr_status']=='no_match' and r['eggnog_status']=='no_reported_seed' for r in merged))
 dump(sp,summary);finish(marker,sig,[mp,dp,sp],**summary);log(a,'04 完成 '+json.dumps(summary,ensure_ascii=False))

if __name__=='__main__':entry(run)
