### 脚本名称：02_interpro.py；功能：按蛋白块运行本地 InterProScan 并验证完整 JSON。
### 输入：--input 为 01 结果目录；输出：--outdir 中分块 TSV/JSON 及完成清单；预期路径：任意。
### 环境：conda activate <config.json 中 ipr_runtime 路径>；示例：python 02_interpro.py --input x/01 --outdir x/02 --config config.json --sample A
### 参数：沿用已验证七库、关闭远程匹配查询；每块失败留档、重试新目录，成功块按校验和复用。
from common import *
from validators import ipr

def run(a):
 prep=require(a.input,'01');cfg=a.cfg
 env=os.environ.copy();env['PATH']=cfg['ipr_runtime']+'/bin:'+env.get('PATH','');env['JAVA_HOME']=cfg['ipr_runtime']
 install=Path(cfg['interproscan']).parent
 identity=tool_identity([cfg['interproscan'],Path(cfg['ipr_runtime'])/'bin/java',install/'interproscan.properties'])
 jars=list(install.glob('interproscan-*.jar'));identity+=tool_identity(jars)
 sig=signature(a,[a.input/'01_complete.json'],identity);marker=a.outdir/'02_complete.json'
 if cached(marker,sig):log(a,'02 已完成且校验一致，跳过');return
 execute(a,[cfg['interproscan'],'-version'],env=env,name='version')
 execute(a,[Path(cfg['ipr_runtime'])/'bin/java','-version'],env=env,name='version')
 completed=[];outputs=[]
 for source in map(Path,prep['chunks']):
  folder=a.outdir/'chunks'/source.stem;folder.mkdir(parents=True,exist_ok=True)
  cm=folder/'02_chunk_complete.json';cs=hashlib.sha256((sig+sha(source)).encode()).hexdigest();old=cached(cm,cs)
  if old:completed.append(old);outputs.append(cm);log(a,'复用 '+source.stem);continue
  trial=attempt(folder);prefix=trial/('02_'+a.sample);temp=trial/'tmp'
  execute(a,[cfg['interproscan'],'-i',source,'-b',prefix,'-f','tsv,json','-goterms','-dp','-cpu',a.threads,'-appl',','.join(cfg['ipr_applications']),'-T',temp],env=env,cwd=trial,name=source.stem)
  jp=Path(str(prefix)+'.json');tp=Path(str(prefix)+'.tsv')
  observed=ipr(jp,source)
  if not tp.is_file():raise ValueError('InterPro 缺少 TSV；无命中时允许空文件')
  finish(cm,cs,[jp,tp],json=str(jp),tsv=str(tp),fasta=str(source),members=len(observed))
  completed.append(json.loads(cm.read_text()));outputs.append(cm)
  log(a,'完成 '+source.stem+'，蛋白数 '+str(len(observed)))
 if sum(x['members'] for x in completed)!=prep['ready']:raise ValueError('InterPro 分块总数不守恒')
 for d in completed:outputs.extend(Path(p) for p in d['outputs'])
 finish(marker,sig,outputs,members=prep['ready'],chunks=completed,software_identity=identity)
 log(a,'02 全部分块完成')

if __name__=='__main__':entry(run)
