### 脚本名称：common.py；功能：公共日志、输入输出、校验及续跑；不是独立步骤入口。
### 输入：各步骤传入的路径；输出：日志、版本与完成清单；预期路径：任意。
### 环境：conda activate base；使用示例：由 01 至 04 步骤导入；标准库实现无随机抽样。
import argparse, csv, datetime, gzip, hashlib, json, os, random, subprocess, sys
from pathlib import Path
random.seed(42)

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def dump(path,data):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n');temp.replace(p)

def tsv(path,rows,fields=None):
 with Path(path).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def table(path):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def fasta(path):
 opener=gzip.open if str(path).endswith('.gz') else open
 with opener(path,'rt') as f:
  name=None;seq=[]
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq)
    name=line[1:].strip();seq=[]
    if not name:raise ValueError('FASTA 标题为空')
   elif line.strip():
    if name is None:raise ValueError('FASTA 序列前缺少标题')
    seq.append(''.join(line.split()).upper())
  if name is not None:yield name,''.join(seq)

def args():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
 p.add_argument('--config',type=Path,required=True);p.add_argument('--sample',required=True);p.add_argument('--threads',type=int,default=2)
 a=p.parse_args();a.input=a.input.resolve();a.outdir=a.outdir.resolve();a.config=a.config.resolve()
 if a.threads<1:raise ValueError('线程数必须为正整数')
 a.outdir.mkdir(parents=True,exist_ok=True);a.cfg=json.loads(a.config.read_text());return a

def log(a,message):
 text=datetime.datetime.now().isoformat(timespec='seconds')+' '+str(message)
 print(text,flush=True)
 with (a.outdir/'run.log').open('a') as f:f.write(text+'\n')

def execute(a,cmd,env=None,cwd=None,name='tool'):
 log(a,'运行 '+json.dumps([str(x) for x in cmd],ensure_ascii=False))
 with (a.outdir/(name+'.log')).open('a') as f:
  subprocess.run([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT,env=env,cwd=cwd,check=True)

def signature(a,inputs,extra=None):
 sources={str(p):sha(p) for p in inputs}
 code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
 data=dict(sample=a.sample,threads=a.threads,config=a.cfg,inputs=sources,code=code,extra=extra)
 return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()

def cached(path,sig):
 if not path.exists():return None
 d=json.loads(path.read_text())
 if d['signature']!=sig:raise ValueError('已有结果的输入/配置/代码已改变，请使用新的输出目录：'+str(path))
 for p,h in d['outputs'].items():
  if not Path(p).is_file() or sha(p)!=h:raise ValueError('已有输出缺失或被改动：'+p)
 return d

def finish(path,sig,outputs,**extra):
 dump(path,dict(signature=sig,outputs={str(p.resolve()):sha(p) for p in outputs},**extra))

def require(folder,stage):
 p=folder/(stage+'_complete.json')
 if not p.is_file():raise ValueError('缺少上游完成清单 '+str(p)+'；请先运行 '+stage+' 或 all')
 d=json.loads(p.read_text())
 cached(p,d['signature']);return d

def attempt(folder):
 p=folder/('attempt_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f'));p.mkdir(parents=True);return p

def tool_identity(paths):
 result=[]
 for p in map(Path,paths):
  if not p.exists():raise ValueError('缺少软件或数据库：'+str(p))
  st=p.stat();result.append(dict(path=str(p.resolve()),size=st.st_size,mtime_ns=st.st_mtime_ns))
 return result

def entry(fn):
 a=None
 try:
  a=args();log(a,dict(cwd=os.getcwd(),python=sys.executable,conda=os.environ.get('CONDA_DEFAULT_ENV','未检测到 conda 环境')))
  fn(a)
  with (a.outdir/'version.log').open('a') as f:f.write(sys.version+'\n')
  with (a.outdir/'requirements.txt').open('w') as f:subprocess.run([sys.executable,'-m','pip','freeze'],stdout=f,stderr=subprocess.STDOUT,check=False)
 except Exception as e:
  print('错误：'+str(e),file=sys.stderr)
  if a:log(a,'失败：'+str(e))
  sys.exit(1)
