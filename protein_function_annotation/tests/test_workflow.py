### 测试名称：test_workflow.py；功能：输入、续跑、完整性与原 ID 恢复的回归测试；不是生产入口。
### 环境：conda activate base；示例：python -m unittest discover -s tests -v；使用临时目录，无原始数据写入。
import csv,io,importlib.util,json,pathlib,subprocess,sys,tempfile,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from common import dump,finish,sha,tsv
from validators import ipr,egg
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.d=pathlib.Path(self.tmp.name)
  self.cfg=json.loads((ROOT/'config.example.json').read_text());self.cfg.update(python=sys.executable,chunk_size=1)
  self.config=self.d/'config.json';dump(self.config,self.cfg)
 def tearDown(self):self.tmp.cleanup()
 def runstep(self,step,source,dest,ok=True):
  p=subprocess.run([sys.executable,str(ROOT/'scripts'/step),'--input',str(source),'--outdir',str(dest),'--config',str(self.config),'--sample','Test'],capture_output=True,text=True)
  self.assertEqual(p.returncode==0,ok,p.stdout+p.stderr);return p
 def prepare(self):
  f=self.d/'input.fa';f.write_text('>cds:a details\nMAEK*\n>b\nMAEK\n')
  out=self.d/'sample'/'01';self.runstep('01_prepare.py',f,out);return f,out
 def test_preserve_identical_proteins_and_resume(self):
  f,out=self.prepare();m=json.loads((out/'01_complete.json').read_text());self.assertEqual(m['ready'],2)
  rows=list(csv.DictReader(io.StringIO(pathlib.Path(m['mapping']).read_text()),delimiter='\t'));self.assertEqual(rows[0]['original_id'],'cds:a')
  self.assertEqual(rows[0]['normalization'],'removed_terminal_stop');self.assertEqual(len(m['chunks']),2)
  p=self.runstep('01_prepare.py',f,out);self.assertIn('跳过',p.stdout)
  f.write_text('>a\nMALK\n');self.runstep('01_prepare.py',f,out,False)
 def test_duplicate_invalid_and_quarantine(self):
  f=self.d/'bad.fa';f.write_text('>a\nMAEK\n>a\nMALK\n');self.runstep('01_prepare.py',f,self.d/'dup',False)
  f.write_text('>a\nMAEK\n>b\nMA*K\n');self.runstep('01_prepare.py',f,self.d/'fail',False)
  self.cfg['invalid_policy']='quarantine';dump(self.config,self.cfg);self.runstep('01_prepare.py',f,self.d/'quarantine')
  self.assertEqual(json.loads((self.d/'quarantine/01_complete.json').read_text())['quarantined'],1)
 def test_no_match_vs_missing(self):
  fa=self.d/'q.fa';fa.write_text('>q1\nMAEK\n>q2\nMAEK\n')
  p=self.d/'i.json';dump(p,{'results':[{'xref':[{'id':'q1'},{'id':'q2'}],'sequence':'MAEK','matches':[]}]})
  self.assertEqual(len(ipr(p,fa)),2)
  dump(p,{'results':[{'xref':[{'id':'q1'}],'sequence':'MAEK','matches':[]}]})
  with self.assertRaises(ValueError):ipr(p,fa)
 def test_eggnog_header_footer_and_unknown(self):
  fa=self.d/'q.fa';fa.write_text('>q1\nMAEK\n');prefix=self.d/'a.dot'
  ann=pathlib.Path(str(prefix)+'.emapper.annotations');seed=pathlib.Path(str(prefix)+'.emapper.seed_orthologs')
  ann.write_text('#query\tseed_ortholog\n## Total time (seconds): 1\n');seed.write_text('#qseqid\tsseqid\n## Total time (seconds): 1\n')
  self.assertEqual(egg(prefix,fa),({},{}))
  seed.write_text('#qseqid\tsseqid\nbad\ts\n## Total time (seconds): 1\n')
  with self.assertRaises(ValueError):egg(prefix,fa)
  seed.write_text('#qseqid\tsseqid\n')
  with self.assertRaises(ValueError):egg(prefix,fa)
 def test_merge_restores_ids_and_no_hit(self):
  f,prep=self.prepare();meta=json.loads((prep/'01_complete.json').read_text());root=prep.parent
  two=root/'02';two.mkdir();jp=two/'ipr.json'
  dump(jp,{'results':[{'xref':[{'id':'Test__p000000001'},{'id':'Test__p000000002'}],'sequence':'MAEK','matches':[]}]})
  finish(two/'02_complete.json','fixture',[jp],chunks=[{'json':str(jp),'fasta':meta['proteins']}])
  three=root/'03';three.mkdir();prefix=three/'egg';ann=three/'egg.emapper.annotations';seed=three/'egg.emapper.seed_orthologs'
  ann.write_text('#query\tseed_ortholog\tDescription\nTest__p000000001\ts1\tExample\n## Total time (seconds): 1\n')
  seed.write_text('#qseqid\tsseqid\tqcov\tscov\nTest__p000000001\ts1\t90\t85\n## Total time (seconds): 1\n')
  finish(three/'03_complete.json','fixture',[ann,seed],prefix=str(prefix))
  self.runstep('04_merge.py',root,root/'04');rows=list(csv.DictReader(io.StringIO((root/'04/04_Test_members.tsv').read_text()),delimiter='\t'))
  self.assertEqual([r['original_id'] for r in rows],['cds:a','b']);self.assertEqual(rows[1]['eggnog_status'],'no_reported_seed')
  self.assertEqual(rows[0]['seed_both_coverage_ge80'],'1')
 def test_dryrun_and_malformed_table(self):
  f,_=self.prepare();sheet=self.d/'samples.txt';sheet.write_text('### comment\nA\t'+str(f)+'\nB\t'+str(f)+'\n')
  out=self.d/'batch';p=subprocess.run(['bash',str(ROOT/'submit_annotation.sh'),str(sheet),str(out),'--config',str(self.config),'--no-submit'],capture_output=True,text=True)
  self.assertEqual(p.returncode,0,p.stderr);jobs=list(out.glob('*/jobs/*.slurm'));self.assertEqual(len(jobs),2)
  for job in jobs:
   subprocess.run(['bash','-n',str(job)],check=True);self.assertNotIn('--mem',job.read_text());self.assertIn('#SBATCH -n 8',job.read_text())
  sheet.write_text('A '+str(f)+'\nB '+str(f)+' unexpected\n');bad=self.d/'badbatch'
  p=subprocess.run(['bash',str(ROOT/'submit_annotation.sh'),str(sheet),str(bad),'--config',str(self.config)],capture_output=True,text=True)
  self.assertNotEqual(p.returncode,0);self.assertFalse(list(bad.glob('*/jobs/*.slurm')))
if __name__=='__main__':unittest.main()
