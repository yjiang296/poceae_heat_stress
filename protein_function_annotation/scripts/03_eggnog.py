### 脚本名称：03_eggnog.py；功能：每材料整体运行 eggNOG，记录数据库与种子覆盖度。
### 输入：--input 为 01 结果目录；输出：--outdir 中原始注释及完成清单；预期路径：任意。
### 环境：conda activate <config.json 中 eggnog_env 路径>；示例：python 03_eggnog.py --input x/01 --outdir x/03 --config config.json --sample A
### 参数：沿用已验证版本；自动分类范围、非电子 GO、禁用 Pfam 重比对，保留原始命中。
from common import *
from validators import egg

def run(a):
 prep=require(a.input,'01');cfg=a.cfg
 env=os.environ.copy();env['PATH']=cfg['eggnog_env']+'/bin:'+env.get('PATH','');env['EGGNOG_DATA_DIR']=cfg['eggnog_data']
 exe=Path(cfg['eggnog_env'])/'bin/emapper.py';db=Path(cfg['eggnog_data'])
 identity=tool_identity([exe,Path(cfg['eggnog_env'])/'bin/diamond']+[db/n for n in ('eggnog.db','eggnog.taxa.db','eggnog.taxa.db.traverse.pkl','eggnog_proteins.dmnd')])
 sig=signature(a,[a.input/'01_complete.json'],identity);marker=a.outdir/'03_complete.json'
 if cached(marker,sig):log(a,'03 已完成且校验一致，跳过');return
 version=subprocess.check_output([str(exe),'--version'],env=env,stderr=subprocess.STDOUT,text=True)
 with (a.outdir/'version.log').open('a') as f:f.write(version)
 if 'Expected eggNOG DB version: 5.0.2' not in version or 'Installed eggNOG DB version: 5.0.2' not in version:raise ValueError('数据库版本未通过已验证的 5.0.2 检查；见 version.log')
 trial=attempt(a.outdir);prefix=trial/('03_'+a.sample);temp=trial/'tmp';temp.mkdir()
 execute(a,[exe,'-i',prep['proteins'],'--itype','proteins','-m','diamond','--cpu',a.threads,'--data_dir',db,'--tax_scope',cfg['tax_scope'],'--target_orthologs','all','--go_evidence',cfg['go_evidence'],'--pfam_realign','none','--output',prefix.name,'--output_dir',trial,'--temp_dir',temp,'--block_size',cfg['diamond_block_size'],'--sensmode',cfg['sensitivity'],'--dmnd_iterate','yes'],env=env,cwd=trial,name='03_'+a.sample)
 ann,seeds=egg(prefix,Path(prep['proteins']))
 files=list(trial.glob(prefix.name+'.emapper.*'))
 finish(marker,sig,files,members=prep['ready'],annotated=len(ann),seed_hits=len(seeds),prefix=str(prefix),software_database_identity=identity,version=version)
 log(a,f'03 完成：输入 {prep["ready"]}，种子 {len(seeds)}，注释记录 {len(ann)}')

if __name__=='__main__':entry(run)
