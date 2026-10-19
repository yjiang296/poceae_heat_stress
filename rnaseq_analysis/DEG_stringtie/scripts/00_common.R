### 脚本名称：00_common.R
### 功能：统一R脚本参数、日志和带错误处理的读写。
### 输入文件：--input 上一步结果目录；输出文件：--outdir 下日志与分析结果。
### 预期运行路径：流程根目录；环境激活：conda activate rnaseq_de。
### 使用示例：由 Rscript scripts/03_deseq2.R --input 上游目录 --outdir 结果目录调用。
suppressPackageStartupMessages(library(optparse))
set.seed(42)
opts <- list(make_option('--input',type='character'), make_option('--outdir',type='character'),
             make_option('--min-count',type='integer',default=10), make_option('--min-samples',type='integer',default=2),
             make_option('--fdr',type='double',default=0.05), make_option('--lfc',type='double',default=1))
args <- parse_args(OptionParser(option_list=opts), convert_hyphens_to_underscores=TRUE)
if (is.null(args$input) || is.null(args$outdir)) stop('必须指定 --input 和 --outdir')
dir.create(args$outdir, recursive=TRUE, showWarnings=FALSE)
log_line <- function(x) {
  cat(x,'\n'); cat(x,'\n',file=file.path(args$outdir,'run.log'),append=TRUE)
}
log_line(paste('运行目录：',getwd(),'；R.home：',R.home(),'；Conda：',Sys.getenv('CONDA_DEFAULT_ENV')))
read_tsv <- function(path) tryCatch(read.delim(path,check.names=FALSE,stringsAsFactors=FALSE),
                                   error=function(e) stop('读取失败：',path,'；',e$message))
write_tsv <- function(x,name) tryCatch(write.table(x,file.path(args$outdir,name),sep='\t',quote=FALSE,row.names=FALSE),
                                      error=function(e) stop('写入失败：',name,'；',e$message))
save_rds <- function(x,name) tryCatch(saveRDS(x,file.path(args$outdir,name)),error=function(e) stop('保存失败：',e$message))
versions <- function() tryCatch(writeLines(capture.output(sessionInfo()),file.path(args$outdir,'session_info.txt')),
                                error=function(e) stop('版本写入失败：',e$message))
