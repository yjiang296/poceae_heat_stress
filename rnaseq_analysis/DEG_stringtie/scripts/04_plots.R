### 脚本名称：04_plots.R
### 功能：绘制全基因PCA、样本相关性和各比较火山图；保存矢量PDF、300dpi PNG、绘图TSV和对象。
### 输入文件：--input 步骤03目录；输出文件：以04及材料名称开头的图件和源数据。
### 预期运行路径：流程根目录；环境激活：conda activate rnaseq_de。
### 使用示例：Rscript scripts/04_plots.R --input results/B73/03_deseq2 --outdir results/B73/04_plots。
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(script),'00_common.R'))
suppressPackageStartupMessages(library(ggplot2))
main <- function() {
  if (!file.exists(file.path(args$input,'03_SUCCESS'))) stop('缺少03结果，请先运行03或all')
  meta <- read_tsv(file.path(args$input,'03_samples.tsv'))
  material <- unique(meta$material); prefix <- paste0('04_',material)
  vst <- read_tsv(file.path(args$input,paste0('03_',material,'_vst.tsv')))
  mat <- as.matrix(vst[,-1]); rownames(mat) <- vst$gene_id
  palette <- c('#000000','#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7')
  theme_set(theme_classic(base_family='Arial',base_size=10))
  emit <- function(p,dat,suffix) {
    name <- paste0(prefix,'_',suffix)
    write_tsv(dat,paste0(name,'.tsv')); save_rds(p,paste0(name,'.ggplot.rds'))
    ggsave(file.path(args$outdir,paste0(name,'.pdf')),p,device=cairo_pdf,width=7,height=5)
    ggsave(file.path(args$outdir,paste0(name,'.png')),p,width=7,height=5,dpi=300)
  }
  ### PCA使用VST后方差最高的至多500个基因，参数参照DESeq2的常用展示方式。
  vv <- apply(mat,1,var); chosen <- order(vv,decreasing=TRUE)[seq_len(min(500,nrow(mat)))]
  pc <- prcomp(t(mat[chosen,,drop=FALSE]))
  frac <- pc$sdev^2/sum(pc$sdev^2)
  pts <- data.frame(sample_id=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2])
  pts$condition <- meta$condition[match(pts$sample_id,meta$sample_id)]
  pts$condition <- sub('^.*_CK_and_','',pts$condition)
  p <- ggplot(pts,aes(PC1,PC2,color=condition))+geom_point(size=3)+geom_text(aes(label=sample_id),vjust=-0.7,size=2.5)+
    scale_color_manual(values=palette)+labs(x=sprintf('PC1 (%.1f%%)',100*frac[1]),y=sprintf('PC2 (%.1f%%)',100*frac[2]))
  emit(p,pts,'pca')
  corr <- as.data.frame(as.table(cor(mat))); names(corr) <- c('sample_1','sample_2','pearson_r')
  p <- ggplot(corr,aes(sample_1,sample_2,fill=pearson_r))+geom_tile()+coord_equal()+
    scale_fill_gradient(low='#FFFFFF',high='#0072B2',limits=c(0,1))+labs(x=NULL,y=NULL)+
    theme(axis.text.x=element_text(angle=60,hjust=1,size=7),axis.text.y=element_text(size=7))
  emit(p,corr,'sample_correlation')
  contrasts <- read_tsv(file.path(args$input,'03_contrasts.tsv'))
  for (cid in contrasts$contrast_id) {
    full <- read_tsv(file.path(args$input,paste0('03_',material,'_',cid,'_all_genes.tsv')))
    valid <- !is.na(full$padj) & is.finite(full$log2FoldChange)
    ### 极小padj仅为绘图截断到1e-300，完整统计表保持原值；排除不可评价点。
    dat <- data.frame(gene_id=full$gene_id[valid],log2FoldChange=full$log2FoldChange[valid],
                      neglog10_padj_capped300=-log10(pmax(full$padj[valid],1e-300)),response=full$response[valid])
    p <- ggplot(dat,aes(log2FoldChange,neglog10_padj_capped300,color=response))+geom_point(size=0.5,alpha=0.5)+
      scale_color_manual(values=c(UP='#D55E00',DOWN='#0072B2',NO_DETECTED_RESPONSE='#999999'))+
      geom_vline(xintercept=c(-args$lfc,args$lfc),linetype=2,color='#777777')+
      geom_hline(yintercept=-log10(args$fdr),linetype=2,color='#777777')+
      labs(x='log2 fold change (treatment / CK)',y='-log10 adjusted P (capped at 300)',title=cid)
    emit(p,dat,paste0(cid,'volcano'))
  }
  versions(); writeLines('绘图完成',file.path(args$outdir,'04_SUCCESS'))
  log_line(paste(material,'图件完成'))
}
tryCatch(main(),error=function(e) {message('错误：',conditionMessage(e)); quit(status=1)})
