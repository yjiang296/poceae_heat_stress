### 脚本名称：03_deseq2.R
### 功能：对单材料全部条件联合拟合DESeq2，输出既定比较、整体条件LRT及完整不可评价状态。
### 输入文件：--input 步骤02目录；输出文件：逐比较全基因与DEG表、LRT、标准化counts、VST及模型。
### 预期运行路径：流程根目录；环境激活：conda activate rnaseq_de。
### 使用示例：Rscript scripts/03_deseq2.R --input results/B73/02_matrices --outdir results/B73/03_deseq2。
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(script),'00_common.R'))
suppressPackageStartupMessages(library(DESeq2))
main <- function() {
  if (!file.exists(file.path(args$input,'02_SUCCESS'))) stop('缺少02结果，请先运行02或all')
  raw <- read_tsv(file.path(args$input,'02_gene_counts.tsv'))
  meta <- read_tsv(file.path(args$input,'02_samples.tsv'))
  contrasts <- read_tsv(file.path(args$input,'02_contrasts.tsv'))
  material <- unique(meta$material); stopifnot(length(material)==1)
  prefix <- paste0('03_',material)
  stopifnot(!anyDuplicated(raw$gene_id), !anyDuplicated(meta$sample_id))
  mat <- as.matrix(raw[,meta$sample_id,drop=FALSE]); rownames(mat) <- raw$gene_id
  stopifnot(is.numeric(mat), all(is.finite(mat)), all(mat>=0), all(mat==floor(mat)))
  stopifnot(args$min_count>=0, args$min_samples>=2, args$fdr>0, args$fdr<1, args$lfc>=0)
  storage.mode(mat) <- 'integer'
  meta$condition <- relevel(factor(meta$condition),ref='CK'); rownames(meta) <- meta$sample_id
  stopifnot(all(table(meta$condition)>=2), all(colSums(mat)>0))
  ### 全材料条件联合预过滤与拟合，参数沿用原流程，不将生物学重复编号作为配对项。
  keep <- rowSums(mat>=args$min_count)>=args$min_samples
  if (sum(keep)<10) stop('过滤后基因过少，检查输入')
  dds <- DESeqDataSetFromMatrix(mat[keep,,drop=FALSE],meta,design=~condition)
  dds <- DESeq(dds,quiet=TRUE)
  vsd <- varianceStabilizingTransformation(dds,blind=TRUE)
  summary <- list()
  for (i in seq_len(nrow(contrasts))) {
    cc <- contrasts[i,]
    res <- as.data.frame(results(dds,contrast=c('condition',cc$numerator,cc$denominator),alpha=args$fdr))
    full <- data.frame(gene_id=raw$gene_id,res[match(raw$gene_id,rownames(res)),],row.names=NULL)
    full$passes_prefilter <- keep
    full$test_status <- ifelse(!keep,'LOW_COUNT_FILTERED',ifelse(is.na(full$pvalue),'NO_VALID_PVALUE',
                               ifelse(is.na(full$padj),'NO_VALID_PADJ','EVALUABLE')))
    full$response <- 'NOT_EVALUABLE'
    valid <- !is.na(full$padj) & is.finite(full$log2FoldChange)
    full$response[valid] <- 'NO_DETECTED_RESPONSE'
    full$response[valid & full$padj<args$fdr & full$log2FoldChange>=args$lfc] <- 'UP'
    full$response[valid & full$padj<args$fdr & full$log2FoldChange<=-args$lfc] <- 'DOWN'
    name <- paste0(prefix,'_',cc$contrast_id)
    write_tsv(full,paste0(name,'_all_genes.tsv'))
    write_tsv(full[full$response %in% c('UP','DOWN'),],paste0(name,'_deg.tsv'))
    write_tsv(full[full$response=='UP',],paste0(name,'_up.tsv'))
    write_tsv(full[full$response=='DOWN',],paste0(name,'_down.tsv'))
    summary[[i]] <- data.frame(material=material,contrast_id=cc$contrast_id,numerator=cc$numerator,
                              denominator=cc$denominator,input_genes=nrow(mat),passed_prefilter=sum(keep),
                              evaluable=sum(valid),up=sum(full$response=='UP'),down=sum(full$response=='DOWN'))
  }
  ### LRT检验全部条件是否存在表达差异，不能解释为材料间交互效应。
  lrt <- DESeq(dds,test='LRT',reduced=~1,quiet=TRUE)
  lr <- as.data.frame(results(lrt,alpha=args$fdr))
  lr <- lr[,c('baseMean','stat','pvalue','padj')]
  write_tsv(data.frame(gene_id=raw$gene_id,lr[match(raw$gene_id,rownames(lr)),],passes_prefilter=keep,row.names=NULL),
            paste0(prefix,'_condition_lrt.tsv'))
  normalized <- sweep(mat,2,sizeFactors(dds),'/')
  write_tsv(data.frame(gene_id=rownames(mat),normalized,check.names=FALSE),paste0(prefix,'_normalized_counts.tsv'))
  write_tsv(data.frame(gene_id=rownames(vsd),assay(vsd),check.names=FALSE),paste0(prefix,'_vst.tsv'))
  meta$size_factor <- sizeFactors(dds)
  write_tsv(meta,'03_samples.tsv'); write_tsv(contrasts,'03_contrasts.tsv')
  write_tsv(do.call(rbind,summary),paste0(prefix,'_summary.tsv'))
  write_tsv(data.frame(parameter=c('count_source','design','min_count','min_samples','fdr','abs_log2fc','lfc_type'),
                      value=c('StringTie estimated counts','~condition',args$min_count,args$min_samples,args$fdr,args$lfc,'unshrunk')),
            paste0(prefix,'_parameters.tsv'))
  save_rds(dds,paste0(prefix,'_wald_dds.rds')); save_rds(lrt,paste0(prefix,'_lrt_dds.rds'))
  versions()
  writeLines('模型完成',file.path(args$outdir,'03_SUCCESS'))
  log_line(paste(material,'完成：',nrow(mat),'个输入基因；',sum(keep),'个通过预过滤'))
}
tryCatch(main(),error=function(e) {message('错误：',conditionMessage(e)); quit(status=1)})
