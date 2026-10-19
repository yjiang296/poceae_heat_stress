# 二代 RNA-seq 基础流程

本流程用于批量完成双端、非链特异性 RNA-seq 的 fastp 质控、HISAT2 比对、StringTie 已知转录本定量、质控与比对率汇总，以及按参考基因组分别生成表达矩阵。

## 1. 环境检查

```bash
mamba activate rnaseq

fastp --version
hisat2 --version
hisat2_extract_splice_sites.py --help
samtools --version
stringtie --version
multiqc --version
python -c 'import numpy; print(numpy.__version__)'
python -c 'import sys; assert sys.version_info >= (3, 8); print(sys.version)'
```

流程汇总脚本要求 Python 3.8 或更高版本，不支持已经停止维护的 Python 3.5。

如果现有环境仍使用 Python 3.5，并且没有其他项目依赖该环境，建议导出备份后重建同名环境：

```bash
mamba env export -n rnaseq --no-builds > rnaseq_before_python_upgrade.yml

mamba env remove -n rnaseq -y

mamba env create -f environment.yml
```

`environment.yml` 使用 `conda-forge` 和 `bioconda` 安装主要生信软件，通过 pip 安装 MultiQC，以避开部分集群镜像中 MultiQC Conda 依赖缺失的问题。创建完成后检查：

```bash
mamba run -n rnaseq python --version
mamba run -n rnaseq fastp --version
mamba run -n rnaseq hisat2 --version
mamba run -n rnaseq samtools --version
mamba run -n rnaseq stringtie --version
mamba run -n rnaseq multiqc --version
mamba run -n rnaseq bash -c 'command -v hisat2_extract_splice_sites.py'
```

## 2. 样本信息表

编辑项目根目录的 `sample_info.txt`。第一行表头必须保持不变，每个样本一行，共六列：

```text
sample  raw_fq_gz_R1  raw_fq_gz_R2  hisat2_index_prefix  annotation_gtf  reference_name
```

示例：

```text
sample_01	/path/sample_01_R1.fastq.gz	/path/sample_01_R2.fastq.gz	/path/hisat2_index/BTx623	/path/BTx623.annotation.gtf	BTx623
```

注意：

- 推荐使用 Tab 分隔，也支持普通空白符分隔。
- 文件路径中不能包含空格。
- `hisat2_index_prefix` 必须是索引前缀，而不是索引目录。例如索引文件为 `BTx623.1.ht2`，该列应填写 `/path/hisat2_index/BTx623`。
- 流程会检查全部 8 个索引文件：普通 `.ht2` 索引直接运行；检测到 `.ht2l` 大型索引时自动向 HISAT2 包装器传入 `--large-index`，从而调用 `hisat2-align-l`。
- `annotation_gtf` 必须填写提前转换并检查好的 GTF。
- 同一个 `reference_name` 必须使用相同的索引前缀和 GTF。
- 不同参考基因组的表达矩阵分别生成，不会直接合并。
- `sample` 和 `reference_name` 只能包含字母、数字、点、下划线和连字符。

## 3. 运行方式

先生成作业并检查，不提交：

```bash
bash main.sh sample_info.txt results all --no-submit
```

确认各样本目录中的 `.slurm` 文件后正式提交：

```bash
bash main.sh sample_info.txt results all
```

每个样本会提交为一个独立的 20 核作业。全部样本成功后，主脚本通过 Slurm `afterok` 依赖自动提交汇总作业。脚本不会添加任何 Slurm 内存限制。

## 4. 选择性运行

```bash
### 只重新运行 fastp
bash main.sh sample_info.txt results 01

### 只重新运行 HISAT2，运行前检查步骤 01 输出
bash main.sh sample_info.txt results 02

### 从 HISAT2 开始，依次运行后续全部步骤
bash main.sh sample_info.txt results --from 02

### 只重新生成统计表和 MultiQC
bash main.sh sample_info.txt results 04

### 只重新生成 StringTie count 矩阵
bash main.sh sample_info.txt results 05

### 只重新生成 gene TPM/FPKM 矩阵
bash main.sh sample_info.txt results 06

### 只重新生成 transcript TPM/FPKM 矩阵
bash main.sh sample_info.txt results 07
```

步骤定义：

1. `01`：fastp，使用软件默认过滤参数，仅指定双端输入输出、20 线程及 JSON/HTML 报告。
2. `02`：从 GTF 提取已知剪接位点，HISAT2 使用 `--dta` 比对，并根据 `.ht2`/`.ht2l` 自动选择 small/large index 模式，随后排序并建立 BAM 索引。
3. `03`：StringTie 使用 `-G -e -B`，只定量 GTF 中已有转录本。
4. `04`：汇总 fastp 和 HISAT2 指标并运行 MultiQC。
5. `05`：按参考基因组使用 prepDE 公式生成 gene/transcript count 矩阵，不依赖外部 `prepDE.py3` 文件。
6. `06`：按参考基因组生成 gene TPM/FPKM 矩阵。
7. `07`：按参考基因组生成 transcript TPM/FPKM 矩阵。

## 5. 主要输出

```text
results/
├── sample_info_validation.log
├── samples/
│   └── sample_01/
│       ├── job/
│       ├── 01_fastp/
│       │   ├── sample_01.clean.R1.fastq.gz
│       │   ├── sample_01.clean.R2.fastq.gz
│       │   ├── sample_01.fastp.json
│       │   └── sample_01.fastp.html
│       ├── 02_hisat2/
│       │   ├── sample_01.sorted.bam
│       │   ├── sample_01.sorted.bam.bai
│       │   ├── sample_01.hisat2_summary.txt
│       │   └── sample_01.flagstat.txt
│       └── 03_stringtie/
│           ├── sample_01.stringtie.gtf
│           └── sample_01.gene_abund.tsv
└── summary/
    ├── 04_stats/
    │   ├── 04_fastp_stats.tsv
    │   ├── 04_hisat2_alignment_stats.tsv
    │   ├── 04_qc_alignment_summary.tsv
    │   └── multiqc/
    ├── 05_counts/
    │   └── reference_name/
    │       ├── reference_name.gene_count_matrix.tsv
    │       └── reference_name.transcript_count_matrix.tsv
    └── 06_abundance/
        └── reference_name/
            ├── reference_name.gene_tpm_matrix.tsv
            ├── reference_name.gene_fpkm_matrix.tsv
            ├── reference_name.transcript_tpm_matrix.tsv
            └── reference_name.transcript_fpkm_matrix.tsv
```

流程先临时生成未排序 BAM，再使用 20 线程排序；排序成功后自动删除临时 BAM。最终只保留排序 BAM、BAM 索引和 gzip 格式的 clean FASTQ。

## 6. count 矩阵说明

步骤 05 的 raw count 是按照 prepDE 公式，根据 StringTie coverage、转录本外显子总长度和平均 read 长度推算的 hypothetical read count。流程从每个参考基因组组内所有样本的 fastp 清洗后数据计算加权平均 read 长度，并将其记录在：

```text
results/summary/05_counts/05_prepde_metadata.tsv
```

该矩阵可用于 DESeq2/edgeR，但它不是 featureCounts 对 BAM 直接计数的结果。正式差异分析前应确认同一比较中的样本属于同一参考基因组，并检查 count 分布和样本设计。
