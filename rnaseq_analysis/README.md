# StringTie表达定量和全基因差异表达

2026-10-20发布，来源为已用于六材料正式分析的超算流程。代码分为两个可单独运行的目录：

| 流程 | 输入 | 输出 | 入口 |
|---|---|---|---|
| map_and_quantify | 双端非链特异RNA-seq FASTQ、预建HISAT2索引、固定参考GTF、六列样本表 | 每样本StringTie GTF与丰度表、每材料gene/transcript counts和TPM/FPKM | [定量说明](map_and_quantify/README_RNASEQ_PIPELINE.md) |
| DEG_stringtie | 同批StringTie全基因counts和TPM、样本分组、比较表及可核验的上游来源 | 每材料联合DESeq2模型、逐比较DEG、LRT、标准化counts与QC图 | [差异分析说明](DEG_stringtie/README.md)、[逐列配置指南](DEG_stringtie/CONFIG_GUIDE.md) |

## 方法

StringTie使用 `-G -e -B -A` 在固定注释下定量。counts由本流程实现prepDE公式：每转录本 `ceil(cov × exon_length / read_length)` 后按gene累加；read长度为同材料样本的fastp过滤后总碱基数/总reads数取整。TPM来自同批StringTie的gene abundance表。这是StringTie估计counts，不是featureCounts直接计数，也没有在此流程中调用官方prepDE脚本。

每材料所有条件和重复联合拟合 `~condition`；再按比较表提取Wald检验，另做条件效应LRT。默认至少2个样本counts≥10，DEG为 `padj<0.05` 且 `|log2FC|≥1`。BH校正在每个比较的全基因背景下完成；随后再提取HSF等关注基因。不能把TPM送入DESeq2，也不对counts人为加1。

同材料使用相同参考GTF；对照条件固定为 `CK`，每组至少2个独立生物学重复。配对、批次或交互设计需要调整模型，不能直接套用。当前定量入口面向双端非链特异文库，链特异数据必须先调整比对和StringTie参数。

## 从头运行

在Linux/Slurm上配置 `map_and_quantify/environment.yml` 和 `DEG_stringtie/environment.yml`。上游环境名为 `rnaseq`，下游可以使用单独的Python和R解释器，详见其runtime配置。

```bash
cd map_and_quantify
cp sample_info.example.txt sample_info.txt
### 将样本表中的示例路径换成真实绝对路径后，预览并提交
bash main.sh sample_info.txt results all --no-submit
bash main.sh sample_info.txt results all
```

上游每样本独立作业，全部成功后以afterok依赖运行矩阵汇总；不添加Slurm内存限制。当前原始模板默认样本20核、汇总2核、分区 `x86_64_4cpu`；若集群不同，在主控生成脚本的配置位置修改。

定量完成后，进入 `DEG_stringtie`，按指南填写 `config/materials/MaterialA/00_*.tsv` 和材料清单，准备Arial字体：

```bash
cd ../DEG_stringtie
conda activate rnaseq_de
### 如复用分开的环境，先export PIPE_PYTHON及PIPE_RSCRIPT为对应解释器路径
bash submit_stringtie.sh config/materials/sample_info.txt runs/check --no-submit
bash submit_stringtie.sh config/materials/sample_info.txt runs/formal
```

已完成StringTie定量的项目可以直接运行第二部分，无需重复比对。审核需要上游BAM、StringTie GTF头部命令、丰度表和fastp JSON，不能只有两张矩阵。自动生成配置要求上游样本表保存在 `map_and_quantify/sample_info*.txt`，结果布局为 `map_and_quantify/results/`；其他布局应手工配置并核对路径。

## 选择步骤与复用

上游支持 `all`、`01`—`07`、`--from NN` 和 `--no-submit`。下游支持 `all`、`01`—`04`、`--from NN`；每个材料一个独立Slurm作业。更换样本、参数或代码时使用新的批次目录。

```bash
cd DEG_stringtie
bash main.sh config/materials/MaterialA results/MaterialA all
bash main.sh config/materials/MaterialA results/MaterialA --from 03
```

完整定量和差异分析代码来自已运行版本；发布时仅调整下游默认解释器路径、加入可移动字体配置与示例、完善说明，并清理部分文件末尾的空白行。没有改变计数公式、ID处理、DESeq2模型或统计阈值。上游原有长脚本保留原实现，便于与实际运行版本核对；发布源文件指纹见 `SOURCE_MANIFEST.tsv`。

仓库保存代码、依赖清单和配置模板。实际FASTQ/BAM、表达矩阵、差异分析结果、用户服务器配置及Arial字体在运行机器上保存。

## 验证

下游提供505基因、10样本的模拟测试，检查已知上调/下调、零基线诱导、全零和低计数基因、输出图件，以及三种错误输入的拒绝行为：

```bash
cd DEG_stringtie
source config/runtime.sh
"$PIPE_PYTHON" tests/01_fixture.py --input . --outdir tests/work_release
sbatch tests/run_validation.sh tests/work_release tests/validation_release
```

测试前准备 `assets/Arial.ttf`，在有Slurm、fontconfig、flock及依赖的Linux环境运行。每次使用新的模拟和验证目录。发布检查记录见 [VALIDATION.md](VALIDATION.md)。
