# StringTie统一定量来源的全基因差异表达流程

本流程部署于 `rnaseq_analysis/DEG_stringtie`，与旧 `DEG_identify` 分开。复用已经完成的 StringTie `-G -e` 定量结果，不调用 featureCounts，不重新比对，也不覆盖旧结果。基因 counts 和 TPM 必须来自同一批 StringTie 定量。

## 本批次与统计定义

- 原24个比较配置对应6材料、60个独立样本；共享CK只保留一次。
- 基于原BAM路径和定量样本表映射样本别名，禁止凭列顺序匹配。
- 每材料的10个样本联合拟合 DESeq2 `~condition`，再提取4个heat/CK Wald比较，并做全部条件的LRT；不按重复编号配对，不跨材料合并建模。
- 与旧流程每次仅取4个样本拟合不同，联合模型共享离散度估计；因此新旧结果的差异不能全部归因于计数工具。
- 全基因预过滤：至少2个样本count≥10；BH校正在全基因检验背景中进行，每个比较单独校正；`padj<0.05`且`|log2FC|≥1`为后筛选，四个比较取并集不具有自动的跨时间总体FDR保证。
- 不对counts人为加1。输出未收缩log2FC与SE，零基线导致的大倍数变化需结合表达量和SE解释。
- `LOW_COUNT_FILTERED`、`NO_VALID_PVALUE`、`NO_VALID_PADJ` 与 `NO_DETECTED_RESPONSE` 分开；所有输入基因保留在全结果表。
- LRT只检验单材料条件效应；当前0h CK设计不能完全排除自然时间变化，不是材料×时间交互检验。

## 当前StringTie counts的准确来源

上游 `map_and_quantify/scripts/05_run_prepde.py` 和 `05_extract_sample_counts.py` 实现prepDE公式：每转录本 `ceil(cov × exon_length / read_length)` 后按gene累加。当前流程复用这些已完成矩阵；不把它们描述为重新运行官方prepDE获得的原始片段计数。read长度来自fastp过滤后碱基数/读数，按材料汇总取整。审核记录这项参数、上游实现SHA256、矩阵SHA256和StringTie命令。

TPM来自同批 `StringTie -A` 的gene abundance汇总。DESeq2标准化counts和TPM仍是不同单位；统一来源不要求其倍数变化完全一致。仅去除已知的`gene:`前缀，并保存双向ID对应表；任何碰撞或未匹配均报错。

## 准备配置

软件路径、分区和统计阈值集中在 `config/runtime.sh`。发布版默认从PATH查找 `python3`、`Rscript`，也可通过 `PIPE_PYTHON`、`PIPE_RSCRIPT` 指定现有环境的绝对路径。已有兼容环境可直接复用；新机器可用 `conda env create -f environment.yml` 创建环境后 `conda activate rnaseq_de`。环境文件是依赖清单，不是版本锁定文件，实际版本由运行结果记录。

将自己有权使用的Arial字体放在 `assets/Arial.ttf`，详见 [字体配置](assets/README.md)。字体文件不随仓库分发；提交器会为冻结快照生成独立字体配置。Shell主控需要Linux的 `flock`，批量提交需要Slurm。

没有旧比较表时，复制 `config/example/` 到 `config/materials/MaterialA/`，修改四张表，并将 `config/sample_info.example.txt` 复制为 `config/materials/sample_info.txt` 后填写配置绝对路径。下述自动生成命令仅适用于已有兼容旧比较文件的项目。

```bash
source config/runtime.sh
"$PIPE_PYTHON" scripts/00_build_config.py \
  --input ../DEG_identify --outdir config/materials \
  --quant-root ../map_and_quantify
```

配置生成器只读取`comparison.*.tsv`正式比较，排除example。它要求同材料所有比较使用相同CK样本集合，适用于独立重复的共同CK多处理设计；其他实验设计需修改模型和配置，不能直接套用。

新项目复用时，可重建配置或填写每材料的四个`00_*.tsv`：样本表、比较表、来源表、原比较校验和表。**逐列解释、完整填表示例、自动生成方法及提交前检查，见 [CONFIG_GUIDE.md：配置填写指南](CONFIG_GUIDE.md)。** 没有旧比较文件时，原比较校验和表只保留表头即可。当前审核依赖上游的目录布局及可核验的StringTie GTF、BAM和fastp JSON；仅提供counts、TPM矩阵不足以直接运行。具体限制见指南。

## 提交与选择步骤

```bash
bash submit_stringtie.sh config/materials/sample_info.txt runs/check_20261014 --no-submit
bash submit_stringtie.sh config/materials/sample_info.txt runs/formal_20261014
```

每个材料一个作业；已完成的样本定量不再提交重复作业。使用`-J/-p/-N/-n/-o/-e`，不设置内存限制。每次使用新的批次目录，冻结代码、配置及字体，作业ID即时写入`00_submission.tsv`。失败时不要重复提交整个批次。

```bash
bash main.sh config/materials/B73 results/B73 all
bash main.sh config/materials/B73 results/B73 03
bash main.sh config/materials/B73 results/B73 --from 03
```

正式批次重跑应使用其`code/main.sh`与`configs/B73`，保持与已冻结版本一致；可将上述命令放入对应材料的Slurm脚本再sbatch。每步检查上游SUCCESS；重跑上游会使旧下游SUCCESS失效。日志含时间和作业ID，不覆盖历史；同一输出目录用flock防止并发写入。

## 输出

每材料输出`results/01_audit`、`02_matrices`、`03_deseq2`、`04_plots`。

- 01：样本与比较、输入SHA256、StringTie命令、版本与read长度；GTF记录大小和修改时间，未对全部大GTF逐个做完整内容哈希。
- 02：全基因counts/TPM、ID映射、逐样本定量QC；count=0且TPM≥1触发失败，不把源数据混用问题带入DESeq2。
- 03：每比较all_genes/deg/up/down、条件LRT、所有基因标准化counts、通过预过滤基因VST、模型RDS、汇总和统计参数。
- 04：PCA、样本相关图、四张火山图；每图有7×5in可编辑PDF、300dpi PNG、同名TSV和ggplot RDS；Arial与Wong色板。火山图只在显示时将-log10(padj)上限设为300。
- 各步骤保存日志和软件版本。03/04使用固定随机种子42。

全基因差异结果可供HSF或其他基因集提取。HSF特定分组、成员名单和表达图谱属于下游专题分析，不是此通用流程的依赖。

## 验证

`tests/01_fixture.py`生成可运行的500余基因模拟输入，包含已知上调、下调、零基线诱导、全零基因。`tests/02_check.py`验证方向、不可评价状态、图件及每图源数据。另以人为制造的counts/TPM不一致检查拒绝行为。正式提交前完成语法、dry-run、Slurm脚本与模拟全流程检查。

方法依据：[StringTie官方手册](https://www.ccb.jhu.edu/software/stringtie/index.shtml?t=manual)、[DESeq2官方教程](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)。
