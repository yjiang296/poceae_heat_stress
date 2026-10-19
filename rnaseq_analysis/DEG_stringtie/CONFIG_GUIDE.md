# StringTie DEG 流程配置填写指南

适用于本目录当前版本的 `00_build_config.py`、`01_audit.py`、`02_matrices.py` 和 `03_deseq2.R`。新项目可以从已有比较配置自动生成，也可以手动填表。当前模型为每材料独立的 `~condition`，对照组名称必须为 `CK`，每组至少两个独立生物学重复；配对、批次校正等设计需要修改模型，不能仅增加配置列。

## 1. 文件放在哪里

建议为新项目另建配置目录，避免覆盖本批六材料配置：

```text
config/new_materials/
├── sample_info.txt
└── MaterialA/
    ├── 00_samples.tsv
    ├── 00_contrasts.tsv
    ├── 00_sources.tsv
    └── 00_original_comparisons.tsv
```

四张 `00_*.tsv` 都必须保留下面给出的英文表头，使用真实制表符 Tab 分隔，不能用空格对齐代替。不要在表内加入说明行、合并单元格或空值占位符。路径填写超算上实际存在的绝对路径，不含空白，不填写本地 Mac 路径或未展开的 `$变量`、`~`。

以下示例假设新定量项目为 `/project/map_and_quantify`，新 DEG 流程为 `/project/DEG_stringtie`，材料名为 `MaterialA`。这些都是示意路径，必须换成实际值。仓库提供 `config/example/` 的四张表和 `config/sample_info.example.txt`，复制后填写自己的样本及绝对路径；实际材料配置不随代码上传。

## 2. 样本表：00_samples.tsv

一行对应一个独立样本。同一对照用于多个比较时，只写一次。

| 列名 | 填写内容 | 要求 |
| --- | --- | --- |
| `sample_id` | StringTie 定量使用的真实样本名 | 唯一；与 counts、TPM 的样本列名及下述文件名完全一致，区分大小写 |
| `condition` | 样本所属条件 | 对照写 `CK`；处理可写 `heat0_5h`、`heat1h`；同组重复使用相同条件名 |
| `bam` | 该样本参与 StringTie 定量的 BAM 绝对路径 | 文件存在且非空；路径必须与 StringTie GTF 头部命令中的 BAM 参数一致 |
| `annotation_gtf` | 定量时 `-G` 指定的参考注释 GTF 绝对路径 | 同材料所有行相同；必须与 StringTie GTF 头部 `-G` 参数的路径字符串一致 |

六样本、两处理条件的完整示例：

```tsv
sample_id	condition	bam	annotation_gtf
A_CK_1	CK	/project/map_and_quantify/results/samples/A_CK_1/02_hisat2/A_CK_1.sorted.bam	/project/reference/MaterialA.gtf
A_CK_2	CK	/project/map_and_quantify/results/samples/A_CK_2/02_hisat2/A_CK_2.sorted.bam	/project/reference/MaterialA.gtf
A_heat0_5h_1	heat0_5h	/project/map_and_quantify/results/samples/A_heat0_5h_1/02_hisat2/A_heat0_5h_1.sorted.bam	/project/reference/MaterialA.gtf
A_heat0_5h_2	heat0_5h	/project/map_and_quantify/results/samples/A_heat0_5h_2/02_hisat2/A_heat0_5h_2.sorted.bam	/project/reference/MaterialA.gtf
A_heat1h_1	heat1h	/project/map_and_quantify/results/samples/A_heat1h_1/02_hisat2/A_heat1h_1.sorted.bam	/project/reference/MaterialA.gtf
A_heat1h_2	heat1h	/project/map_and_quantify/results/samples/A_heat1h_2/02_hisat2/A_heat1h_2.sorted.bam	/project/reference/MaterialA.gtf
```

不能把旧比较表中的样本别名直接当作 `sample_id`。例如本项目部分旧名称与定量矩阵名称不同，自动生成器通过 BAM 路径找回真实样本名。

### 当前审核依赖的上游文件布局

审核程序从 `bam` 路径向上取两级作为样本目录，然后查找以下文件。例如 `A_CK_1`：

```text
/project/map_and_quantify/results/samples/A_CK_1/
├── 01_fastp/A_CK_1.fastp.json
├── 02_hisat2/A_CK_1.sorted.bam
└── 03_stringtie/
    ├── A_CK_1.stringtie.gtf
    └── A_CK_1.gene_abund.tsv
```

GTF 第一行须保留真实 StringTie 命令，含 `-G`、`-e` 和相符的 BAM 路径；第二行用于记录版本。fastp JSON 必须包含 `summary.after_filtering.total_bases` 和 `total_reads`，用于记录读长。`annotation_gtf` 指参考注释，不是这里每个样本的 `.stringtie.gtf`。

**只有矩阵、缺少这些文件，或目录组织不同的新项目，不能仅靠填四张表直接运行当前版本。** 应先适配审核程序的路径读取方式；不要为通过审核伪造 GTF 命令或改写来源记录。迁移原始定量文件后若绝对路径改变，也需要明确处理来源路径映射。

## 3. 比较表：00_contrasts.tsv

一行对应一个要提取的差异表达比较。

| 列名 | 填写内容 | 要求 |
| --- | --- | --- |
| `contrast_id` | 比较名称，用于结果文件命名 | 表内唯一；建议仅用英文字母、数字、下划线、点和短横线 |
| `numerator` | 分子条件，通常为处理组 | 必须与样本表的某个 `condition` 完全一致 |
| `denominator` | 分母条件，通常为对照组 | 必须与样本表的某个 `condition` 完全一致；共同 CK 设计写 `CK` |

```tsv
contrast_id	numerator	denominator
MaterialA_heat0_5h_vs_CK	heat0_5h	CK
MaterialA_heat1h_vs_CK	heat1h	CK
```

`log2FoldChange` 的方向为 **numerator / denominator**；上例正值表示热处理相对 CK 上调。条件名建议以字母开头、仅使用字母数字和下划线，避免空格及特殊字符。

比较表出现的全部条件的集合，必须等于样本表中条件的集合。不要保留没有参与任何比较的额外条件。模型使用样本表中的所有样本联合拟合；删去某条件的样本后重跑会改变模型，不能视为仅少输出一张比较表。

## 4. 来源表：00_sources.tsv

每材料仅一行，指定全基因矩阵及其来源。

| 列名 | 填写内容 | 要求 |
| --- | --- | --- |
| `material` | 材料名称，例如 `MaterialA` | 与材料目录名、外层材料清单一致；用于输出文件前缀 |
| `counts` | 全基因 StringTie 估计 counts 矩阵的绝对路径 | 未做 DESeq2 标准化，非 TPM，非 featureCounts，非仅 HSF 子集 |
| `tpm` | 同批 StringTie 定量的全基因 TPM 矩阵绝对路径 | 与 counts 使用相同样本来源和基因集合 |
| `quant_root` | 上游定量流程根目录 | 例如 `/project/map_and_quantify`；不是 `results/summary`，也不是新 DEG 目录 |
| `method` | 当前实现要求的来源标识 | 固定为 `STRINGTIE_PREPDE_FORMULA`，前提是实际来源确实符合下面的方法说明 |

```tsv
material	counts	tpm	quant_root	method
MaterialA	/project/map_and_quantify/results/summary/05_counts/MaterialA/MaterialA.gene_count_matrix.tsv	/project/map_and_quantify/results/summary/06_abundance/MaterialA/MaterialA.gene_tpm_matrix.tsv	/project/map_and_quantify	STRINGTIE_PREPDE_FORMULA
```

本项目 counts 由上游脚本按每转录本 `ceil(cov × exon_length / read_length)` 再按 gene 累加得到，TPM 来自同批 StringTie gene abundance。`method` 是来源声明，不能把其他工具的矩阵填入后仅修改这个字符串。当前审核记录来源与一致性，但不会从 GTF 重新计算矩阵，不能代替对实际生成方法的确认。

审核会记录 `quant_root/scripts/` 中找到的 `05_*counts.py`、`05_run_prepde.py` 的 SHA256；这些脚本未找到时，当前代码不会因此单独报错，因此新项目应保留实际生成矩阵的实现及参数以便追溯。

两个矩阵须满足：

- 第一列名为 `gene_id`，其余包含样本表全部 `sample_id` 对应的列；额外样本列不会纳入该材料模型。TPM 中额外的 `gene_name` 列可以保留。
- 基因 ID 唯一；只允许程序去除已知 `gene:` 前缀，去除后不能产生重复。两矩阵的基因集合必须一致，不会为缺失基因补零。
- 行列顺序可以不同，程序按 ID 匹配。不要用调换列标题的方法强行对应样本。
- counts 为非负有限整数，最大值不超过 2147483647，每个样本总 counts 大于零；TPM 为非负有限数，不使用 `NA`、空白等替代缺失值。
- 任一样本出现同一基因 `count=0` 且 `TPM≥1` 时，步骤 02 会中止并输出 QC。不要修改数值规避检查，应核对是否来自同批定量。

## 5. 原比较校验和表：00_original_comparisons.tsv

这张表用于追溯旧比较文件，不是新模型的分组依据。分组与比较分别由 `00_samples.tsv` 和 `00_contrasts.tsv` 决定。

| 列名 | 填写内容 |
| --- | --- |
| `path` | 旧 `comparison.*.tsv` 文件的绝对路径 |
| `sha256` | 该旧文件内容对应的 SHA256，64 位十六进制字符串 |

**全新项目没有旧比较文件时，保存只有下面一行表头的文件即可，不要虚构路径或校验和：**

```tsv
path	sha256
```

有旧比较文件时优先使用配置生成器自动记录。手动填写可在超算用 `sha256sum /实际路径/comparison.xxx.tsv` 获取校验和；将输出中的路径与校验和值分别填入 `path`、`sha256` 列，注意命令输出的顺序与表格列顺序不同。Mac 可用 `shasum -a 256`。

当前步骤 01 对这张表自身计算校验和，并未逐行重新打开旧文件验证其 SHA256。表内值是生成配置时的历史记录，不能将其描述为运行时完成了所有旧文件的完整性复核。

## 6. 外层材料清单：sample_info.txt

这里一行是一个材料，不是一个 RNA-seq 样本。提交器按材料创建独立作业。

```tsv
material	config_dir
MaterialA	/project/DEG_stringtie/config/new_materials/MaterialA
```

`material` 必须唯一；`config_dir` 指向包含上述四张表的目录。增加材料时增加一行，并为它单独准备四张表。材料之间不会合并拟合。

## 7. 有旧比较配置时自动生成

在新 DEG 流程根目录执行，替换成真实路径；使用尚未存在材料子目录的新输出位置：

```bash
source config/runtime.sh
"$PIPE_PYTHON" scripts/00_build_config.py \
  --input /project/DEG_identify \
  --outdir config/new_materials \
  --quant-root /project/map_and_quantify
```

生成器读取旧目录的 `comparison.*.tsv`，排除 `comparison.example.tsv`。旧表要求恰好包含五列：`sample_id`、`condition`、`bam`、`annotation_gtf`、`reference_name`，其中 `condition` 取值为 `control` 或 `treatment`，每组至少两样本。

同时读取定量根目录的 `sample_info*.txt`。当前上游样本表列为 `sample`、`raw_fq_gz_R1`、`raw_fq_gz_R2`、`hisat2_index_prefix`、`annotation_gtf`、`reference_name`，用空白或 Tab 分隔；生成器使用其中真实样本名、注释和材料名，并按既定 `results/samples/<sample>/02_hisat2/<sample>.sorted.bam` 路径匹配旧表。

同一材料各比较的 CK 必须是同一组 BAM。生成后，对照条件改为 `CK`，处理条件命名为旧比较文件去掉 `comparison.` 和 `.tsv` 后的名称。材料名采用上游定量样本表的 `reference_name`，而非旧表别名。全基因矩阵路径按第 4 节示例的目录约定生成，仍需核对实际文件存在。

生成器拒绝覆盖已存在的材料配置目录。没有旧比较表时，直接按第 2—6 节手动填写，无需为调用生成器而虚构旧配置。

## 8. 填完后如何检查和提交

先核对 `config/runtime.sh` 的软件路径、分区、统计阈值，以及 Arial 字体配置。下面命令在流程根目录运行；按集群要求在获分配的计算资源内执行读取矩阵等实际检查。

```bash
bash main.sh config/new_materials/MaterialA results/config_check_MaterialA 01
bash main.sh config/new_materials/MaterialA results/config_check_MaterialA 02
```

步骤 01 核验文件和来源，步骤 02 核验矩阵与 counts/TPM 一致性。这两个步骤均通过后，可生成提交脚本检查，再提交正式批次：

```bash
bash submit_stringtie.sh config/new_materials/sample_info.txt runs/new_project_check --no-submit
bash submit_stringtie.sh config/new_materials/sample_info.txt runs/new_project_formal
```

`--no-submit` 只生成冻结配置、代码和 Slurm 脚本，不运行完整的步骤 01—02 输入审核。每次提交或 dry-run 使用新的批次目录；不要将已有正式结果目录作为新项目输出。实际正式作业仍会再次从步骤 01 开始核验。

常见报错对应关系：

| 报错或现象 | 优先检查 |
| --- | --- |
| 缺少样本列 | `00_samples.tsv` 的 `sample_id` 是否与两矩阵列名相同，是否误用了旧别名 |
| 比较与条件不匹配 | `numerator`、`denominator` 是否拼写一致，是否存在未参与比较的样本条件 |
| 注释路径不符 / BAM 来源不符 | 配置路径与 GTF 头部真实命令是否一致，是否迁移过目录 |
| 输入缺失 | BAM 上两级推导的样本目录及 fastp、StringTie 文件名是否符合第 2 节 |
| counts 与 TPM 基因集合不同 | 是否使用相同注释与同一批全基因定量，是否提前筛选过其中一个矩阵 |
| 找不到参考水平 CK | 样本表是否将对照明确命名为 `CK`，而不是 `control` 或其他名称 |
