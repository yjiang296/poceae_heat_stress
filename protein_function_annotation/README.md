# 全蛋白集功能注释：InterProScan 5 + eggNOG-mapper

为已有基因模型对应的蛋白 FASTA 添加功能注释。输入是蛋白，不是基因组 DNA 或 CDS；本流程不预测基因、不生成新版 GFF3，不重新选择最长转录本。相同序列但不同 ID 的蛋白仍分别保留。

## 1. 配置与软件

已验证组合：InterProScan **5.78-109.0**，Java **11**，Python **3.11**；eggNOG-mapper **2.1.13**，数据库 **5.0.2**，DIAMOND **2.0.15**。工作流 Python 仅依赖标准库，计算节点需要 Bash 和 `flock`，批量调度需要 Slurm。版本更换后应重新验证。

```bash
cd protein_function_annotation
cp config.example.json config.json
cp sample_info.example.txt sample_info.txt
```

编辑 `config.json` 的软件、Python 和数据库绝对路径，以及分区和线程数。模板中的 `/absolute/path/to/` 必须替换为真实路径，不自动安装软件或下载数据库。安装方式及环境原理见 [INSTALLATION.md](INSTALLATION.md)。

脚本为工具子进程设置独立 PATH；InterProScan 还设置 JAVA_HOME，eggNOG 设置 EGGNOG_DATA_DIR。因此无需手动 `conda activate`，不会更改当前终端环境。主控环境名可能显示 base，应以工具日志和实际可执行路径为准。

## 2. 样本表

两列，以空白或 Tab 分隔，路径必须绝对且不能包含空白：

```text
### 材料名 蛋白FASTA
B73 /absolute/path/to/B73.protein.fa
Bowman /absolute/path/to/Bowman.protein.fa
```

支持普通 FASTA 和 gzip 文件；蛋白 ID 取标题第一个空白前字段，同材料内必须唯一。材料名仅含英文字母、数字、点、短横线或下划线，并以字母/数字开头。跳过空行及 `#` 开头行。

六材料模板只是示例，不内含真实序列。可添加其他材料；每材料独立提交。用于基因层面比较时，应提供已选定代表转录本的蛋白集，本流程不会擅自去重。TEsorter 过滤前后两套输入需使用不同输出目录。

## 3. 批量运行

在服务器登录节点执行外层提交器，不要将 `submit_annotation.sh` 本身再交给 sbatch。

```bash
# 生成作业，检查输入表与 Slurm 脚本；不会提交
bash submit_annotation.sh sample_info.txt results/run_v1 --no-submit

# 正式提交表中全部材料
bash submit_annotation.sh sample_info.txt results/run_v1 --submit

# 也可以仅选择部分材料
bash submit_annotation.sh sample_info.txt results/subset_v1 \
  --samples Bowman,B73 --threads 8 --partition x86_64_4cpu --submit
```

默认每材料 8 线程，每个材料一个独立作业，作业内按 01 → 02 → 03 → 04 执行。Slurm 头部为单节点 `-N 1 -n THREADS`，不设置 `--mem`。默认分区名来自原验证集群，在其他集群必须按实际配置修改。

每次生成独立作业脚本和配置快照，提交记录在 `results/run_v1/submissions.tsv`。输入、代码、线程数、配置改变后应使用新的输出目录。仅出现“仅生成，未提交”不是已经提交。

## 4. 四个步骤

| 步骤 | 工具/工作 | 输出 |
|---|---|---|
| 01 | 输入检查、ID 对照、分块 | 原 ID 映射、规范蛋白 FASTA、chunks |
| **02** | **InterProScan**：扫描家族、结构域和相关 GO | 分块原始 TSV/JSON |
| **03** | **eggNOG-mapper**：同源搜索和功能转移 | annotations、seed_orthologs、hits |
| 04 | 恢复原 ID，合并两套注释 | members.tsv、domains.tsv、summary.json |

InterProScan 默认扫描 Pfam、ProSiteProfiles、ProSitePatterns、Gene3D、SUPERFAMILY、CDD、PANTHER；每块 1,000 条蛋白，块间串行、块内使用配置线程数。参数 `-f tsv,json -goterms -dp`；`-dp` 关闭远程预计算查询。JSON 用于核验包括无命中在内的所有 query 均被处理。未默认启用所有可选库、SignalP、TMHMM 或通路映射。

eggNOG 对每材料全部可用蛋白一次搜索：`--tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --sensmode very-sensitive --dmnd_iterate yes --block_size 1`。其余阈值沿用固定版本默认值，实际命令保存在日志和工具输出头部。eggNOG 的 PFAMs 是转移注释，与 InterProScan 的直接 Pfam 扫描分开保存。

01 统一大小写、去除一个末端 `*` 并记录；允许标准 20 种氨基酸和 X。默认 `invalid_policy: fail`：发现空序列、内部终止符、U/O/B/Z/J 等非法字符时输出问题记录并停止。检查后可在首次运行前改为 `quarantine`，隔离成员仍保留在最终表。重复 ID 总是失败；输入全为核酸字母时拒绝运行。

内部 ID 为 `<sample>__p000000001`；最终结果提供 `original_id` 和原始标题。回接 GFF3 或其他分析使用实际 ID 映射，尤其注意 OrthoFinder 可能替换冒号等字符。

## 5. 结果与解释

每材料结果位于 `<输出根目录>/<材料>/04/`：

- `04_<sample>_members.tsv`：每条输入蛋白一行，包含 original_id、输入状态、InterPro/Pfam/GO、eggNOG 描述和 GO/KEGG 等原始字段。
- `04_<sample>_domains.tsv`：直接扫描的命中坐标、数据库、分值；重叠签名不能简单视为不同结构域。
- `04_<sample>_summary.json`：处理数量和命中计数。
- `04_complete.json`：成功完成并通过校验的记录。

`ipr_status` 区分 hit、no_match、not_run_quarantined；`eggnog_status` 区分 annotation_row、seed_without_annotation、no_reported_seed 和未运行。无命中不等于新基因或没有功能；注释名称也不等于实验验证。`seed_both_coverage_ge80` 是双侧覆盖≥80%的描述性标记，不是功能可靠性的自动判定。本流程不做富集检验，不根据注释自动判定 TE、谱系限制或适应性创新。

## 6. 续跑与单步执行

```bash
# 修复错误后续跑；已完成且校验一致的步骤自动跳过
bash submit_annotation.sh sample_info.txt results/run_v1 --samples Bowman --resume --submit
# 从 eggNOG 开始，要求已有上游结果
bash submit_annotation.sh sample_info.txt results/run_v1 --samples Bowman --from 03 --resume --submit
# 只合并，要求 01、02、03 完整
bash submit_annotation.sh sample_info.txt results/run_v1 --samples Bowman --steps 04 --resume --submit
```

提交器检查旧作业是否仍在队列，材料目录还使用运行锁。InterProScan 可复用已完成且校验一致的块；失败尝试保留。eggNOG 失败后整步在新 attempt 目录重跑，不支持 DIAMOND 内部进度续跑。结果校验失败时停止，不自动覆盖。代码、配置或位置改变时不得强改完成签名。

交互计算节点也可运行：

```bash
bash main.sh /absolute/path/to/A.protein.fa /absolute/path/to/results/A all \
  --sample A --config config.json --threads 8
```

不要在登录节点运行全量 `main.sh`。每材料的 jobs 保存配置和 Slurm 脚本，logs 保存调度输出，每步骤有 run.log、version.log、requirements.txt 和完成标记。完成标记包含绝对路径，下载结果仅供查看，不直接在新位置续跑。

## 7. 测试与状态

```bash
bash -n main.sh
bash -n submit_annotation.sh
python3 -m unittest discover -s tests -v
squeue -u "$USER"
sacct -j 作业号 --format=JobID,State,ExitCode,Elapsed -X
```

测试使用合成的小输入，不需要安装 InterProScan 或 eggNOG。真实工具端到端验证记录见 [VALIDATION.md](VALIDATION.md)。完整蛋白集的耗时和内存需求需要实测，不能按两条测试蛋白线性外推。
