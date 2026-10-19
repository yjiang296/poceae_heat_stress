# LAI 批处理 t6：诊断、使用与存储清理

核查日期：2026-10-20。服务器目录：
`/path/to/poceae_heat_stress/LAI`。

本次仓库版沿用当前本地文件名 `submit_ltr_pipeline.sh` 与 `submit_ltr_pipeline.md`。
服务器此前部署的文件名带 `_t6` 后缀；执行下面的无后缀命令前，应先将本目录文件下载到服务器目标目录。

### 参考样本表

`sample_info.txt` 包含 B73、05W002、0L503002 三个参考样本名，FASTA 路径使用占位符。以 Tab 分隔两列，注释行以 `#` 开头，不设置普通表头。
运行前必须将占位路径替换为实际 FASTA 绝对路径。
该文件只是三样本示例，不是完整 210 样本清单：

```bash
bash submit_ltr_pipeline.sh sample_info.txt --no-submit
```

批量分析仍使用服务器现有的 `sample_info_filter_organelle_contig.txt`。

## 1. B73 的失败原因已确认

旧 `results_LTR/B73/step3.err` 明确记录：

```text
ERROR: Fail to convert seq IDs to ≤ 13 characters! Please provide a genome with shorter seq IDs.
```

旧 `step3.out` 显示 LTR_retriever v3.0.4，历史运行报告最长 ID 为 32 字符；本次重新读取当前输入得到 350 条序列、162 条 ID 超过 13 字符、最长 30 字符。两者是不同时间的观测值，不强行视为一致；当前输入确实仍然超限。例子：`ptg000009l_162837965_162931838`。

解决方法是在 LTR_FINDER、LTRharvest、LTR_retriever **之前统一重命名**。新版生成 `s000000001` 形式的 ID，记录 `contig_id_map.tsv`，保留原 ID 与完整 header。序列内容及顺序不变，原始输入文件不修改。不能直接截断名称，否则可能重名；不能只在第三步改 FASTA，而继续沿用旧名称的候选文件。

旧脚本还有一个独立问题：缺少遇错退出机制，程序失败后仍执行时间统计，作业可能最终以 0 退出，导致 `afterok` 无法可靠识别失败。新版对子命令检查返回码，并验证最终文件和完成日志。

## 2. 已交付的文件与默认行为

- `submit_ltr_pipeline.sh`：独立新版，保留原 t5/t20。
- `submit_ltr_pipeline.md`：本说明。
- `test_lai_t6.py`：可重复运行的小型功能测试，不调用真实基因组分析。
- 结果默认写入当前目录下的 `results_LTR_t6`，不混入旧结果。
- 默认只生成作业，必须明确添加 `--submit` 才会提交。
- `--max-running` 默认 50，允许 1–50；每个数组元素处理一个样本，三个步骤依次运行。
- 同一输出目录有活动数组时，脚本拒绝再次提交，避免两次各 50 个叠加。该限制不覆盖其他输出目录、旧 t5 作业或人工直接 `sbatch`。
- 默认保留原分析参数：FINDER 30 线程，retriever 5 线程，分区 `x86_64_2cpu`。每个作业申请 1 个 task、30 CPUs；第三步实际使用 5 线程，因此会有 CPU 预留浪费。可调整 `--finder-threads` / `--retriever-threads`。
- 默认不指定 Slurm `--mem`，使用集群分区默认值；可通过 `--mem 64G` 等显式申请。
- 分析环境默认 `LTR_retriever`，默认使用 `~/miniconda3/etc/profile.d/conda.sh` 激活，可通过 `--conda-sh` 指定；兼容环境中的 Python 3.6。

并发上限控制总体同时运行数量，不能保证单个作业不会 OOM。若出现 OOM，应根据 `sacct` 的内存峰值调整 `--mem` 或选择合适分区。

## 3. 建议执行顺序

```bash
cd /path/to/poceae_heat_stress/LAI

# 先用 B73 验证真实软件完整执行；该命令会提交一个样本
bash submit_ltr_pipeline.sh sample_info_filter_organelle_contig.txt \
  --only B73 --max-running 1 --submit

# B73 完成后提交全表，成功且输入未变的 B73 会自动跳过
# 按要求最多同时 50 个样本
bash submit_ltr_pipeline.sh sample_info_filter_organelle_contig.txt \
  --max-running 50 --submit
```

如果存储空间有限，首轮可把第二条命令的 `50` 改为 `10` 或 `20`，观察真实峰值再调整后续批次。已运行数组也可通过 Slurm 的 `scontrol update JobId=数组ID ArrayTaskThrottle=20` 调低上限。

预检与本地功能测试：

```bash
bash submit_ltr_pipeline.sh sample_info_filter_organelle_contig.txt --no-submit
python3 test_lai_t6.py
```

210 个样本的文件路径与样本名重复检查已通过。样本表为两列：`sample_id FASTA绝对路径`，支持空行、注释行，路径中不能包含空白；相对 FASTA 路径按样本表目录解释。`.gz` 输入可直接解压流式重命名。

批次下保存脚本副本、样本表副本和 manifest，确保作业使用提交时的参数和代码。不要手工提交旧预览批次的 `run.slurm`；用主脚本 `--submit` 生成并登记当前版本。

## 4. 哪些文件清除，哪些保留

| 时点 | 清理文件 | 理由 |
|---|---|---|
| LTRharvest 返回成功且两份候选文件存在后 | `.suf .lcp .llv .esq .des .sds .ssp .prj .md5` 等已知 GenomeTools 索引，以及 `.mod.ltrTE.fa` | 后续 retriever 使用基因组与 SCN 候选，不使用这些索引或候选序列副本 |
| 最终结果验证通过后 | `.fa.masked`、`.fa.defalse`、`.fa.out.<数字>.LAI.LTR.ava.out`、已知 BED 交集临时表 | 当前 LAI 已完成，这些工作文件可重建 |
| 新版最终成功后 | 本次生成的 `genome.fa` 短 ID 副本 | 原始输入、映射及输入信息保留，可重新生成；旧流程原 FASTA 软链接不删除 |

**保留**：`.fa.out.LAI`、`.fa.pass.list`、相关 GFF3、`.fa.LTRlib.fa`、冗余 LTR 库、`.fa.out`、LTR 注释及汇总、所有 SCN 候选、日志、`contig_id_map.tsv`、输入来源信息、成功/失败标记、清理记录。没有列在清理清单里的文件全部保留。不会递归删除目录或跟随软链接。

成功判定要求：LTR_retriever 正常返回；日志包含 `All analyses were finished!`；LAI、pass.list、非冗余库和 RepeatMasker `.out` 都存在且非空；LAI 列结构正确且 `whole_genome` 行数值有限。只看到 `Runtime` 或一个存在但空的 LAI 文件不算成功。

每次清理前写清单，完成后更新状态。任何失败不会触发最终清理；但若 LTRharvest 已完成，其索引可能已按阶段清除。`--keep-intermediates` 可关闭新版两阶段清理，适合调试。

### muzhiyumi 实测

大小均按十进制 GB（10^9 bytes）表示，和 `du -h` 的 GiB 显示略有区别。

| 内容 | 大小 |
|---|---:|
| 清理前目录普通文件总量 | 29.881 GB |
| GenomeTools 索引 | 23.396 GB |
| masked FASTA | 2.256 GB |
| `.mod.ltrTE.fa` | 0.579 GB |
| 三个 LAI 比对临时文件合计 | 0.725 GB |
| 其余列入清理的临时表及 `.defalse` | 0.233 GB |
| **可释放合计** | **27.189 GB** |
| **清理后保留** | **2.692 GB** |

该样本 whole-genome raw_LAI=31.36、LAI=26.75；清理不会更改这些结果。

以该样本粗略外推，210 个样本最终保留约 565 GB；这不是峰值空间保证。分析期间还有基因组副本、索引和工具临时文件，而且失败样本不会自动清理。峰值需要用本批真实样本测量。

运行前请检查实际账号配额。`df -h /share` 显示共享文件系统总体空间，不能替代账号配额：

```bash
lfs quota -u "$USER" /share
```

本次没有为整个旧结果树得出完整占用统计；不要把 muzhiyumi 的大小当作全部样本的实测值。

## 5. 对旧结果做清理

新版只自动清理自己成功产生的结果。已有旧目录可先预览：

```bash
bash submit_ltr_pipeline.sh --cleanup-existing results_LTR \
  --only muzhiyumi B73 > cleanup_preview.json
```

实测预览中 muzhiyumi 可以释放 27,188,908,441 bytes；B73 因没有有效 LAI 自动跳过。

在没有排队/运行作业时应用清理：

```bash
bash submit_ltr_pipeline.sh --cleanup-existing results_LTR \
  --only muzhiyumi --apply
```

不指定 `--only` 会检查旧结果目录下全部样本，仅清理符合上述严格成功条件者。旧流程缺少文件锁，因此 `--apply` 会保守拒绝账号存在任何活动作业的情况；还应避免同时从其他终端启动分析。

本次只执行了预览，**未删除任何旧结果**。B73 等旧失败目录继续保留；即使新版成功，也不会自动删除旧目录，避免输入版本或分析结果混淆。

## 6. 输出、重跑与追溯

```text
results_LTR_t6/
  _batches/<批次ID>/
    manifest.json
    samples.original.txt
    submit_ltr_pipeline.sh
    run.slurm
    <数组ID>_<任务ID>.out/.err
  B73/
    complete.json
    run_<批次ID>_<随机后缀>/
      genome.fa.out.LAI
      genome.fa.pass.list
      genome.fa.LTRlib.fa
      genome.fa.out
      contig_id_map.tsv
      input_provenance.json
      step*.out/.err
      SUCCESS.json 或 FAILED.json
      cleanup.*.json
```

`complete.json` 指向最终成功目录并记录来源文件、大小、mtime 和全基因组 LAI。再次运行相同命令会跳过输入路径/大小/mtime 相同且成功结果仍有效的样本。这里不是输入文件完整 checksum 比对；若人工替换输入同时保留大小与 mtime，应另设输出目录强制重算。

失败样本重跑会创建新目录，不混用旧中间文件；这避免旧 ID、旧输入或部分输出污染结果，但不会从失败步骤断点续跑，旧失败目录可能仍占空间。需在完成诊断后单独处置。输入来源和序列 hash、参数、名称映射均保留；若原始 FASTA 将来被移动/删除，则不能仅凭映射恢复序列。

结果 GFF/LAI 中 contig 名称为短 ID。需要和原注释对接时，使用 `contig_id_map.tsv` 做明确字段映射；坐标未改变，不能对整个文件执行不加边界的文本替换。

本次验证覆盖重命名、重复/空序列拒绝、严格成功判定、定向清理及软链接保护、失败保留、成功跳过、真实子进程日志、重复提交拦截、数组上限与 Bash 语法；尚未重新运行 B73 的完整生物信息学计算。

## 7. 依据

- 服务器 B73 的 `step3.err` / `step3.out` 与当前输入 FASTA。
- muzhiyumi 完成日志、LAI 表和文件大小逐项核查。
- [LTR_retriever 官方说明](https://github.com/oushujun/LTR_retriever)：输入命名建议及最终输出。
- [LTR_retriever 发布说明](https://github.com/oushujun/LTR_retriever/releases)：13 字符限制的说明。
- [Slurm job array 文档](https://slurm.schedmd.com/job_array.html)：`%50` 同时运行任务限制。
- [Slurm squeue 文档](https://slurm.schedmd.com/squeue.html)：使用 `%F` 查询数组父作业 ID。
