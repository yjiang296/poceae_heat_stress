# FIRE TSS 批量绘图

只需要三个文件：`plot_fire_tss.sh`、本说明、`sample_info_config.txt`。

## 1. 填写样本信息表

每行四列，Tab 或空白分隔，表头以 `#` 开头：

```text
#sample_name  bigwig_path          tss_bed_path          output_dir
sample1       /data/sample1.bw     /ref/sample1.tss.bed  /results/sample1
sample2       /data/sample2.bw     /ref/sample2.tss.bed  /results/sample2
```

样本名唯一；输入文件必须存在，输出目录自动创建。建议使用绝对路径，路径不能含空白。相对路径以执行命令时的目录为基准。示例使用占位路径，运行前请填写自己的真实输入路径。

## 2. 运行

```bash
# 先检查路径并生成作业脚本，不提交。
bash plot_fire_tss.sh sample_info_config.txt --no-submit

# 正式提交，每个样本一个独立 Slurm 作业。
bash plot_fire_tss.sh sample_info_config.txt
```

默认：`x86_64_4cpu` 分区，每样本 2 CPU，TSS 上下游各 2000 bp，bin 为 20 bp，不设置 Slurm 内存限制。可通过 `--partition`、`--threads`、`--before`、`--after`、`--bin-size` 调整。

环境需要 deepTools、Python3/pyBigWig 和 Slurm 已加入 PATH。脚本内已限制数值库线程；computeMatrix 使用指定 CPU 数。

## 3. TSS 检查与缺失 contig

每个计算作业开始时检查 BED6、链方向、染色体和坐标。原始 BED 不修改。默认遇到 bigWig 中不存在的 contig 会报错；若决定排除对应 TSS，则运行：

```bash
bash plot_fire_tss.sh sample_info_config.txt --exclude-missing-contigs
```

排除行与数量会保存到该样本输出目录，需检查其是否符合预期。

`--no-submit` 只检查路径和生成脚本，BED 内容检查在计算作业中执行。使用 BED 前六列，不保留 BED12 或分组标签；保留 `--missingDataAsZero`，已存在 contig 内未记录信号按零处理。

## 4. 统一绘图尺度

使用同一组 `--y-min/--y-max` 统一平均曲线及热图顶部曲线，`--z-min/--z-max` 统一热图色标。例如：

```bash
# 10、50 仅示范语法，请根据实际信号范围修改。
bash plot_fire_tss.sh sample_info_config.txt \
  --exclude-missing-contigs --y-min 0 --y-max 10 --z-min 0 --z-max 50
```

默认下限为 0；不设置上限时仍是各样本自动缩放。已有矩阵后，可加 `--plot-only` 重绘，避免重新计算矩阵；应等上一批作业结束后再重绘。

图中是平均 FIRE coverage，不额外进行测序深度归一化，也不是背景富集倍数。

## 5. 输出

所有结果写入样本表第四列目录，文件以样本名为前缀：

- `.plot_fire_tss.slurm/.out/.err`：作业脚本和日志。
- `.TSS.used.bed/.excluded.bed/.input_qc.tsv`：实际使用的 TSS、缺失 contig 排除行和计数。
- `.TSS.matrix.gz/.matrix.tsv/.sorted.bed`：信号矩阵和排序后的区域。
- `.TSS.profile.pdf/.png/.tsv`：平均曲线及 source data。
- `.TSS.heatmap.pdf/.png`：热图。
- `.TSS.run.log/.version.log`：运行及软件版本记录。

提交命令会输出作业号。不要对同一配置重复执行提交命令；重复运行会再次提交，并可能覆盖同名结果。

## 6. 多样本汇总

单样本作业完成后，可以用同一张样本表运行 `summarize_tss.py`，将各样本的峰值/侧翼背景比汇总为散点叠加箱线图，详见 [汇总说明](summarize_tss_usage.md)。汇总脚本要求使用默认的 ±2 kb 窗口和 20 bp 分箱。
