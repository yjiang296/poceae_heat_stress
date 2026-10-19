# TSS 峰值富集汇总

输入为四列样本配置表（样本名、bigWig、TSS BED、结果目录，Tab 分隔），样本数量不限。脚本直接读取每个结果目录的 `.TSS.profile.tsv` 和 `.TSS.matrix.gz` 首行参数。要求已有结果的窗口为 TSS 上下游各 2 kb、20 bp 分箱，缺失信号按 0 处理。

## 计算口径

原曲线是各样本所有纳入 TSS 的 **平均 FIRE coverage**，原始纵轴没有背景归一化。这里将“最高富集倍数”定义为：

`peak_enrichment_fold = max_mean_FIRE_coverage / flank_mean_FIRE_coverage`

- 峰值：平均曲线在 `[-2000, 2000)` bp 内 200 个 bin 的最大值，无额外平滑；不是单个基因的最大值，也不是各基因峰值的平均值。
- 背景：同一平均曲线 `[-2000, -1500)` 和 `[1500, 2000)` bp 两段共 50 个 bin 的算术平均值。
- 坐标相对于按链方向对齐的 TSS；区间为左闭右开。若并列最高，记录最左侧 bin。
- 这是本次定义的 FIRE 峰值/侧翼背景比，不应与其他软件定义的 ATAC-seq TSS enrichment score 直接比较；侧翼背景也不等同于全基因组背景。
- 倍数降低了整体信号量级差异，但仍可能受覆盖度、TSS 注释及背景区域生物学差异影响。不同材料并不因此获得完全相同的生物学基线。

## 交付文件

绘图脚本为 `summarize_tss.py`；结果保存在 `--outdir` 指定的目录。文件名中的 `N` 自动取实际有效样本数，不从目录名猜测：6 个样本用 `6`，142 个样本用 `142`。

- `tss_peak_enrichment_N.png` / `.pdf`：散点叠加箱线图，每个点对应一个样本，使用轻微横向抖动及半透明显示。箱体为 Q1–Q3，中线为中位数，须线到 1.5 IQR 范围内最远观测值；所有样本均显示，包括离群点。
- `tss_peak_summary.N.tsv`：N 行，每个样本一行，保持配置顺序。包含原始峰值、背景均值、峰值倍数、峰值区间、纳入 TSS 数、倍数降序排名及源文件 SHA256。样本名须按文本导入，以保留前导 0。
- `tss_profiles.N.long.tsv`：N × 200 行，保存所有样本每个 bin 的原始平均 coverage、归一化倍数及背景标记。
- `source_profiles.json`：所有输入曲线原文、SHA256、矩阵参数和原配置的轻量快照，可离线复现，无需下载矩阵主体或 bigWig。
- `sample_info_config.N.txt`：本次输入配置副本。`n_TSS` 读取实际结果矩阵；配置中的原始 BED 路径不代表此前过滤后的实际保留行数。
- `summarize_tss.py`：汇总和绘图代码。依赖 Python 3、NumPy、Matplotlib。

## 运行

先按 [单样本绘图说明](plot_fire_tss_usage.md) 完成每个样本的曲线计算，再在脚本目录执行：

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python3 summarize_tss.py --config sample_info_config.txt --outdir summary
```

`sample_info_config.txt` 使用与单样本脚本相同的四列配置，建议使用 Tab 分隔。脚本只汇总配置中列出的样本，不自动扫描目录。输入含 6 个样本时输出文件带 `6`，含 142 个样本时带 `142`。

使用已保存快照重新绘制（无需访问原始结果路径）：

```bash
python3 summarize_tss.py --snapshot summary/source_profiles.json --outdir summary_replot
```

重新执行会更新输出目录中的同名汇总文件。脚本会检查样本唯一性、矩阵窗口/分箱一致性、曲线完整性及背景为正；遇到异常即退出，不会静默丢弃样本。

建议不同样本集合使用不同输出目录。更改样本数时脚本不会自动删除该目录中历史批次的文件。随机抖动使用固定种子 142，仅用于保证散点位置可复现，不限制样本数量。
