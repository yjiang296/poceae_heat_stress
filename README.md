# poceae_heat_stress

禾本科比较基因组研究的可复用分析流程。

## StringTie表达定量与差异表达

[rnaseq_analysis](rnaseq_analysis/README.md)：从双端RNA-seq的fastp、HISAT2、StringTie固定注释定量，到同来源counts/TPM矩阵和DESeq2全基因差异分析。保留逐样本定量、逐材料DEG作业、断点重跑、配置填写指南与模拟测试。

## 蛋白功能注释

[protein_function_annotation](protein_function_annotation/README.md)：对已有基因模型对应的蛋白 FASTA 运行 InterProScan 5 和 eggNOG-mapper，按材料生成独立 Slurm 作业，支持输入检查、分块续跑和原蛋白 ID 恢复。

仓库保存代码、测试、配置模板和说明。实际软件环境、数据库、输入序列、服务器配置和计算结果保存在运行机器上。

## FIRE TSS 绘图与多样本汇总

[FIRE_TSS 单样本绘图](FIRE_TSS/plot_fire_tss_usage.md)：按四列样本表为每个样本提交独立 Slurm 作业，生成 TSS 平均 FIRE coverage 曲线和热图，支持缺失 contig 检查及统一绘图尺度。

[FIRE_TSS 多样本汇总](FIRE_TSS/summarize_tss_usage.md)：汇总平均曲线的峰值/侧翼背景比，绘制散点叠加箱线图，保留绘图表格及可复现快照，输出文件按实际样本数命名。两步共用 [示例样本表](FIRE_TSS/sample_info_config.txt)。

## LAI 基因组组装质量评估

[LAI 批量评估](LAI/README.md)：统一短 contig ID，按样本执行 LTR_FINDER、LTRharvest 和 LTR_retriever，通过 Slurm 数组限制最多 50 个样本并行，并在验证成功后定向清理大型中间文件。附 [参考样本表](LAI/sample_info.txt)、使用说明和功能测试。

## 一条命令更新 GitHub

首次使用需要 Git、Python 3，以及本机可用的 GitHub Git 写入认证（SSH key 或 HTTPS 凭据管理器）。Codex 的 GitHub 连接授权不会自动变成本机 Git 凭据。不要把令牌写入脚本或仓库。

```bash
git clone https://github.com/yjiang296/poceae_heat_stress.git
cd poceae_heat_stress
# 修改本克隆中的代码或说明后：
bash update_github.sh "说明本次修改"
```

脚本在当前克隆的 `main` 分支运行，检查 origin 地址、远端进度、Shell 语法和回归测试，然后只暂存 `publish_files.txt` 列出的文件、提交并推送。新增源文件需要先将其加入清单。不会自动复制超算或其他目录中的改动，不会自动合并远端更新，也不会强制推送。

如果本机尚未配置 Git 写入认证，可继续让 Codex 通过已连接的 GitHub 工具上传。脚本中断于推送时，本地提交仍保留；解决认证后再次运行即可重试。Git 提交身份使用本机已有的 `user.name` 与 `user.email`。

同步运行中的工作流时请保留原代码和配置：已完成阶段使用代码/配置指纹校验，更新后应另建输出目录。
