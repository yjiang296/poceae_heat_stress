# LAI 批量评估

使用 LTR_FINDER_parallel、GenomeTools LTRharvest 和 LTR_retriever 计算基因组 LAI。流程先生成唯一短 contig ID，按样本提交 Slurm 数组任务，最多同时运行 50 个样本，并在各阶段成功后定向清理大型中间文件。

- [完整使用说明、B73 诊断及清理清单](submit_ltr_pipeline.md)
- [主脚本](submit_ltr_pipeline.sh)
- [参考样本表](sample_info.txt)：三样本示例，路径为占位符，运行前必须修改。
- [功能测试](test_lai_t6.py)
- [当前本地测试记录](validation.txt)
- [此前服务器 Python 3.6 验证摘要](remote_validation.txt)：不包含服务器账号、绝对路径和原始日志。

```bash
# 先编辑 sample_info.txt，将占位路径改为真实路径
# 在仓库的 LAI 目录中执行；不会提交作业
bash submit_ltr_pipeline.sh sample_info.txt --no-submit

# 使用小型模拟数据测试，不调用真实生物信息学工具
python3 test_lai_t6.py

# 实际提交示例样本；可改为完整样本表
bash submit_ltr_pipeline.sh sample_info.txt --max-running 50 --submit
```

默认 FINDER 路径为 `~/software/LTR_FINDER_parallel/LTR_FINDER_parallel`，conda 初始化脚本为 `~/miniconda3/etc/profile.d/conda.sh`；环境名与分区沿用原流程。迁移前通过 `--finder`、`--conda-sh`、`--env`、`--partition` 调整，详见主脚本 `--help`。

完整基因组数据和分析结果不包含在仓库内。功能测试及 210 样本生成预检通过，不代表已经完成全部真实样本计算。
