#!/bin/bash
### 脚本名称：runtime.sh
### 功能：集中保存可覆盖的软件命令和默认统计参数。
### 输入：可覆盖环境变量；输出：主控脚本运行环境。
### 预期运行路径：流程根目录；环境名称：poaceae_hsf_core、rnaseq_de。
### 使用示例：由main.sh自动读取；路径可在此处修改。
export PIPE_PYTHON="${PIPE_PYTHON:-python3}"
export PIPE_RSCRIPT="${PIPE_RSCRIPT:-Rscript}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MIN_COUNT="${MIN_COUNT:-10}" MIN_SAMPLES="${MIN_SAMPLES:-2}"
export FDR="${FDR:-0.05}" ABS_LOG2FC="${ABS_LOG2FC:-1}"
export PARTITION="${PARTITION:-x86_64_4cpu}" THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-2}}"
