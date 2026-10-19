#!/bin/bash
#SBATCH -J stDE_test
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o tests/validation.%j.out
#SBATCH -e tests/validation.%j.err
### 脚本名称：run_validation.sh
### 功能：模拟端到端和错误输入拒绝测试。
### 输入：$1为模拟工作目录；$2为验证目录；输出：测试日志和通过标记。
### 预期运行路径：新流程根目录；环境名称：poaceae_hsf_core、rnaseq_de。
### 使用示例：sbatch tests/run_validation.sh tests/work tests/validation。
set -euo pipefail
[[ $# -eq 2 && -d "$1" ]] || exit 2
WORK="$1"; OUTDIR="$2"; set --
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -f "$ROOT/config/runtime.sh" ]] || ROOT="${SLURM_SUBMIT_DIR:?无法定位流程根目录}"
source "$ROOT/config/runtime.sh"
mkdir -p "$OUTDIR"
echo "运行目录=$(pwd)；Conda=${CONDA_DEFAULT_ENV:-未检测到 conda 环境}" | tee "$OUTDIR/run.log"
bash "$ROOT/main.sh" "$WORK/config" "$WORK/results" all
CONDA_DEFAULT_ENV=poaceae_hsf_core "$PIPE_PYTHON" "$ROOT/tests/02_check.py" --input "$WORK/results" --outdir "$OUTDIR/positive"
CONDA_DEFAULT_ENV=poaceae_hsf_core "$PIPE_PYTHON" "$ROOT/tests/03_negative.py" --input "$WORK" --outdir "$OUTDIR/negative"
"$PIPE_PYTHON" --version > "$OUTDIR/version.log" 2>&1
"$PIPE_RSCRIPT" --version >> "$OUTDIR/version.log" 2>&1
