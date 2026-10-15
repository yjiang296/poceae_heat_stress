#!/bin/bash
#SBATCH -J protein_function
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o protein_function.%j.out
#SBATCH -e protein_function.%j.err
### 脚本名称：main.sh；功能：单材料四步功能注释；输入：$1 蛋白 FASTA；输出：$2 材料目录。
### 预期路径：任意；环境：配置文件指定的独立工具环境；示例：bash main.sh a.fa results/A all --sample A
### 批量任务请用 submit_annotation.sh；不在此循环多个物种；无随机步骤，不设置内存限制。
set -euo pipefail
[[ $# -ge 2 ]] || { echo '用法：main.sh 蛋白FASTA 输出目录 [all|01|02|03|04] --sample 名称 [--from 03] [--config 文件] [--threads 8]' >&2; exit 2; }
INPUT=$1; OUTDIR=$2; shift 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONFIG="$ROOT/config.json"; SAMPLE=''; STEPS=all; THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-2}}
if [[ $# -gt 0 && "$1" != --* ]]; then STEPS=$1; shift; fi
while [[ $# -gt 0 ]]; do
 case "$1" in
 --sample) SAMPLE=$2; shift 2;;
 --config) CONFIG=$2; shift 2;;
 --threads) THREADS=$2; shift 2;;
 --from) case "$2" in 01) STEPS=all;; 02) STEPS=02,03,04;; 03) STEPS=03,04;; 04) STEPS=04;; *) exit 2;; esac; shift 2;;
 *) echo "未知参数：$1" >&2; exit 2;;
 esac
done
[[ -f "$INPUT" && -f "$CONFIG" ]] || { echo '输入 FASTA 或配置文件不存在' >&2; exit 2; }
[[ "$SAMPLE" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo '样本名或线程数非法' >&2; exit 2; }
case "$STEPS" in all) STEPS=01,02,03,04;; esac
IFS=',' read -r -a CHOSEN <<< "$STEPS"
for STEP in "${CHOSEN[@]}"; do [[ "$STEP" =~ ^0[1-4]$ ]] || { echo '步骤必须为 01 至 04' >&2; exit 2; }; done
PYTHON=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["python"])' "$CONFIG")
[[ -x "$PYTHON" ]] || { echo "Python 不可执行：$PYTHON" >&2; exit 2; }
INPUT=$($PYTHON -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$INPUT")
CONFIG=$($PYTHON -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$CONFIG")
mkdir -p "$OUTDIR"; OUTDIR=$(cd "$OUTDIR" && pwd)
exec 9>"$OUTDIR/.workflow.lock"
flock -n 9 || { echo "该材料已有任务运行：$OUTDIR" >&2; exit 1; }
exec > >(tee -a "$OUTDIR/run.log") 2>&1
log_time() { printf '%s [%s] %s\n' "$(date '+%F %T')" "$SAMPLE" "$*"; }
STEP=初始化
trap 'log_time "失败：步骤 $STEP，行 $LINENO"' ERR
log_time "工作目录 $(pwd)；环境 ${CONDA_DEFAULT_ENV:-未检测到 conda 环境}；线程 $THREADS；步骤 $STEPS"
### 所有参数在工具环境设置前保存；工具通过独立 PATH 启动，不调用旧式 conda 激活。
### 单步/中途启动也核对原始输入、配置与 01 指纹，避免旧目录混入新输入。
if [[ ",$STEPS," != *,01,* ]]; then
 "$PYTHON" - "$ROOT/scripts" "$INPUT" "$OUTDIR/01" "$SAMPLE" "$CONFIG" "$THREADS" <<'PY_CHECK'
import json,sys
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0,sys.argv[1])
from common import signature,cached
try:
 a=SimpleNamespace(sample=sys.argv[4],threads=int(sys.argv[6]),cfg=json.loads(Path(sys.argv[5]).read_text()))
 marker=Path(sys.argv[3])/'01_complete.json'
 if not cached(marker,signature(a,[Path(sys.argv[2])])):raise ValueError('缺少 01 完成记录，请先运行 01 或 all')
except Exception as e:
 print('上游输入检查失败：'+str(e),file=sys.stderr);sys.exit(1)
PY_CHECK
fi
for STEP in "${CHOSEN[@]}"; do
 START=$SECONDS; log_time "开始 $STEP"
 case "$STEP" in
 01) SCRIPT=01_prepare.py; UPSTREAM=$INPUT;;
 02) SCRIPT=02_interpro.py; UPSTREAM="$OUTDIR/01";;
 03) SCRIPT=03_eggnog.py; UPSTREAM="$OUTDIR/01";;
 04) SCRIPT=04_merge.py; UPSTREAM=$OUTDIR;;
 esac
 "$PYTHON" "$ROOT/scripts/$SCRIPT" --input "$UPSTREAM" --outdir "$OUTDIR/$STEP" --config "$CONFIG" --sample "$SAMPLE" --threads "$THREADS"
 log_time "结束 $STEP，耗时 $((SECONDS-START)) 秒"
done
printf '%s\n' "$BASH_VERSION" >> "$OUTDIR/version.log"
log_time '所选步骤完成'
