#!/bin/bash
#SBATCH -J stringtie_deg
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o ./logs/stringtie_deg.%j.out
#SBATCH -e ./logs/stringtie_deg.%j.err
### 脚本名称：main.sh
### 功能：独立运行一个材料的审核、矩阵整理、DESeq2与绘图，支持按步骤重跑。
### 输入：$1材料配置目录；$2输出目录；输出：四步骤结果与独立日志。
### 预期运行路径：新流程根目录；环境名称：poaceae_hsf_core、rnaseq_de。
### 使用示例：bash main.sh config/materials/B73 results/B73 all；末尾可换03或--from 03。
set -euo pipefail
[[ $# -ge 2 ]] || { echo "用法：bash main.sh 输入配置目录 输出目录 [all|01|02|03|04|--from NN]" >&2; exit 2; }
INPUT="$1"; OUTDIR="$2"; shift 2
STEPS=(01 02 03 04)
if [[ $# -eq 0 || "${1:-}" == all ]]; then
    [[ $# -le 1 ]] || exit 2
elif [[ "$1" == --from && $# -eq 2 ]]; then
    case "$2" in 01) ;; 02) STEPS=(02 03 04);; 03) STEPS=(03 04);; 04) STEPS=(04);; *) exit 2;; esac
elif [[ $# -eq 1 && "$1" =~ ^0[1-4]$ ]]; then
    STEPS=("$1")
else
    echo "步骤参数错误" >&2; exit 2
fi
### 先保存并清空位置参数，避免环境配置读取步骤编号。
set --
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$ROOT/config/runtime.sh" ]] || ROOT="${SLURM_SUBMIT_DIR:?无法定位流程根目录}"
[[ -d "$INPUT" ]] || { echo "输入目录不存在：$INPUT" >&2; exit 2; }
source "$ROOT/config/runtime.sh"
mkdir -p "$OUTDIR/logs"
INPUT="$(cd "$INPUT" && pwd)"; OUTDIR="$(cd "$OUTDIR" && pwd)"
exec 9>"$OUTDIR/.run.lock"
flock -n 9 || { echo "同一输出目录已有作业运行" >&2; exit 2; }
STAMP="$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}_$$"
log_time() { echo "[$(date '+%F %T')] 材料=$(basename "$INPUT") $*" | tee -a "$OUTDIR/run.log"; }
log_time "运行目录=$(pwd)；Conda=${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
STEP=初始化
trap 'code=$?; log_time "步骤${STEP}失败，退出码=${code}"; exit "$code"' ERR
"$PIPE_PYTHON" --version >> "$OUTDIR/version.log" 2>&1
"$PIPE_RSCRIPT" --version >> "$OUTDIR/version.log" 2>&1
NAMES=(01_audit 02_matrices 03_deseq2 04_plots)
for STEP in "${STEPS[@]}"; do
    index=$((10#$STEP-1)); DEST="$OUTDIR/${NAMES[$index]}"
    UPSTREAM="$INPUT"
    if (( index>0 )); then
        UPSTREAM="$OUTDIR/${NAMES[$((index-1))]}"
        prior="${NAMES[$((index-1))]:0:2}"
        [[ -s "$UPSTREAM/${prior}_SUCCESS" ]] || { echo "缺少步骤${prior}结果，请先运行${prior}或all" >&2; exit 2; }
    fi
    ### 重跑上游时使旧下游成功标记失效，防止误用旧结果；结果文件仍保留以便审计。
    for ((j=index; j<4; j++)); do
        rm -f "$OUTDIR/${NAMES[$j]}/${NAMES[$j]:0:2}_SUCCESS"
    done
    log_time "步骤${STEP}开始"; STARTED=$SECONDS
    if [[ "$STEP" == 01 || "$STEP" == 02 ]]; then
        CONDA_DEFAULT_ENV=poaceae_hsf_core "$PIPE_PYTHON" "$ROOT/scripts/${NAMES[$index]}.py" --input "$UPSTREAM" --outdir "$DEST" 2>&1 | tee "$OUTDIR/logs/${STEP}_${STAMP}.log"
    else
        if [[ "$STEP" == 04 ]]; then
            export FONTCONFIG_FILE="$ROOT/assets/fonts.conf"
            fc-match Arial | grep -qi Arial || { echo "Arial字体不可用" >&2; exit 2; }
        fi
        CONDA_DEFAULT_ENV=rnaseq_de "$PIPE_RSCRIPT" "$ROOT/scripts/${NAMES[$index]}.R" --input "$UPSTREAM" --outdir "$DEST" --min-count "$MIN_COUNT" --min-samples "$MIN_SAMPLES" --fdr "$FDR" --lfc "$ABS_LOG2FC" 2>&1 | tee "$OUTDIR/logs/${STEP}_${STAMP}.log"
    fi
    [[ -s "$DEST/${STEP}_SUCCESS" ]] || { echo "步骤${STEP}未生成成功标记" >&2; exit 1; }
    log_time "步骤${STEP}结束，耗时$((SECONDS-STARTED))秒"
done
