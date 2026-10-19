#!/bin/bash
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -J rnaseq_global
#SBATCH -o ./logs/rnaseq_global_%j.out
#SBATCH -e ./logs/rnaseq_global_%j.err

### 脚本名称：04_07_run_global.sh
### 功能：运行统计汇总、StringTie count 矩阵和 TPM/FPKM 矩阵构建
### 输入文件：样本信息表和样本级结果目录
### 输出文件：统计 TSV、MultiQC 报告以及按参考基因组分组的表达矩阵
### 预期运行路径：由 main.sh 提交
### 环境名称：rnaseq
### 环境激活命令：conda activate rnaseq
### 使用示例：由 main.sh 自动调用，不建议手动拼接参数

set -euo pipefail

if [[ "$#" -ne 3 ]]; then
    echo "错误：本脚本需要 sample_info、outdir 和步骤列表" >&2
    exit 1
fi

SAMPLE_INFO="$1"
RESULTS_ROOT="$2"
STEPS="$3"
set --

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${RESULTS_ROOT}/summary"
mkdir -p "${OUTDIR}"
RUN_LOG="${OUTDIR}/run.log"

activate_env() {
    if command -v mamba >/dev/null 2>&1; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate rnaseq
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate rnaseq
    else
        echo "错误：找不到 mamba 或 conda" >&2
        exit 1
    fi
}

has_step() {
    [[ ",${STEPS}," == *",$1,"* ]]
}

run_step() {
    local STEP="$1"
    shift
    echo "[$(date '+%F %T')] ${STEP} 开始" | tee -a "${RUN_LOG}"
    "$@" 2>&1 | tee "${OUTDIR}/${STEP}.log"
    echo "[$(date '+%F %T')] ${STEP} 完成" | tee -a "${RUN_LOG}"
}

activate_env
PYTHON_VERSION_OK="$(
    python -c 'import sys; print(int(sys.version_info >= (3, 8)))'
)"
if [[ "${PYTHON_VERSION_OK}" != "1" ]]; then
    echo "错误：汇总脚本要求 Python >= 3.8，当前版本为 $(python --version 2>&1)" >&2
    echo "请先升级 rnaseq 环境中的 Python，然后重新运行 --from 04" >&2
    exit 1
fi
{
    echo "运行目录：$(pwd)"
    echo "Conda 环境：${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
    echo "执行步骤：${STEPS}"
} | tee "${RUN_LOG}"

if has_step "04"; then
    run_step "04_collect_stats" \
        python "${SCRIPT_DIR}/04_collect_stats.py" \
        --input "${SAMPLE_INFO}" \
        --outdir "${OUTDIR}/04_stats" \
        --results-root "${RESULTS_ROOT}"
    mkdir -p "${OUTDIR}/04_stats/multiqc"
    multiqc "${RESULTS_ROOT}/samples" \
        --outdir "${OUTDIR}/04_stats/multiqc" \
        --filename "04_rnaseq_multiqc_report.html" \
        --force \
        2>&1 | tee "${OUTDIR}/04_stats/multiqc.log"
fi

if has_step "05"; then
    run_step "05_run_prepde" \
        python "${SCRIPT_DIR}/05_run_prepde.py" \
        --input "${SAMPLE_INFO}" \
        --outdir "${OUTDIR}/05_counts" \
        --results-root "${RESULTS_ROOT}"
fi

if has_step "06"; then
    run_step "06_build_gene_matrices" \
        python "${SCRIPT_DIR}/06_build_gene_matrices.py" \
        --input "${SAMPLE_INFO}" \
        --outdir "${OUTDIR}/06_abundance" \
        --results-root "${RESULTS_ROOT}"
fi

if has_step "07"; then
    run_step "07_build_transcript_matrices" \
        python "${SCRIPT_DIR}/07_build_transcript_matrices.py" \
        --input "${SAMPLE_INFO}" \
        --outdir "${OUTDIR}/06_abundance" \
        --results-root "${RESULTS_ROOT}"
fi

{
    fastp --version
    hisat2 --version | head -n 1
    samtools --version | head -n 2
    stringtie --version
    multiqc --version
    python --version
} > "${OUTDIR}/version.log" 2>&1
