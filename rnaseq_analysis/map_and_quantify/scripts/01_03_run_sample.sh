#!/bin/bash
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -J rnaseq_sample
#SBATCH -o ./logs/rnaseq_sample_%j.out
#SBATCH -e ./logs/rnaseq_sample_%j.err

### 脚本名称：01_03_run_sample.sh
### 功能：依次完成单样本 fastp、HISAT2 比对和 StringTie 参考注释定量
### 输入文件：样本名、双端 FASTQ、HISAT2 索引前缀、GTF、参考基因组名称
### 输出文件：清洗 FASTQ、排序 BAM、比对统计、StringTie GTF 和基因丰度表
### 预期运行路径：由 main.sh 提交
### 环境名称：rnaseq
### 环境激活命令：conda activate rnaseq
### 使用示例：由 main.sh 自动调用，不建议手动拼接参数

set -euo pipefail

if [[ "$#" -ne 8 ]]; then
    echo "错误：本脚本需要 8 个参数，应通过 main.sh 调用" >&2
    exit 1
fi

SAMPLE="$1"
R1="$2"
R2="$3"
INDEX_PREFIX="$4"
GTF="$5"
REFERENCE="$6"
RESULTS_ROOT="$7"
STEPS="$8"
set --

THREADS="${SLURM_CPUS_PER_TASK:-20}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLE_DIR="${RESULTS_ROOT}/samples/${SAMPLE}"
FASTP_DIR="${SAMPLE_DIR}/01_fastp"
ALIGN_DIR="${SAMPLE_DIR}/02_hisat2"
STRINGTIE_DIR="${SAMPLE_DIR}/03_stringtie"
mkdir -p "${SAMPLE_DIR}"
RUN_LOG="${SAMPLE_DIR}/run.log"

activate_env() {
    if command -v mamba >/dev/null 2>&1; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate rnaseq
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate rnaseq
    else
        echo "错误：找不到 mamba 或 conda，请先在批处理环境中初始化 Conda" >&2
        exit 1
    fi
}

log_time() {
    echo "[$(date '+%F %T')] [${SAMPLE}] $*" | tee -a "${RUN_LOG}"
}

has_step() {
    [[ ",${STEPS}," == *",$1,"* ]]
}

require_file() {
    if [[ ! -s "$1" ]]; then
        echo "错误：缺少或为空：$1；请先运行步骤 $2 或 all" | tee -a "${RUN_LOG}" >&2
        exit 1
    fi
}

activate_env
{
    echo "运行目录：$(pwd)"
    echo "Conda 环境：${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
    echo "样本：${SAMPLE}"
    echo "参考基因组：${REFERENCE}"
    echo "执行步骤：${STEPS}"
    echo "线程数：${THREADS}"
} | tee "${RUN_LOG}"

for INPUT_FILE in "${R1}" "${R2}" "${GTF}"; do
    require_file "${INPUT_FILE}" "输入准备"
done

CLEAN_R1="${FASTP_DIR}/${SAMPLE}.clean.R1.fastq.gz"
CLEAN_R2="${FASTP_DIR}/${SAMPLE}.clean.R2.fastq.gz"
SORTED_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
UNSORTED_BAM="${ALIGN_DIR}/${SAMPLE}.unsorted.bam"
HISAT2_SUMMARY="${ALIGN_DIR}/${SAMPLE}.hisat2_summary.txt"
SPLICE_SITES="${ALIGN_DIR}/${SAMPLE}.known_splice_sites.txt"
STRINGTIE_GTF="${STRINGTIE_DIR}/${SAMPLE}.stringtie.gtf"
GENE_ABUND="${STRINGTIE_DIR}/${SAMPLE}.gene_abund.tsv"

if has_step "01"; then
    mkdir -p "${FASTP_DIR}"
    log_time "Step 01 fastp 开始"
    rm -f "${CLEAN_R1}" "${CLEAN_R2}" \
        "${FASTP_DIR}/${SAMPLE}.fastp.json" "${FASTP_DIR}/${SAMPLE}.fastp.html"
    fastp \
        --thread "${THREADS}" \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${CLEAN_R1}" \
        --out2 "${CLEAN_R2}" \
        --json "${FASTP_DIR}/${SAMPLE}.fastp.json" \
        --html "${FASTP_DIR}/${SAMPLE}.fastp.html" \
        2>&1 | tee "${FASTP_DIR}/${SAMPLE}.fastp.log"
    require_file "${CLEAN_R1}" "01"
    require_file "${CLEAN_R2}" "01"
    require_file "${FASTP_DIR}/${SAMPLE}.fastp.json" "01"
    fastp --version > "${FASTP_DIR}/version.log" 2>&1
    log_time "Step 01 fastp 完成"
fi

if has_step "02"; then
    require_file "${CLEAN_R1}" "01"
    require_file "${CLEAN_R2}" "01"
    mkdir -p "${ALIGN_DIR}"
    HISAT2_INDEX_OPTIONS=()
    INDEX_SUFFIX="ht2"
    if [[ -s "${INDEX_PREFIX}.1.ht2l" ]]; then
        HISAT2_INDEX_OPTIONS=(--large-index); INDEX_SUFFIX="ht2l"
    fi
    for INDEX_PART in {1..8}; do require_file "${INDEX_PREFIX}.${INDEX_PART}.${INDEX_SUFFIX}" "HISAT2 索引构建"; done
    log_time "HISAT2 索引类型：${INDEX_SUFFIX}"
    log_time "Step 02 HISAT2 比对开始"
    if ! command -v hisat2_extract_splice_sites.py >/dev/null 2>&1; then
        echo "错误：找不到 hisat2_extract_splice_sites.py" | tee -a "${RUN_LOG}" >&2
        exit 1
    fi
    hisat2_extract_splice_sites.py "${GTF}" > "${SPLICE_SITES}"
    rm -f "${UNSORTED_BAM}" "${SORTED_BAM}" "${SORTED_BAM}.bai" "${SORTED_BAM}.csi" "${HISAT2_SUMMARY}"
    trap 'rm -f "${UNSORTED_BAM}"' EXIT
    hisat2 \
        "${HISAT2_INDEX_OPTIONS[@]}" \
        --threads "${THREADS}" \
        --dta \
        --known-splicesite-infile "${SPLICE_SITES}" \
        -x "${INDEX_PREFIX}" \
        -1 "${CLEAN_R1}" \
        -2 "${CLEAN_R2}" \
        --summary-file "${HISAT2_SUMMARY}" \
        2> "${ALIGN_DIR}/${SAMPLE}.hisat2.log" |
        samtools view -b -o "${UNSORTED_BAM}" -
    require_file "${UNSORTED_BAM}" "02"
    samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" "${UNSORTED_BAM}"
    rm -f "${UNSORTED_BAM}"
    bash "${SCRIPT_DIR}/02_index_bam.sh" "${SORTED_BAM}" "${ALIGN_DIR}"
    samtools flagstat -@ "${THREADS}" "${SORTED_BAM}" \
        > "${ALIGN_DIR}/${SAMPLE}.flagstat.txt"
    require_file "${SORTED_BAM}" "02"
    require_file "${HISAT2_SUMMARY}" "02"
    {
        hisat2 --version | head -n 1
        samtools --version | head -n 2
    } > "${ALIGN_DIR}/version.log"
    log_time "Step 02 HISAT2 比对完成"
fi

if has_step "03"; then
    require_file "${SORTED_BAM}" "02"
    mkdir -p "${STRINGTIE_DIR}"
    log_time "Step 03 StringTie 定量开始"
    rm -f "${STRINGTIE_GTF}" "${GENE_ABUND}"
    stringtie "${SORTED_BAM}" \
        -p "${THREADS}" \
        -G "${GTF}" \
        -e \
        -B \
        -o "${STRINGTIE_GTF}" \
        -A "${GENE_ABUND}" \
        2>&1 | tee "${STRINGTIE_DIR}/${SAMPLE}.stringtie.log"
    require_file "${STRINGTIE_GTF}" "03"
    require_file "${GENE_ABUND}" "03"
    stringtie --version > "${STRINGTIE_DIR}/version.log" 2>&1
    log_time "Step 03 StringTie 定量完成"
fi

log_time "所选样本步骤全部完成"
