#!/bin/bash
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -J index_bam
#SBATCH -o ./logs/index_bam_%j.out
#SBATCH -e ./logs/index_bam_%j.err

### 脚本名称：02_index_bam.sh
### 功能：根据 BAM header 中最长参考序列自动生成 BAI 或 CSI 索引
### 输入文件：坐标排序后的 BAM
### 输出文件：BAM 旁边的 .bai 或 .csi 索引及索引类型记录
### 预期运行路径：任意目录或由 01_03_run_sample.sh 调用
### 环境名称：rnaseq
### 环境激活命令：conda activate rnaseq
### 使用示例：bash scripts/02_index_bam.sh sample.sorted.bam sample_output_dir

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "用法：bash $0 <sorted.bam> <outdir>" >&2
    exit 1
fi

INPUT_BAM="$1"
OUTDIR="$2"
THREADS="${SLURM_CPUS_PER_TASK:-2}"
BAI_MAX_REFERENCE_LENGTH=536870912

if [[ ! -f "${INPUT_BAM}" || ! -s "${INPUT_BAM}" ]]; then
    echo "错误：输入 BAM 不存在或为空：${INPUT_BAM}" >&2
    exit 1
fi
if ! samtools quickcheck "${INPUT_BAM}"; then
    echo "错误：BAM 未完整或格式异常：${INPUT_BAM}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
RUN_LOG="${OUTDIR}/02_bam_index_run.log"
{
    echo "运行目录：$(pwd)"
    echo "Conda 环境：${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
    echo "输入 BAM：${INPUT_BAM}"
    echo "线程数：${THREADS}"
} | tee "${RUN_LOG}"

MAX_REFERENCE_LENGTH="$(
    samtools view -H "${INPUT_BAM}" |
        awk 'BEGIN{max=0} $1=="@SQ"{for(i=1;i<=NF;i++){if($i~/^LN:/){v=substr($i,4)+0;if(v>max)max=v}}} END{print max}'
)"

if [[ ! "${MAX_REFERENCE_LENGTH}" =~ ^[0-9]+$ || "${MAX_REFERENCE_LENGTH}" -eq 0 ]]; then
    echo "错误：无法从 BAM header 读取参考序列长度" | tee -a "${RUN_LOG}" >&2
    exit 1
fi

rm -f "${INPUT_BAM}.bai" "${INPUT_BAM}.csi"
if (( MAX_REFERENCE_LENGTH >= BAI_MAX_REFERENCE_LENGTH )); then
    INDEX_TYPE="CSI"
    INDEX_PATH="${INPUT_BAM}.csi"
    samtools index -@ "${THREADS}" -c "${INPUT_BAM}"
else
    INDEX_TYPE="BAI"
    INDEX_PATH="${INPUT_BAM}.bai"
    samtools index -@ "${THREADS}" -b "${INPUT_BAM}"
fi

if [[ ! -s "${INDEX_PATH}" ]]; then
    echo "错误：${INDEX_TYPE} 索引未正常生成：${INDEX_PATH}" | tee -a "${RUN_LOG}" >&2
    exit 1
fi

BAM_BASENAME="$(basename "${INPUT_BAM}" .bam)"
{
    echo -e "bam\t${INPUT_BAM}"
    echo -e "max_reference_length\t${MAX_REFERENCE_LENGTH}"
    echo -e "index_type\t${INDEX_TYPE}"
    echo -e "index_path\t${INDEX_PATH}"
} > "${OUTDIR}/${BAM_BASENAME}.index_type.txt"

{
    samtools --version | head -n 2
} > "${OUTDIR}/02_bam_index_version.log"
echo "BAM 索引完成：${INDEX_TYPE}，最长参考序列 ${MAX_REFERENCE_LENGTH} bp" | tee -a "${RUN_LOG}"
