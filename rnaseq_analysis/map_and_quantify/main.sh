#!/bin/bash
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -J rnaseq_submit
#SBATCH -o ./logs/rnaseq_submit.log
#SBATCH -e ./logs/rnaseq_submit.err

### 脚本名称：main.sh
### 功能：读取样本信息表，为每个样本生成独立 Slurm 作业，并在样本作业完成后提交全局汇总作业
### 输入文件：六列样本信息表
### 输出文件：样本级结果、参考基因组分组表达矩阵、统计表和 MultiQC 报告
### 预期运行路径：项目根目录
### 环境名称：rnaseq
### 环境激活命令：conda activate rnaseq
### 使用示例：bash main.sh sample_info.txt results all
### 使用示例：bash main.sh sample_info.txt results 02
### 使用示例：bash main.sh sample_info.txt results --from 03 --no-submit

set -euo pipefail

usage() {
    echo "用法：bash $0 <sample_info.txt> <outdir> <all|01|02|03|04|05|06|07|--from NN> [--no-submit]" >&2
}

if [[ "$#" -lt 3 || "$#" -gt 5 ]]; then
    usage
    exit 1
fi

SAMPLE_INFO="$1"
OUTDIR="$2"
shift 2

MODE="$1"
shift
FROM_STEP=""
if [[ "${MODE}" == "--from" ]]; then
    if [[ "$#" -lt 1 ]]; then
        usage
        exit 1
    fi
    FROM_STEP="$1"
    shift
fi

NO_SUBMIT=false
if [[ "$#" -gt 0 ]]; then
    if [[ "$#" -eq 1 && "$1" == "--no-submit" ]]; then
        NO_SUBMIT=true
    else
        usage
        exit 1
    fi
fi
set --

if [[ ! -f "${SAMPLE_INFO}" ]]; then
    echo "错误：样本信息表不存在：${SAMPLE_INFO}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
SAMPLE_INFO="$(cd "$(dirname "${SAMPLE_INFO}")" && pwd)/$(basename "${SAMPLE_INFO}")"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLE_RUNNER="${SCRIPT_DIR}/scripts/01_03_run_sample.sh"
GLOBAL_RUNNER="${SCRIPT_DIR}/scripts/04_07_run_global.sh"
VALIDATOR="${SCRIPT_DIR}/scripts/00_validate_sample_info.sh"

for REQUIRED_FILE in "${SAMPLE_RUNNER}" "${GLOBAL_RUNNER}" "${VALIDATOR}"; do
    if [[ ! -f "${REQUIRED_FILE}" ]]; then
        echo "错误：缺少流程脚本：${REQUIRED_FILE}" >&2
        exit 1
    fi
done

bash "${VALIDATOR}" "${SAMPLE_INFO}" "${OUTDIR}"

SAMPLE_STEPS=""
GLOBAL_STEPS=""
case "${MODE}:${FROM_STEP}" in
    all:) SAMPLE_STEPS="01,02,03"; GLOBAL_STEPS="04,05,06,07" ;;
    01:) SAMPLE_STEPS="01" ;;
    02:) SAMPLE_STEPS="02" ;;
    03:) SAMPLE_STEPS="03" ;;
    04:) GLOBAL_STEPS="04" ;;
    05:) GLOBAL_STEPS="05" ;;
    06:) GLOBAL_STEPS="06" ;;
    07:) GLOBAL_STEPS="07" ;;
    --from:01) SAMPLE_STEPS="01,02,03"; GLOBAL_STEPS="04,05,06,07" ;;
    --from:02) SAMPLE_STEPS="02,03"; GLOBAL_STEPS="04,05,06,07" ;;
    --from:03) SAMPLE_STEPS="03"; GLOBAL_STEPS="04,05,06,07" ;;
    --from:04) GLOBAL_STEPS="04,05,06,07" ;;
    --from:05) GLOBAL_STEPS="05,06,07" ;;
    --from:06) GLOBAL_STEPS="06,07" ;;
    --from:07) GLOBAL_STEPS="07" ;;
    *)
        echo "错误：无效步骤选择：${MODE} ${FROM_STEP}" >&2
        usage
        exit 1
        ;;
esac

mkdir -p "${OUTDIR}/samples" "${OUTDIR}/summary" "${OUTDIR}/jobs"
JOB_IDS=()

if [[ -n "${SAMPLE_STEPS}" ]]; then
    while read -r SAMPLE R1 R2 INDEX_PREFIX GTF REFERENCE; do
        [[ -z "${SAMPLE}" || "${SAMPLE}" == "sample" || "${SAMPLE}" == \#* ]] && continue
        for INPUT_FILE in "${R1}" "${R2}" "${GTF}"; do
            if [[ ! -f "${INPUT_FILE}" ]]; then
                echo "错误：${SAMPLE} 的输入文件不存在：${INPUT_FILE}" >&2
                exit 1
            fi
        done
        if [[ ! -s "${INDEX_PREFIX}.1.ht2" && ! -s "${INDEX_PREFIX}.1.ht2l" ]]; then
            echo "错误：${SAMPLE} 的 HISAT2 索引前缀无效：${INDEX_PREFIX}" >&2
            exit 1
        fi

        SAMPLE_DIR="${OUTDIR}/samples/${SAMPLE}"
        JOB_DIR="${SAMPLE_DIR}/job"
        mkdir -p "${JOB_DIR}"
        JOB_SCRIPT="${JOB_DIR}/${SAMPLE}.slurm"

        {
            echo '#!/bin/bash'
            echo '#SBATCH -p x86_64_4cpu'
            echo '#SBATCH -N 1'
            echo '#SBATCH -c 20'
            echo "#SBATCH -J rnaseq_${SAMPLE}"
            echo "#SBATCH -o ${JOB_DIR}/${SAMPLE}_%j.out"
            echo "#SBATCH -e ${JOB_DIR}/${SAMPLE}_%j.err"
            printf '\nbash %q %q %q %q %q %q %q %q %q\n' \
                "${SAMPLE_RUNNER}" "${SAMPLE}" "${R1}" "${R2}" "${INDEX_PREFIX}" \
                "${GTF}" "${REFERENCE}" "${OUTDIR}" "${SAMPLE_STEPS}"
        } > "${JOB_SCRIPT}"
        chmod +x "${JOB_SCRIPT}"

        if [[ "${NO_SUBMIT}" == true ]]; then
            printf 'sbatch %q\n' "${JOB_SCRIPT}"
        else
            JOB_ID="$(sbatch --parsable "${JOB_SCRIPT}")"
            JOB_IDS+=("${JOB_ID}")
            echo "已提交样本 ${SAMPLE}：job ${JOB_ID}"
        fi
    done < <(awk 'BEGIN{FS="[[:space:]]+"} NR>1 && $0 !~ /^[[:space:]]*(#|$)/ {print $1,$2,$3,$4,$5,$6}' "${SAMPLE_INFO}")
fi

if [[ -n "${GLOBAL_STEPS}" ]]; then
    GLOBAL_JOB="${OUTDIR}/jobs/04_07_global.slurm"
    {
        echo '#!/bin/bash'
        echo '#SBATCH -p x86_64_4cpu'
        echo '#SBATCH -N 1'
        echo '#SBATCH -c 2'
        echo '#SBATCH -J rnaseq_summary'
        echo "#SBATCH -o ${OUTDIR}/jobs/rnaseq_summary_%j.out"
        echo "#SBATCH -e ${OUTDIR}/jobs/rnaseq_summary_%j.err"
        printf '\nbash %q %q %q %q\n' "${GLOBAL_RUNNER}" "${SAMPLE_INFO}" "${OUTDIR}" "${GLOBAL_STEPS}"
    } > "${GLOBAL_JOB}"
    chmod +x "${GLOBAL_JOB}"

    if [[ "${NO_SUBMIT}" == true ]]; then
        if [[ -n "${SAMPLE_STEPS}" ]]; then
            printf 'sbatch --dependency=afterok:<sample_job_ids> %q\n' "${GLOBAL_JOB}"
        else
            printf 'sbatch %q\n' "${GLOBAL_JOB}"
        fi
    elif [[ "${#JOB_IDS[@]}" -gt 0 ]]; then
        DEPENDENCY="$(IFS=:; echo "${JOB_IDS[*]}")"
        SUMMARY_JOB_ID="$(sbatch --parsable --dependency="afterok:${DEPENDENCY}" "${GLOBAL_JOB}")"
        echo "已提交汇总作业：job ${SUMMARY_JOB_ID}，等待全部样本作业成功"
    else
        SUMMARY_JOB_ID="$(sbatch --parsable "${GLOBAL_JOB}")"
        echo "已提交汇总作业：job ${SUMMARY_JOB_ID}"
    fi
fi
