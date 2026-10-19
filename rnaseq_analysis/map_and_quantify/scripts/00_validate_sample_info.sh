#!/bin/bash
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -J validate_samples
#SBATCH -o ./logs/validate_samples.log
#SBATCH -e ./logs/validate_samples.err

### 脚本名称：00_validate_sample_info.sh
### 功能：检查样本信息表的表头、列数、名称唯一性和参考基因组字段一致性
### 输入文件：六列样本信息表
### 输出文件：验证日志
### 预期运行路径：项目根目录或由 main.sh 调用
### 环境名称：rnaseq
### 环境激活命令：conda activate rnaseq
### 使用示例：bash scripts/00_validate_sample_info.sh sample_info.txt results

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "用法：bash $0 <sample_info.txt> <outdir>" >&2
    exit 1
fi

INPUT="$1"
OUTDIR="$2"

if [[ ! -f "${INPUT}" ]]; then
    echo "错误：样本信息表不存在：${INPUT}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
RUN_LOG="${OUTDIR}/sample_info_validation.log"
{
    echo "运行目录：$(pwd)"
    echo "Conda 环境：${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
    echo "输入文件：${INPUT}"
} | tee "${RUN_LOG}"

awk '
BEGIN { FS="[[:space:]]+" }
NR == 1 {
    expected = "sample raw_fq_gz_R1 raw_fq_gz_R2 hisat2_index_prefix annotation_gtf reference_name"
    actual = $1" "$2" "$3" "$4" "$5" "$6
    if (NF != 6 || actual != expected) {
        print "错误：表头或列数不正确。要求：" expected > "/dev/stderr"
        exit 1
    }
    next
}
/^[[:space:]]*$/ || /^[[:space:]]*#/ { next }
NF != 6 {
    print "错误：第 " NR " 行不是 6 列" > "/dev/stderr"
    exit 1
}
$1 !~ /^[A-Za-z0-9._-]+$/ || $6 !~ /^[A-Za-z0-9._-]+$/ {
    print "错误：sample 和 reference_name 含不允许字符，第 " NR " 行" > "/dev/stderr"
    exit 1
}
seen[$1]++ {
    print "错误：sample 重复：" $1 > "/dev/stderr"
    exit 1
}
ref_index[$6] && (ref_index[$6] != $4 || ref_gtf[$6] != $5) {
    print "错误：同一 reference_name 对应不同索引或 GTF：" $6 > "/dev/stderr"
    exit 1
}
{
    ref_index[$6] = $4
    ref_gtf[$6] = $5
    count++
}
END {
    if (count == 0) {
        print "错误：样本表中没有数据行" > "/dev/stderr"
        exit 1
    }
    print "验证通过：共 " count " 个样本"
}
' "${INPUT}" | tee -a "${RUN_LOG}"

awk --version 2>&1 | head -n 1 > "${OUTDIR}/version.log" || true
