#!/bin/bash
### 脚本名称: plot_fire_tss.sh
### 功能: 读取样本配置表，为每个样本独立生成并提交 FIRE TSS deepTools Slurm 作业
### 输入文件: sample_info_config.txt，依次为样本名、bigWig、TSS BED、输出目录
### 输出文件: 每个样本的矩阵、平均曲线、热图、source data、日志及 Slurm 脚本
### 预期运行路径: 任意具有 sbatch 命令且可访问输入文件的提交节点目录
### 环境名称: 无需激活 conda，deepTools 已加入 PATH
### 使用示例: bash plot_fire_tss.sh sample_info_config.txt
### 使用示例: bash plot_fire_tss.sh sample_info_config.txt --no-submit

set -euo pipefail

usage() {
    cat <<'EOF'
用法:
  bash plot_fire_tss.sh <sample_info_config.txt> [选项]

配置表（空白或 Tab 分隔）:
  #sample  bigwig_path  tss_bed_path  output_dir
  sample1  /path/a.bw   /path/a.tss.bed  /path/sample1_result

每行必须有四列，建议使用绝对路径；相对路径以运行命令时的目录为基准。
输出目录不存在时自动创建；全表检查通过后才提交作业。
不指定 --y-max / --z-max 时仍为独立自动缩放。

选项:
  --partition STR   Slurm 分区，默认: x86_64_4cpu
  --threads INT     每个样本使用的 CPU 数，默认: 2
  --before INT      TSS 上游距离，默认: 2000
  --after INT       TSS 下游距离，默认: 2000
  --bin-size INT    bin 大小，默认: 20
  --y-min NUM       所有平均曲线（含热图顶部）的下限，默认: 0
  --y-max NUM       所有平均曲线的统一上限；省略时各样本自动缩放
  --z-min NUM       所有热图的色标下限，默认: 0
  --z-max NUM       所有热图的统一色标上限；省略时各样本自动缩放
  --plot-only       使用已有矩阵重新绘图，不重算 computeMatrix
  --exclude-missing-contigs  排除 bigWig 中不存在的 contig 上的 TSS；默认报错
  --no-submit       仅生成 Slurm 脚本，不执行 sbatch
  -h, --help        显示帮助
EOF
}

die() { echo "错误: $*" >&2; exit 1; }
log_time() { echo "[$(date '+%F %T')] $*"; }

absolute_file() {
    local path=$1
    [[ -f "$path" ]] || die "文件不存在: $path"
    printf '%s/%s\n' "$(cd "$(dirname "$path")" && pwd -P)" "$(basename "$path")"
}

absolute_dir() {
    local path=$1
    mkdir -p -- "$path" || die "无法创建输出目录: $path"
    (cd "$path" && pwd -P)
}

run_worker() {
    local sample=$1 bigwig=$2 tss_bed=$3 outdir=$4 before=$5 after=$6 bin_size=$7
    local y_min=$8 y_max=$9 z_min=${10} z_max=${11} plot_only=${12}
    local exclude_missing=${13}
    local -a profile_scale=(--yMin "$y_min") heatmap_scale=(--zMin "$z_min")
    [[ $y_max == auto ]] || profile_scale+=(--yMax "$y_max")
    [[ $z_max == auto ]] || heatmap_scale+=(--zMax "$z_max")
    local threads=${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-2}}
    # 数值库保持单线程；computeMatrix 的并行数仍由 --numberOfProcessors 控制。
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    local matrix="$outdir/${sample}.TSS.matrix.gz"
    local sorted_bed="$outdir/${sample}.TSS.sorted.bed"
    local run_log="$outdir/${sample}.TSS.run.log"

    exec > >(tee -a "$run_log") 2>&1
    log_time "样本 $sample 开始运行"
    echo "pwd=$(pwd)"
    echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-未检测到 conda 环境}"
    echo "SLURM_JOB_ID=${SLURM_JOB_ID:-非 Slurm 环境}"
    echo "profile_y=[$y_min,$y_max]; heatmap_z=[$z_min,$z_max]; plot_only=$plot_only"
    echo "numeric_library_threads=1; computeMatrix_processors=$threads"

    for tool in computeMatrix plotProfile plotHeatmap; do
        command -v "$tool" >/dev/null 2>&1 || die "计算节点找不到 $tool"
    done

    if (( plot_only == 0 )); then
    # TSS 检查和处理内嵌于本脚本，不需要额外的 Python 文件。
    python3 - "$bigwig" "$tss_bed" "$outdir/$sample.TSS" "$exclude_missing" <<'PY'
import sys
from pathlib import Path
import pyBigWig
bw, bed, prefix, exclude = sys.argv[1:]
with pyBigWig.open(bw) as handle:
    chroms = handle.chroms()
total = kept = missing = 0
with open(bed) as src, open(prefix+'.used.bed','w') as good, open(prefix+'.excluded.bed','w') as dropped:
    for n, line in enumerate(src,1):
        if not line.strip() or line.startswith(('#','track ','browser ')):
            continue
        x=line.split()
        if len(x)<6 or not x[1].isdigit() or not x[2].isdigit() or x[5] not in ('+','-'):
            raise SystemExit(f'Invalid BED6: {bed}, line {n}')
        start,end=int(x[1]),int(x[2])
        if not 0 <= start < end:
            raise SystemExit(f'Invalid coordinates: {bed}, line {n}')
        total += 1
        record='\t'.join(x[:6])+'\n'
        if x[0] not in chroms:
            missing += 1
            dropped.write(record)
        elif end > chroms[x[0]]:
            raise SystemExit(f'Coordinates exceed bigWig chromosome: {bed}, line {n}')
        else:
            kept += 1
            good.write(record)
Path(prefix+'.input_qc.tsv').write_text('total_tss\tretained_tss\tmissing_contig_tss\n'+f'{total}\t{kept}\t{missing}\n')
print(f'TSS: total={total}, retained={kept}, missing_contig={missing}', flush=True)
if not kept:
    raise SystemExit('No TSS can be located in this bigWig')
if missing and exclude != '1':
    raise SystemExit('Missing contigs: inspect .excluded.bed; use --exclude-missing-contigs only if exclusion is intended')
PY
    tss_bed="$outdir/$sample.TSS.used.bed"
    log_time "Step 1: 计算 TSS 信号矩阵"
    computeMatrix reference-point \
        --referencePoint TSS -R "$tss_bed" -S "$bigwig" \
        -b "$before" -a "$after" --binSize "$bin_size" \
        --averageTypeBins mean --missingDataAsZero \
        --sortRegions descend --sortUsing mean \
        --samplesLabel "$sample" --numberOfProcessors "$threads" \
        --outFileName "$matrix" \
        --outFileNameMatrix "$outdir/${sample}.TSS.matrix.tsv" \
        --outFileSortedRegions "$sorted_bed"
    fi
    [[ -s "$matrix" && -s "$sorted_bed" ]] || die "Step 1 输出缺失"

    log_time "Step 2: 绘制平均曲线"
    plotProfile -m "$matrix" -out "$outdir/${sample}.TSS.profile.pdf" \
        "${profile_scale[@]}" \
        --plotType lines --averageType mean --colors '#D55E00' \
        --plotTitle "FIRE signal around TSS: $sample" --refPointLabel TSS \
        --yAxisLabel "Mean FIRE coverage" --regionsLabel Genes \
        --outFileNameData "$outdir/${sample}.TSS.profile.tsv"
    plotProfile -m "$matrix" -out "$outdir/${sample}.TSS.profile.png" \
        "${profile_scale[@]}" \
        --plotType lines --averageType mean --colors '#D55E00' \
        --plotTitle "FIRE signal around TSS: $sample" --refPointLabel TSS \
        --yAxisLabel "Mean FIRE coverage" --regionsLabel Genes --dpi 300

    log_time "Step 3: 绘制热图"
    plotHeatmap -m "$matrix" -out "$outdir/${sample}.TSS.heatmap.pdf" \
        "${profile_scale[@]}" "${heatmap_scale[@]}" --averageTypeSummaryPlot mean \
        --colorMap YlOrRd --plotTitle "FIRE signal around TSS: $sample" \
        --refPointLabel TSS --yAxisLabel Genes --regionsLabel Genes \
        --whatToShow "plot, heatmap and colorbar"
    plotHeatmap -m "$matrix" -out "$outdir/${sample}.TSS.heatmap.png" \
        "${profile_scale[@]}" "${heatmap_scale[@]}" --averageTypeSummaryPlot mean \
        --colorMap YlOrRd --plotTitle "FIRE signal around TSS: $sample" \
        --refPointLabel TSS --yAxisLabel Genes --regionsLabel Genes \
        --whatToShow "plot, heatmap and colorbar" --dpi 300

    [[ -s "$outdir/${sample}.TSS.profile.pdf" ]] || die "平均曲线 PDF 缺失"
    [[ -s "$outdir/${sample}.TSS.heatmap.pdf" ]] || die "热图 PDF 缺失"
    {
        computeMatrix --version
        plotProfile --version
        plotHeatmap --version
    } > "$outdir/${sample}.TSS.version.log" 2>&1
    log_time "样本 $sample 运行完成"
}

if [[ ${1:-} == "--worker" ]]; then
    [[ $# -eq 14 ]] || die "worker 参数数量错误"
    shift
    run_worker "$@"
    exit 0
fi

[[ $# -ge 1 ]] || { usage; exit 1; }
[[ $1 != "-h" && $1 != "--help" ]] || { usage; exit 0; }

config=$1
shift
partition=x86_64_4cpu
threads=2
before=2000
after=2000
bin_size=20
no_submit=0
y_min=0
y_max=auto
z_min=0
z_max=auto
plot_only=0
exclude_missing=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --partition) [[ $# -ge 2 ]] || die "--partition 缺少参数"; partition=$2; shift 2 ;;
        --threads) [[ $# -ge 2 ]] || die "--threads 缺少参数"; threads=$2; shift 2 ;;
        --before) [[ $# -ge 2 ]] || die "--before 缺少参数"; before=$2; shift 2 ;;
        --after) [[ $# -ge 2 ]] || die "--after 缺少参数"; after=$2; shift 2 ;;
        --bin-size) [[ $# -ge 2 ]] || die "--bin-size 缺少参数"; bin_size=$2; shift 2 ;;
        --y-min) [[ $# -ge 2 ]] || die "--y-min 缺少参数"; y_min=$2; shift 2 ;;
        --y-max) [[ $# -ge 2 ]] || die "--y-max 缺少参数"; y_max=$2; shift 2 ;;
        --z-min) [[ $# -ge 2 ]] || die "--z-min 缺少参数"; z_min=$2; shift 2 ;;
        --z-max) [[ $# -ge 2 ]] || die "--z-max 缺少参数"; z_max=$2; shift 2 ;;
        --plot-only) plot_only=1; shift ;;
        --exclude-missing-contigs) exclude_missing=1; shift ;;
        --no-submit) no_submit=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "无法识别的参数: $1" ;;
    esac
done
for value in "$y_min" "$z_min" "$y_max" "$z_max"; do
    [[ $value == auto || $value =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]] || die "坐标范围必须为数值: $value"
done
[[ $y_min != auto && $z_min != auto ]] || die "坐标下限必须为数值"
[[ $y_max == auto ]] || awk -v lo="$y_min" -v hi="$y_max" 'BEGIN {exit !(hi>lo)}' || die "y-max 必须大于 y-min"
[[ $z_max == auto ]] || awk -v lo="$z_min" -v hi="$z_max" 'BEGIN {exit !(hi>lo)}' || die "z-max 必须大于 z-min"
if [[ $y_max == auto || $z_max == auto ]]; then
    echo '提示: 未同时设置 --y-max 和 --z-max；跨样本的曲线或色标仍会自动缩放。' >&2
fi

[[ -f "$config" ]] || die "配置文件不存在: $config"
for value in "$threads" "$before" "$after" "$bin_size"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "线程数和窗口参数必须为正整数"
done
(( before % bin_size == 0 && after % bin_size == 0 )) || die "上下游窗口必须是 bin-size 的整数倍"
[[ $partition =~ ^[A-Za-z0-9._-]+$ ]] || die "分区名格式不正确"
(( no_submit == 1 )) || command -v sbatch >/dev/null 2>&1 || die "找不到 sbatch；检查环境或使用 --no-submit"

script_path=$(absolute_file "$0")
sample_count=0
seen_samples="|"
slurm_files=()

while read -r sample bigwig tss_bed outdir extra || [[ -n ${sample:-} ]]; do
    [[ -n ${sample:-} ]] || continue
    [[ $sample == \#* ]] && continue
    [[ -z ${extra:-} ]] || die "配置行字段超过 4 列: $sample"
    [[ $sample =~ ^[A-Za-z0-9._-]+$ ]] || die "样本名只能包含字母、数字、点、下划线和短横线: $sample"
    [[ $sample != . && $sample != .. ]] || die "非法样本名: $sample"
    [[ -n ${bigwig:-} && -n ${tss_bed:-} && -n ${outdir:-} ]] || die "配置行必须有 4 列: $sample"
    [[ $seen_samples != *"|${sample}|"* ]] || die "样本名重复: $sample"
    seen_samples+="${sample}|"

    bigwig=$(absolute_file "$bigwig")
    tss_bed=$(absolute_file "$tss_bed")
    outdir=$(absolute_dir "$outdir")
    if (( plot_only == 1 )); then
        [[ -s "$outdir/${sample}.TSS.matrix.gz" && -s "$outdir/${sample}.TSS.sorted.bed" ]] || die "重绘所需矩阵或 sorted BED 不存在: $sample"
    fi
    slurm_file="$outdir/${sample}.plot_fire_tss.slurm"
    stdout_log="$outdir/${sample}.plot_fire_tss.out"
    stderr_log="$outdir/${sample}.plot_fire_tss.err"

    (
        export LC_ALL=C
        echo '#!/bin/bash'
        echo "#SBATCH -J fire_tss_${sample}"
        echo "#SBATCH -p ${partition}"
        echo '#SBATCH -N 1'
        echo '#SBATCH --ntasks=1'
        echo "#SBATCH --cpus-per-task=${threads}"
        echo "#SBATCH -o ${stdout_log}"
        echo "#SBATCH -e ${stderr_log}"
        echo '### 功能: 在单个样本中绘制 FIRE TSS 平均曲线和热图'
        echo '### 环境: deepTools 已加入 PATH，无需激活 conda'
        echo 'set -euo pipefail'
        printf 'bash %q --worker' "$script_path"
        printf ' %q' "$sample" "$bigwig" "$tss_bed" "$outdir" "$before" "$after" "$bin_size" \
            "$y_min" "$y_max" "$z_min" "$z_max" "$plot_only" "$exclude_missing"
        printf '\n'
    ) > "$slurm_file"
    chmod +x "$slurm_file"
    slurm_files+=("$slurm_file")
    sample_count=$((sample_count + 1))
done < "$config"

(( sample_count > 0 )) || die "配置文件中没有有效样本"
# 全部样本通过路径检查后再提交，避免读到坏行前已提交部分作业。
for slurm_file in "${slurm_files[@]}"; do
    if (( no_submit == 1 )); then
        echo "[仅生成] sbatch $slurm_file"
    else
        log_time "提交 $slurm_file"
        sbatch "$slurm_file"
    fi
done
log_time "处理完成，共 ${sample_count} 个样本"
