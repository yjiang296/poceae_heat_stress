#!/bin/bash
#SBATCH -J submit_function
#SBATCH -p x86_64_4cpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o submit_function.%j.out
#SBATCH -e submit_function.%j.err
### 脚本名称：submit_annotation.sh；功能：读取样本表，为每个材料独立生成及提交作业。
### 输入：$1 两列样本表；输出：$2 根目录，每样本含作业、日志、01 至 04 结果。
### 预期路径：任意；环境：conda activate base；示例：bash submit_annotation.sh sample_info.txt results --no-submit
### 固定软件路径在 config.json；默认只生成作业，只有 --submit 才提交；无需内存参数。
set -euo pipefail
[[ $# -ge 2 ]] || { echo '用法：submit_annotation.sh 样本表 输出根目录 [--submit|--no-submit] [--steps all|02] [--from 03] [--samples A,B] [--resume] [--config 文件] [--threads 8] [--partition 名称]' >&2; exit 2; }
INPUT=$1; OUTDIR=$2; shift 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONFIG="$ROOT/config.json"; SUBMIT=0; RESUME=0; STEPS=all; SELECT=''; THREADS=''; PARTITION=''
while [[ $# -gt 0 ]]; do
 case "$1" in
 --submit) SUBMIT=1; shift;; --no-submit) SUBMIT=0; shift;; --resume) RESUME=1; shift;;
 --config) CONFIG=$2; shift 2;; --threads) THREADS=$2; shift 2;; --partition) PARTITION=$2; shift 2;;
 --steps) STEPS=$2; shift 2;; --samples) SELECT=$2; shift 2;;
 --from) case "$2" in 01) STEPS=all;; 02) STEPS=02,03,04;; 03) STEPS=03,04;; 04) STEPS=04;; *) exit 2;; esac; shift 2;;
 *) echo "未知参数：$1" >&2; exit 2;;
 esac
done
[[ -f "$INPUT" && -f "$CONFIG" ]] || { echo '样本表或配置不存在' >&2; exit 2; }
CONFIG=$(python3 -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' "$CONFIG")
[[ -n "$THREADS" ]] || THREADS=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["threads"])' "$CONFIG")
[[ -n "$PARTITION" ]] || PARTITION=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["partition"])' "$CONFIG")
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$PARTITION" =~ ^[A-Za-z0-9_]+$ ]] || exit 2
case "$STEPS" in all|01|02|03|04|01,02,03,04|02,03,04|03,04) ;; *) echo '步骤不合法' >&2; exit 2;; esac
mkdir -p "$OUTDIR"; OUTDIR=$(cd "$OUTDIR" && pwd)
[[ "$OUTDIR" != *[[:space:]]* ]] || { echo 'Slurm 输出根目录请勿含空白字符' >&2; exit 2; }
printf '%s 工作目录=%s 环境=%s\n' "$(date '+%F %T')" "$(pwd)" "${CONDA_DEFAULT_ENV:-未检测到 conda 环境}" | tee -a "$OUTDIR/run.log"
### 先验证整张表，避免读到后部错误时已提交前面的样本；输入路径必须为绝对路径且不含空白。
SAMPLES=(); FASTAS=(); FOUND=','
while read -r SAMPLE FASTA EXTRA || [[ -n "${SAMPLE:-}" ]]; do
 [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
 [[ "$SAMPLE" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && -z "${EXTRA:-}" && "$FASTA" == /* && -f "$FASTA" ]] || { echo "非法样本行：$SAMPLE $FASTA" >&2; exit 2; }
 [[ "$FOUND" != *,"$SAMPLE",* ]] || { echo "重复样本名：$SAMPLE" >&2; exit 2; }
 FOUND+="$SAMPLE,"
 if [[ -n "$SELECT" && ",$SELECT," != *,"$SAMPLE",* ]]; then continue; fi
 SAMPLES+=("$SAMPLE"); FASTAS+=("$FASTA")
done < "$INPUT"
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo '没有选择到样本' >&2; exit 2; }
if [[ -n "$SELECT" ]]; then
 IFS=',' read -r -a REQUESTED <<< "$SELECT"
 for SAMPLE in "${REQUESTED[@]}"; do [[ "$FOUND" == *,"$SAMPLE",* ]] || { echo "未找到样本：$SAMPLE" >&2; exit 2; }; done
fi
if [[ "$SUBMIT" == 1 ]]; then
 command -v sbatch >/dev/null; command -v squeue >/dev/null
 for SAMPLE in "${SAMPLES[@]}"; do
  LAST="$OUTDIR/$SAMPLE/jobs/last_submitted.txt"
  if [[ -f "$LAST" ]]; then
   [[ "$RESUME" == 1 ]] || { echo "已有提交记录：$SAMPLE；确认需要续跑时加 --resume" >&2; exit 1; }
   OLD=$(cat "$LAST"); [[ "$OLD" =~ ^[0-9]+$ ]] || exit 2
   ACTIVE=$(squeue -h -u "$(id -un)" -o '%A' | awk -v target="$OLD" '$1==target {print}')
   [[ -z "$ACTIVE" ]] || { echo "$SAMPLE 的作业 $OLD 仍在队列中" >&2; exit 1; }
  fi
 done
fi
STAMP=$(date '+%Y%m%d_%H%M%S')_$$
for I in "${!SAMPLES[@]}"; do
 SAMPLE=${SAMPLES[$I]}; FASTA=${FASTAS[$I]}; DIR="$OUTDIR/$SAMPLE"; mkdir -p "$DIR/jobs" "$DIR/logs"
 JOB="$DIR/jobs/annotation_${SAMPLE}_${STAMP}.slurm"; SNAP="$DIR/jobs/config_${STAMP}.json"; cp "$CONFIG" "$SNAP"
 {
  printf '#!/bin/bash\n#SBATCH -J anno_%s\n#SBATCH -p %s\n#SBATCH -N 1\n#SBATCH -n %s\n' "$SAMPLE" "$PARTITION" "$THREADS"
  printf '#SBATCH -o %s/logs/annotation_%s_%%j.out\n#SBATCH -e %s/logs/annotation_%s_%%j.err\n' "$DIR" "$STAMP" "$DIR" "$STAMP"
  printf '### 单材料作业：步骤 %s；通过独立配置快照复现。\nset -euo pipefail\n' "$STEPS"
  printf 'bash %q %q %q %q --sample %q --config %q --threads %q\n' "$ROOT/main.sh" "$FASTA" "$DIR" "$STEPS" "$SAMPLE" "$SNAP" "$THREADS"
 } > "$JOB"
 bash -n "$JOB"
 if [[ "$SUBMIT" == 1 ]]; then
  ID=$(sbatch --parsable "$JOB"); ID=${ID%%;*}; printf '%s\n' "$ID" > "$DIR/jobs/last_submitted.txt"
  printf '%s\t%s\t%s\t%s\n' "$(date '+%F %T')" "$SAMPLE" "$ID" "$JOB" | tee -a "$OUTDIR/submissions.tsv"
 else printf '仅生成，未提交：sbatch %q\n' "$JOB"; fi
done
printf 'bash %s\n' "$BASH_VERSION" >> "$OUTDIR/version.log"
