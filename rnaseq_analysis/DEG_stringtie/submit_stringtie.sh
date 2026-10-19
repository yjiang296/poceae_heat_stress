#!/bin/bash
### 脚本名称：submit_stringtie.sh
### 功能：读取材料信息表，冻结代码并为每材料生成独立Slurm作业；默认提交。
### 输入：$1为sample_info.txt；$2为新批次输出目录；输出：代码快照、Slurm脚本和提交清单。
### 预期运行路径：流程根目录；环境名称：poaceae_hsf_core。
### 使用示例：bash submit_stringtie.sh config/materials/sample_info.txt runs/批次名 --no-submit。
set -euo pipefail
[[ $# -ge 2 ]] || { echo "用法：bash submit_stringtie.sh sample_info.txt 新批次目录 [--no-submit]" >&2; exit 2; }
INPUT="$1"; OUTDIR="$2"; shift 2
EXTRA=("$@")
set --
[[ -f "$INPUT" ]] || { echo "材料表不存在：$INPUT" >&2; exit 2; }
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$ROOT/config/runtime.sh"
CONDA_DEFAULT_ENV=poaceae_hsf_core "$PIPE_PYTHON" "$ROOT/scripts/00_submit.py" --input "$INPUT" --outdir "$OUTDIR" "${EXTRA[@]}"
