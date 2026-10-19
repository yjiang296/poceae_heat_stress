#!/bin/bash
### 脚本名称：update_github.sh；功能：检查并提交清单内源码，再推送当前仓库 main。
### 输入：$1 提交说明；输出：本地 Git 提交及 origin/main；运行路径：仓库任意位置。
### 环境：Git、Python 3、Bash；无需 conda activate；示例：bash update_github.sh "更新流程说明"
### 此脚本为本地版本同步入口，不提交 Slurm；不使用强制推送，不包含任何认证密钥。
set -euo pipefail
[[ $# == 1 && -n "$1" ]] || { echo '用法：bash update_github.sh "提交说明"' >&2; exit 2; }
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
[[ "$(git rev-parse --show-toplevel)" == "$ROOT" ]] || { echo '请在此仓库的完整 Git 克隆中运行' >&2; exit 1; }
[[ "$(git symbolic-ref --short HEAD)" == main ]] || { echo '请先切换到 main 分支' >&2; exit 1; }
case "$(git remote get-url --push origin)" in
 https://github.com/yjiang296/poceae_heat_stress.git|git@github.com:yjiang296/poceae_heat_stress.git) ;;
 *) echo 'origin 推送地址与目标仓库不符，请检查' >&2; exit 1;;
esac
git diff --cached --quiet || { echo '暂存区已有修改，请先自行处理，避免混入本次提交' >&2; exit 1; }
if git ls-remote --exit-code origin refs/heads/main >/dev/null; then
 git fetch origin main
 git merge-base --is-ancestor origin/main HEAD || { echo '远端已有新提交，请先检查并合并，再重试' >&2; exit 1; }
else
 STATUS=$?
 [[ "$STATUS" == 2 ]] || { echo '无法读取远端分支，请检查认证和网络' >&2; exit "$STATUS"; }
 echo '远端尚无 main 分支，将进行首次推送。'
fi
FILES=()
while IFS= read -r FILE || [[ -n "$FILE" ]]; do
 [[ -z "$FILE" || "$FILE" == \#* ]] && continue
 [[ "$FILE" != /* && "$FILE" != *..* && "$FILE" != *[[:space:]]* ]] || { echo "非法清单路径：$FILE" >&2; exit 1; }
 [[ ! -L "$FILE" ]] || { echo "不提交符号链接：$FILE" >&2; exit 1; }
 FILES+=("$FILE")
done < publish_files.txt
[[ ${#FILES[@]} -gt 0 ]] || { echo '发布清单为空' >&2; exit 1; }
bash -n update_github.sh
bash -n protein_function_annotation/main.sh
bash -n protein_function_annotation/submit_annotation.sh
(cd protein_function_annotation && python3 -m unittest discover -s tests -v)
bash -n LAI/submit_ltr_pipeline.sh
python3 LAI/test_lai_t6.py
git add -- "${FILES[@]}"
git diff --cached --check
if git diff --cached --quiet; then
 echo '清单内没有新修改；尝试推送已有本地提交。'
else
 git diff --cached --stat
 git commit -m "$1"
fi
git push origin HEAD:main
echo '已同步到 https://github.com/yjiang296/poceae_heat_stress'
