# poceae_heat_stress

禾本科比较基因组研究的可复用分析流程。

## 蛋白功能注释

[protein_function_annotation](protein_function_annotation/README.md)：对已有基因模型对应的蛋白 FASTA 运行 InterProScan 5 和 eggNOG-mapper，按材料生成独立 Slurm 作业，支持输入检查、分块续跑和原蛋白 ID 恢复。

仓库保存代码、测试、配置模板和说明。实际软件环境、数据库、输入序列、服务器配置和计算结果保存在运行机器上。

## 一条命令更新 GitHub

首次使用需要 Git、Python 3，以及本机可用的 GitHub Git 写入认证（SSH key 或 HTTPS 凭据管理器）。Codex 的 GitHub 连接授权不会自动变成本机 Git 凭据。不要把令牌写入脚本或仓库。

```bash
git clone https://github.com/yjiang296/poceae_heat_stress.git
cd poceae_heat_stress
# 修改本克隆中的代码或说明后：
bash update_github.sh "说明本次修改"
```

脚本在当前克隆的 `main` 分支运行，检查 origin 地址、远端进度、Shell 语法和回归测试，然后只暂存 `publish_files.txt` 列出的文件、提交并推送。新增源文件需要先将其加入清单。不会自动复制超算或其他目录中的改动，不会自动合并远端更新，也不会强制推送。

如果本机尚未配置 Git 写入认证，可继续让 Codex 通过已连接的 GitHub 工具上传。脚本中断于推送时，本地提交仍保留；解决认证后再次运行即可重试。Git 提交身份使用本机已有的 `user.name` 与 `user.email`。

同步运行中的工作流时请保留原代码和配置：已完成阶段使用代码/配置指纹校验，更新后应另建输出目录。
