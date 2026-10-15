# 软件安装与环境配置

本仓库不包含软件安装包或数据库。实际生产验证使用以下固定版本；本页说明安装来源和环境调用方式，不自动重新安装。

## InterProScan

主程序和配套数据库来自 EBI 官方 [5.78-109.0 Linux 64-bit 包](https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz)。原部署下载官方 MD5 文件校验后解压，运行随包测试。InterProScan 本体并非 Conda 安装，Java/Python 由独立 Conda 环境提供：

```bash
conda create -y -p /absolute/path/to/software/miniforge3/envs/ipr_runtime_5.78 \
  --override-channels -c conda-forge --solver libmamba python=3.11 openjdk=11
```

将 `config.json` 中 `interproscan` 指向实际 `interproscan.sh`，`ipr_runtime` 指向上述环境根目录，`python` 指向工作流使用的 Python 可执行文件。

## eggNOG-mapper

原部署复用了独立的 Conda 环境，安装记录为：

```bash
conda create -y -p /absolute/path/to/software/miniforge3/envs/eggnog_mapper_2.1.13 \
  -c conda-forge -c bioconda eggnog-mapper=2.1.13
```

数据库独立于软件环境。原部署从 [eggNOG emapperdb 5.0.2](http://eggnog5.embl.de/download/emapperdb-5.0.2/) 下载 `eggnog.db.gz`、`eggnog_proteins.dmnd.gz`、`eggnog.taxa.tar.gz`，检查压缩完整性后解压。经工具检查数据库版本为 5.0.2，DIAMOND 可读索引。

`eggnog_env` 设置为环境根目录，`eggnog_data` 设置为包含 eggnog.db、eggnog_proteins.dmnd 等文件的目录。新建环境时仍需确认依赖版本；安装命令本身不保证以后解析出与历史环境完全相同的依赖。

## 不手动激活也能运行的原因

`02_interpro.py` 为工具子进程设置 `PATH=<ipr_runtime>/bin:原PATH`、`JAVA_HOME=<ipr_runtime>`；`03_eggnog.py` 设置 `PATH=<eggnog_env>/bin:原PATH`、`EGGNOG_DATA_DIR=<eggnog_data>`，以绝对路径调用 emapper.py，入口 shebang 使用其环境的 Python。

这些设置不改变交互终端，不等价于执行所有 Conda 激活钩子；该方式已针对这里的两套环境实际验证。手动使用工具时仍可 `conda activate /absolute/path/to/env`。

环境需要迁移时，Conda 前缀可能包含绝对路径，可使用 [conda create --clone --prefix](https://docs.conda.io/projects/conda/en/stable/commands/create.html) 重建前缀并重新验证；不要直接移动环境后假定可用。大型数据库可以独立集中存放。修改配置后使用新输出目录进行小规模测试，不改写旧完成签名。
