### 脚本名称：00_common.py
### 功能：统一参数、中文日志、表格读写和软件版本记录，供各步骤调用。
### 输入文件：调用脚本的 --input 和 --outdir；输出文件：run.log、requirements.txt。
### 预期运行路径：流程根目录；环境激活：conda activate poaceae_hsf_core。
### 使用示例：由 python scripts/02_matrices.py --input 上游目录 --outdir 结果目录导入。
import argparse
import hashlib
import os
from pathlib import Path
import random
import subprocess
import sys
import numpy as np
import pandas as pd

random.seed(42)
np.random.seed(42)

def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', required=True)
    return p

def require(ok, message):
    if not ok:
        raise ValueError(message)

def read(path):
    try:
        return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, comment='#')
    except Exception as e:
        raise OSError(f'读取失败：{path}：{e}') from e

def write(frame, path):
    try:
        frame.to_csv(path, sep='\t', index=False, na_rep='NA')
    except Exception as e:
        raise OSError(f'写入失败：{path}：{e}') from e

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for part in iter(lambda: f.read(1048576), b''):
            h.update(part)
    return h.hexdigest()

def start(args):
    out = Path(args.outdir)
    require('data/raw' not in out.as_posix(), '禁止写入 data/raw')
    out.mkdir(parents=True, exist_ok=True)
    def log(message):
        print(message, flush=True)
        with (out/'run.log').open('a') as f:
            f.write(str(message)+'\n')
    log(f'运行目录：{os.getcwd()}；Python：{sys.executable}；Conda：{os.environ.get("CONDA_DEFAULT_ENV")}')
    return out, log

def run(main):
    try:
        args = main()
        with open(Path(args.outdir)/'requirements.txt', 'w') as f:
            subprocess.run([sys.executable, '-m', 'pip', 'freeze'], stdout=f, check=True)
    except Exception as e:
        print(f'错误：{e}', file=sys.stderr)
        sys.exit(1)
