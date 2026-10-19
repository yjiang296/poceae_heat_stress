#!/usr/bin/env python3

### 脚本名称：05_extract_sample_counts.py
### 功能：按照 prepDE 公式从单个 StringTie GTF 计算 gene 和 transcript count
### 输入文件：StringTie 使用 -G -e 生成的 GTF
### 输出文件：单样本 gene_counts.tsv 和 transcript_counts.tsv
### 预期运行路径：任意目录或由 05_run_prepde.py 调用
### Conda 环境激活命令：conda activate rnaseq
### 使用示例：python scripts/05_extract_sample_counts.py --input sample.stringtie.gtf --outdir sample_counts --sample sample --read-length 150

import argparse
import csv
import math
import os
import random
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(description="从单样本 StringTie GTF 推算 count")
    parser.add_argument("--input", required=True, help="StringTie GTF")
    parser.add_argument("--outdir", required=True, help="单样本 count 输出目录")
    parser.add_argument("--sample", required=True, help="样本名称")
    parser.add_argument("--read-length", required=True, type=int, help="平均 read 长度")
    return parser.parse_args()


def parse_gtf(path, read_length):
    transcripts = {}
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] not in {"transcript", "exon"}:
                    continue
                attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
                transcript_id = attrs.get("transcript_id")
                if not transcript_id:
                    continue
                record = transcripts.setdefault(
                    transcript_id,
                    {"gene_id": attrs.get("gene_id", transcript_id), "coverage": 0.0, "length": 0},
                )
                if fields[2] == "transcript":
                    record["gene_id"] = attrs.get("gene_id", transcript_id)
                    record["coverage"] = max(0.0, float(attrs.get("cov", 0)))
                else:
                    record["length"] += abs(int(fields[4]) - int(fields[3])) + 1
    except (OSError, ValueError) as exc:
        print(f"错误：解析 StringTie GTF 失败：{path}：{exc}", file=sys.stderr)
        raise

    transcript_counts = {}
    gene_counts = defaultdict(int)
    for transcript_id, record in transcripts.items():
        count = int(math.ceil(record["coverage"] * record["length"] / read_length))
        transcript_counts[transcript_id] = count
        gene_counts[record["gene_id"]] += count
    return dict(gene_counts), transcript_counts


def write_counts(path, id_column, counts):
    try:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([id_column, "count"])
            writer.writerows((feature, counts[feature]) for feature in sorted(counts))
    except OSError as exc:
        print(f"错误：写出单样本 count 失败：{path}：{exc}", file=sys.stderr)
        raise


def main():
    args = parse_args()
    if args.read_length < 1:
        raise ValueError("--read-length 必须大于 0")
    print(f"运行目录：{os.getcwd()}")
    print(f"Python：{sys.executable}")
    print(f"Conda 环境：{os.environ.get('CONDA_DEFAULT_ENV')}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_counts, transcript_counts = parse_gtf(args.input, args.read_length)
    write_counts(outdir / f"{args.sample}.gene_counts.tsv", "gene_id", gene_counts)
    write_counts(
        outdir / f"{args.sample}.transcript_counts.tsv",
        "transcript_id", transcript_counts,
    )
    with open(outdir / "requirements.txt", "w", encoding="utf-8") as handle:
        subprocess.run(["pip", "freeze"], stdout=handle, check=True)
    print(
        f"验证：{args.sample} 输出 {len(gene_counts)} 个基因和 "
        f"{len(transcript_counts)} 个转录本 count"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"错误：单样本 count 提取失败：{exc}", file=sys.stderr)
        sys.exit(1)
