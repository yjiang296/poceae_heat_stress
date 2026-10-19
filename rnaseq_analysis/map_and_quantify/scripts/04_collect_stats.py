#!/usr/bin/env python3

### 脚本名称：04_collect_stats.py
### 功能：汇总 fastp JSON 和 HISAT2 summary 中的关键质控与比对指标
### 输入文件：六列样本信息表及 results-root 下的样本结果
### 输出文件：fastp、HISAT2 和合并统计 TSV
### 预期运行路径：项目根目录或由全局 Slurm 作业调用
### Conda 环境激活命令：conda activate rnaseq
### 使用示例：python scripts/04_collect_stats.py --input sample_info.txt --outdir results/summary/04_stats --results-root results

import argparse
import csv
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(description="汇总 fastp 和 HISAT2 统计信息")
    parser.add_argument("--input", required=True, help="六列样本信息表")
    parser.add_argument("--outdir", required=True, help="统计结果输出目录")
    parser.add_argument("--results-root", required=True, help="流程结果根目录")
    return parser.parse_args()


def read_samples(path):
    samples = []
    try:
        with open(path, encoding="utf-8") as handle:
            header = handle.readline().strip().split()
            expected = [
                "sample", "raw_fq_gz_R1", "raw_fq_gz_R2",
                "hisat2_index_prefix", "annotation_gtf", "reference_name",
            ]
            if header != expected:
                raise ValueError("样本表表头不符合要求")
            for line_number, line in enumerate(handle, start=2):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.strip().split()
                if len(fields) != 6:
                    raise ValueError(f"第 {line_number} 行不是 6 列")
                samples.append(dict(zip(expected, fields)))
    except (OSError, ValueError) as exc:
        print(f"错误：读取样本表失败：{exc}", file=sys.stderr)
        raise
    return samples


def percent(value):
    return round(float(value) * 100.0, 4)


def parse_fastp(path):
    try:
        with open(path, encoding="utf-8") as handle:
            report = json.load(handle)
        before = report["summary"]["before_filtering"]
        after = report["summary"]["after_filtering"]
        filtering = report.get("filtering_result", {})
        duplication = report.get("duplication", {})
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(f"错误：解析 fastp JSON 失败：{path}：{exc}", file=sys.stderr)
        raise
    return {
        "raw_reads": before["total_reads"],
        "raw_bases": before["total_bases"],
        "raw_q20_percent": percent(before["q20_rate"]),
        "raw_q30_percent": percent(before["q30_rate"]),
        "raw_gc_percent": percent(before["gc_content"]),
        "clean_reads": after["total_reads"],
        "clean_bases": after["total_bases"],
        "clean_q20_percent": percent(after["q20_rate"]),
        "clean_q30_percent": percent(after["q30_rate"]),
        "clean_gc_percent": percent(after["gc_content"]),
        "passed_filter_reads": filtering.get("passed_filter_reads", ""),
        "low_quality_reads": filtering.get("low_quality_reads", ""),
        "too_many_N_reads": filtering.get("too_many_N_reads", ""),
        "too_short_reads": filtering.get("too_short_reads", ""),
        "duplication_percent": percent(duplication.get("rate", 0)),
    }


def matched_count(text, pattern):
    match = re.search(pattern, text, flags=re.MULTILINE)
    return int(match.group(1)) if match else 0


def parse_hisat2(path):
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        print(f"错误：读取 HISAT2 summary 失败：{path}：{exc}", file=sys.stderr)
        raise
    total_match = re.search(r"^\s*(\d+)\s+reads; of these:", text, flags=re.MULTILINE)
    rate_match = re.search(r"([\d.]+)%\s+overall alignment rate", text)
    if not total_match or not rate_match:
        raise ValueError(f"无法识别 HISAT2 summary 格式：{path}")
    total_pairs = int(total_match.group(1))
    unique_pairs = matched_count(
        text, r"^\s*(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time",
    )
    multi_pairs = matched_count(
        text, r"^\s*(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times",
    )
    return {
        "input_pairs": total_pairs,
        "concordant_unique_pairs": unique_pairs,
        "concordant_unique_percent": round(unique_pairs * 100 / total_pairs, 4),
        "concordant_multi_pairs": multi_pairs,
        "concordant_multi_percent": round(multi_pairs * 100 / total_pairs, 4),
        "overall_alignment_percent": float(rate_match.group(1)),
    }


def write_tsv(path, rows, columns):
    try:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    except OSError as exc:
        print(f"错误：写出 TSV 失败：{path}：{exc}", file=sys.stderr)
        raise


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"运行目录：{os.getcwd()}")
    print(f"Python：{sys.executable}")
    print(f"Conda 环境：{os.environ.get('CONDA_DEFAULT_ENV')}")

    fastp_rows = []
    hisat2_rows = []
    combined_rows = []
    for row in read_samples(args.input):
        sample = row["sample"]
        reference = row["reference_name"]
        sample_dir = Path(args.results_root) / "samples" / sample
        fastp_path = sample_dir / "01_fastp" / f"{sample}.fastp.json"
        hisat2_path = sample_dir / "02_hisat2" / f"{sample}.hisat2_summary.txt"
        if not fastp_path.is_file() or not hisat2_path.is_file():
            raise FileNotFoundError(f"{sample} 缺少 fastp JSON 或 HISAT2 summary")
        fastp_data = {"sample": sample, "reference_name": reference, **parse_fastp(fastp_path)}
        hisat2_data = {"sample": sample, "reference_name": reference, **parse_hisat2(hisat2_path)}
        fastp_rows.append(fastp_data)
        hisat2_rows.append(hisat2_data)
        combined_rows.append({**fastp_data, **{k: v for k, v in hisat2_data.items() if k not in fastp_data}})

    fastp_columns = list(fastp_rows[0])
    hisat2_columns = list(hisat2_rows[0])
    combined_columns = list(combined_rows[0])
    write_tsv(outdir / "04_fastp_stats.tsv", fastp_rows, fastp_columns)
    write_tsv(outdir / "04_hisat2_alignment_stats.tsv", hisat2_rows, hisat2_columns)
    write_tsv(outdir / "04_qc_alignment_summary.tsv", combined_rows, combined_columns)
    with open(outdir / "requirements.txt", "w", encoding="utf-8") as handle:
        subprocess.run(["pip", "freeze"], stdout=handle, check=True)
    print(f"验证：共汇总 {len(combined_rows)} 个样本，无遗漏")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"错误：统计汇总失败：{exc}", file=sys.stderr)
        sys.exit(1)
