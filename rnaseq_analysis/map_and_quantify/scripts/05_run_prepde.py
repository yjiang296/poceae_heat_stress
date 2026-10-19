#!/usr/bin/env python3

### 脚本名称：05_run_prepde.py
### 功能：按 prepDE 公式从 StringTie GTF 计算基因和转录本 count 矩阵
### 输入文件：六列样本信息表、fastp JSON 和 StringTie GTF
### 输出文件：按参考基因组分组的 gene/transcript count TSV 及 read 长度记录
### 预期运行路径：项目根目录或由全局 Slurm 作业调用
### Conda 环境激活命令：conda activate rnaseq
### 使用示例：python scripts/05_run_prepde.py --input sample_info.txt --outdir results/summary/05_counts --results-root results

import argparse
import csv
import json
import os
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)


def parse_args():
    parser = argparse.ArgumentParser(description="按参考基因组生成 StringTie count 矩阵")
    parser.add_argument("--input", required=True, help="六列样本信息表")
    parser.add_argument("--outdir", required=True, help="count 矩阵输出目录")
    parser.add_argument("--results-root", required=True, help="流程结果根目录")
    return parser.parse_args()


def read_samples(path):
    expected = ["sample", "raw_fq_gz_R1", "raw_fq_gz_R2",
                "hisat2_index_prefix", "annotation_gtf", "reference_name"]
    samples = []
    try:
        with open(path, encoding="utf-8") as handle:
            if handle.readline().strip().split() != expected:
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


def get_clean_bases_reads(path):
    try:
        with open(path, encoding="utf-8") as handle:
            after = json.load(handle)["summary"]["after_filtering"]
        return int(after["total_bases"]), int(after["total_reads"])
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(f"错误：无法从 fastp JSON 获取 read 长度：{path}：{exc}", file=sys.stderr)
        raise


def read_count_table(path):
    counts = {}
    try:
        with open(path, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "count" not in reader.fieldnames:
                raise ValueError("count 表缺少表头")
            id_column = reader.fieldnames[0]
            for row in reader:
                counts[row[id_column]] = int(row["count"])
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"错误：读取单样本 count 失败：{path}：{exc}", file=sys.stderr)
        raise
    return counts


def write_matrix(path, feature_name, sample_names, sample_counts):
    features = sorted({feature for counts in sample_counts.values() for feature in counts})
    try:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([feature_name, *sample_names])
            for feature in features:
                writer.writerow([feature, *[sample_counts[sample].get(feature, 0) for sample in sample_names]])
    except OSError as exc:
        print(f"错误：写出 count 矩阵失败：{path}：{exc}", file=sys.stderr)
        raise


def main():
    args = parse_args()
    for label, value in [("运行目录", os.getcwd()), ("Python", sys.executable),
                         ("Conda 环境", os.environ.get("CONDA_DEFAULT_ENV"))]:
        print(f"{label}：{value}")
    groups = defaultdict(list)
    for sample in read_samples(args.input):
        groups[sample["reference_name"]].append(sample)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    metadata = []
    total_samples = 0
    extractor = Path(__file__).with_name("05_extract_sample_counts.py")
    if not extractor.is_file():
        raise FileNotFoundError(f"缺少单样本 count 提取脚本：{extractor}")
    for reference, samples in sorted(groups.items()):
        reference_dir = outdir / reference
        reference_dir.mkdir(parents=True, exist_ok=True)
        total_bases = 0
        total_reads = 0
        gtf_paths = {}
        for sample in samples:
            sample_name = sample["sample"]
            sample_dir = Path(args.results_root) / "samples" / sample_name
            gtf = sample_dir / "03_stringtie" / f"{sample_name}.stringtie.gtf"
            fastp_json = sample_dir / "01_fastp" / f"{sample_name}.fastp.json"
            if not gtf.is_file() or not fastp_json.is_file():
                raise FileNotFoundError(f"{sample_name} 缺少 StringTie GTF 或 fastp JSON")
            bases, reads = get_clean_bases_reads(fastp_json)
            total_bases += bases
            total_reads += reads
            gtf_paths[sample_name] = gtf

        if total_reads == 0:
            raise ValueError(f"{reference} 的 clean read 数为 0")
        average_read_length = max(1, round(total_bases / total_reads))
        gene_counts = {}
        transcript_counts = {}
        for sample_name, gtf in gtf_paths.items():
            sample_outdir = reference_dir / "sample_counts" / sample_name
            command = [
                sys.executable, str(extractor),
                "--input", str(gtf),
                "--outdir", str(sample_outdir),
                "--sample", sample_name,
                "--read-length", str(average_read_length),
            ]
            subprocess.run(command, check=True)
            gene_counts[sample_name] = read_count_table(sample_outdir / f"{sample_name}.gene_counts.tsv")
            transcript_counts[sample_name] = read_count_table(sample_outdir / f"{sample_name}.transcript_counts.tsv")
        sample_names = [sample["sample"] for sample in samples]
        write_matrix(reference_dir / f"{reference}.gene_count_matrix.tsv",
                     "gene_id", sample_names, gene_counts)
        write_matrix(reference_dir / f"{reference}.transcript_count_matrix.tsv",
                     "transcript_id", sample_names, transcript_counts)
        metadata.append({"reference_name": reference, "sample_count": len(samples),
                         "average_clean_read_length": average_read_length,
                         "count_type": "prepDE formula: ceil(coverage * exon_length / read_length)"})
        total_samples += len(samples)

    try:
        with open(outdir / "05_prepde_metadata.tsv", "w", encoding="utf-8", newline="") as handle:
            columns = ["reference_name", "sample_count", "average_clean_read_length", "count_type"]
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(metadata)
        with open(outdir / "requirements.txt", "w", encoding="utf-8") as handle:
            subprocess.run(["pip", "freeze"], stdout=handle, check=True)
    except OSError as exc:
        print(f"错误：写出 prepDE 元数据失败：{exc}", file=sys.stderr)
        raise
    print(f"验证：{len(groups)} 个参考基因组、{total_samples} 个样本均已生成 count 矩阵")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"错误：prepDE 汇总失败：{exc}", file=sys.stderr)
        sys.exit(1)
