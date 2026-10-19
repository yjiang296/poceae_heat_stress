#!/usr/bin/env python3

### 脚本名称：06_build_gene_matrices.py
### 功能：按参考基因组汇总 StringTie 基因 TPM/FPKM 矩阵
### 输入文件：六列样本信息表和 gene_abund.tsv
### 输出文件：按参考基因组分组的 gene TPM/FPKM 矩阵
### 预期运行路径：项目根目录或由全局 Slurm 作业调用
### Conda 环境激活命令：conda activate rnaseq
### 使用示例：python scripts/06_build_gene_matrices.py --input sample_info.txt --outdir results/summary/06_abundance --results-root results

import argparse
import csv
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
    parser = argparse.ArgumentParser(description="构建 StringTie gene TPM 和 FPKM 矩阵")
    parser.add_argument("--input", required=True, help="六列样本信息表")
    parser.add_argument("--outdir", required=True, help="丰度矩阵输出目录")
    parser.add_argument("--results-root", required=True, help="流程结果根目录")
    return parser.parse_args()


def read_samples(path):
    expected = [
        "sample", "raw_fq_gz_R1", "raw_fq_gz_R2",
        "hisat2_index_prefix", "annotation_gtf", "reference_name",
    ]
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


def parse_gene_abundance(path):
    records = {}
    try:
        with open(path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"Gene ID", "FPKM", "TPM"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(f"gene_abund.tsv 缺少字段：{required}")
            for row in reader:
                gene_id = row["Gene ID"]
                records[gene_id] = {
                    "gene_name": row.get("Gene Name", ""),
                    "TPM": float(row["TPM"]),
                    "FPKM": float(row["FPKM"]),
                }
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误：解析基因丰度文件失败：{path}：{exc}", file=sys.stderr)
        raise
    return records


def write_matrix(path, feature_id, annotations, sample_names, values, metric):
    feature_ids = sorted({item for sample_data in values.values() for item in sample_data})
    columns = [feature_id, *annotations, *sample_names]
    try:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for item in feature_ids:
                annotation_values = []
                for annotation in annotations:
                    value = ""
                    for sample in sample_names:
                        if item in values[sample]:
                            value = values[sample][item].get(annotation, "")
                            if value:
                                break
                    annotation_values.append(value)
                row = [item, *annotation_values]
                row.extend(values[sample].get(item, {}).get(metric, 0) for sample in sample_names)
                writer.writerow(row)
    except OSError as exc:
        print(f"错误：写出表达矩阵失败：{path}：{exc}", file=sys.stderr)
        raise
    return len(feature_ids)


def main():
    args = parse_args()
    print(f"运行目录：{os.getcwd()}")
    print(f"Python：{sys.executable}")
    print(f"Conda 环境：{os.environ.get('CONDA_DEFAULT_ENV')}")
    groups = defaultdict(list)
    for sample in read_samples(args.input):
        groups[sample["reference_name"]].append(sample)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    verification = []
    for reference, samples in sorted(groups.items()):
        reference_dir = outdir / reference
        reference_dir.mkdir(parents=True, exist_ok=True)
        gene_values = {}
        sample_names = []
        for row in samples:
            sample = row["sample"]
            sample_names.append(sample)
            stringtie_dir = Path(args.results_root) / "samples" / sample / "03_stringtie"
            gene_path = stringtie_dir / f"{sample}.gene_abund.tsv"
            if not gene_path.is_file():
                raise FileNotFoundError(f"{sample} 缺少 StringTie gene_abund.tsv")
            gene_values[sample] = parse_gene_abundance(gene_path)

        gene_count = write_matrix(
            reference_dir / f"{reference}.gene_tpm_matrix.tsv",
            "gene_id", ["gene_name"], sample_names, gene_values, "TPM",
        )
        write_matrix(
            reference_dir / f"{reference}.gene_fpkm_matrix.tsv",
            "gene_id", ["gene_name"], sample_names, gene_values, "FPKM",
        )
        verification.append([reference, len(sample_names), gene_count])

    try:
        with open(outdir / "06_matrix_verification.tsv", "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["reference_name", "sample_count", "gene_count"])
            writer.writerows(verification)
        with open(outdir / "requirements.txt", "w", encoding="utf-8") as handle:
            subprocess.run(["pip", "freeze"], stdout=handle, check=True)
    except OSError as exc:
        print(f"错误：写出验证信息失败：{exc}", file=sys.stderr)
        raise
    print(f"验证：{len(groups)} 个参考基因组的 gene TPM/FPKM 矩阵已生成，无样本遗漏")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"错误：gene 丰度矩阵构建失败：{exc}", file=sys.stderr)
        sys.exit(1)
