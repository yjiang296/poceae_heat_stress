#!/usr/bin/env python3
"""Summarize mean FIRE profiles; use --config on HPC or --snapshot locally."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def collect(config):
    text = config.read_text()
    samples = []
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        name, bw, bed, directory = line.split('\t')
        out = Path(directory)
        profile = out / f'{name}.TSS.profile.tsv'
        raw = profile.read_bytes()
        with gzip.open(out / f'{name}.TSS.matrix.gz', 'rt') as handle:
            header = json.loads(handle.readline().lstrip('@'))
        samples.append(dict(sample=name, bigwig=bw, tss_bed=bed, output_dir=directory,
                            profile_path=str(profile), profile_text=raw.decode(),
                            profile_sha256=hashlib.sha256(raw).hexdigest(), matrix_header=header))
    return dict(config_path=str(config.resolve()), config_text=text, samples=samples)


def summarize(source):
    samples = source['samples']
    if not samples or len({r['sample'] for r in samples}) != len(samples):
        raise ValueError('Empty input or duplicate sample names')
    result, long_rows = [], []
    for record in samples:
        name = record['sample']
        if 'error' in record:
            raise ValueError(f'{name}: {record["error"]}')
        h = record['matrix_header']
        required = {'upstream': [2000], 'downstream': [2000], 'body': [0],
                    'bin size': [20], 'ref point': ['TSS'], 'sample_boundaries': [0, 200],
                    'sample_labels': [name], 'bin avg type': 'mean',
                    'missing data as zero': True, 'scale': 1, 'skip zeros': False}
        for key, value in required.items():
            if h.get(key) != value:
                raise ValueError(f'{name}: unexpected matrix parameter {key}={h.get(key)}')
        lines = list(csv.reader(record['profile_text'].splitlines(), delimiter='\t'))
        if len(lines) != 3 or lines[1][:2] != ['bins', ''] or lines[2][0] != name:
            raise ValueError(f'{name}: expected one sample and one region-group profile')
        y = np.array(lines[2][2:], dtype=float)
        if len(y) != 200 or not np.isfinite(y).all() or (y < 0).any():
            raise ValueError(f'{name}: invalid coverage values')
        if hashlib.sha256(record['profile_text'].encode()).hexdigest() != record['profile_sha256']:
            raise ValueError(f'{name}: profile checksum mismatch')
        starts = np.arange(-2000, 2000, 20)
        ends = starts + 20
        background_mask = (ends <= -1500) | (starts >= 1500)
        background = float(y[background_mask].mean())
        if background <= 0:
            raise ValueError(f'{name}: background must be positive')
        peak_i = int(np.argmax(y))
        boundaries = h['group_boundaries']
        if len(boundaries) != 2 or boundaries[0] != 0:
            raise ValueError(f'{name}: expected one region group')
        result.append(dict(sample_name=name, max_mean_FIRE_coverage=float(y[peak_i]),
                           flank_mean_FIRE_coverage=background,
                           peak_enrichment_fold=float(y[peak_i] / background),
                           peak_bin_start_bp=int(starts[peak_i]), peak_bin_end_bp=int(ends[peak_i]),
                           peak_bin_center_bp=int(starts[peak_i]+10),
                           n_TSS=boundaries[-1], n_profile_bins=len(y), n_background_bins=50,
                           profile_path=record['profile_path'], profile_sha256=record['profile_sha256']))
        for i in range(len(y)):
            long_rows.append(dict(sample_name=name, bin_start_bp=int(starts[i]),
                                  bin_end_bp=int(ends[i]), bin_center_bp=int(starts[i]+10),
                                  mean_FIRE_coverage=float(y[i]), enrichment_fold=float(y[i]/background),
                                  is_background=int(background_mask[i])))
    ranks = {r['sample_name']: i+1 for i, r in enumerate(sorted(result, key=lambda r: -r['peak_enrichment_fold']))}
    for r in result:
        r['rank_by_enrichment'] = ranks[r['sample_name']]
    return result, long_rows


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def draw(rows, out):
    values = np.array([r['peak_enrichment_fold'] for r in rows])
    q1, median, q3 = np.percentile(values, [25, 50, 75])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(figsize=(6.7, 6.5))
    fig.subplots_adjust(left=.17, right=.95, bottom=.26, top=.79)
    ax.boxplot(values, positions=[0], widths=.30, patch_artist=True, showfliers=False,
               boxprops={'facecolor': '#B8DCD8', 'edgecolor': '#22665F', 'linewidth': 1.5},
               medianprops={'color': '#143D39', 'linewidth': 2, 'zorder': 4},
               whiskerprops={'color': '#22665F', 'linewidth': 1.3},
               capprops={'color': '#22665F', 'linewidth': 1.3})
    # Fixed random seed only (not a sample count); all samples are drawn.
    x = np.random.default_rng(142).uniform(-.12, .12, len(values))
    ax.scatter(x, values, s=23, color='#297F77', alpha=.65, linewidths=.35, edgecolors='white', zorder=3)
    ax.axhline(1, color='#A5A5A5', linewidth=1, linestyle=(0, (4, 4)), zorder=0)
    ax.text(-.39, 1.06, 'Flank background = 1', fontsize=9, color='#777777')
    ax.set(xlim=(-.43, .43), ylim=(0, max(values)*1.10), xticks=[0],
           xticklabels=[f'All samples (n = {len(rows)})'], ylabel='Peak / flank background (fold)')
    ax.yaxis.grid(True, color='#E7EBEC', linewidth=.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis='x', length=0, pad=12)
    fig.suptitle('Peak TSS enrichment across samples', x=.56, y=.94, fontsize=16, fontweight='bold')
    fig.text(.56, .872, f'Median {median:.2f}×  |  IQR {q1:.2f}–{q3:.2f}×\nRange {min(values):.2f}–{max(values):.2f}×',
             ha='center', va='center', color='#435057', fontsize=11, linespacing=1.6)
    fig.text(.17, .055, 'Each dot = one sample. Box = median and interquartile range;\nwhiskers extend to the most extreme values within 1.5 × IQR.\nPeak: maximum 20-bp bin within TSS ±2 kb.\nBackground: mean of −2 to −1.5 kb and +1.5 to +2 kb.',
             fontsize=9, color='#56616A', linespacing=1.45)
    for ext in ('png', 'pdf'):
        fig.savefig(out/f'tss_peak_enrichment_{len(rows)}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    return dict(n_samples=len(rows), median=float(median), q1=float(q1), q3=float(q3),
                minimum=min(values), maximum=max(values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--config', type=Path, help='Four-column sample config on the server')
    group.add_argument('--snapshot', type=Path, help='Saved source_profiles.json; no original paths required')
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    source = collect(args.config) if args.config else json.loads(args.snapshot.read_text())
    rows, long_rows = summarize(source)
    n = len(rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir/'source_profiles.json').write_text(json.dumps(source, ensure_ascii=False)+'\n')
    (args.outdir/f'sample_info_config.{n}.txt').write_text(source['config_text'])
    write_tsv(args.outdir/f'tss_peak_summary.{n}.tsv', rows)
    write_tsv(args.outdir/f'tss_profiles.{n}.long.tsv', long_rows)
    stats = draw(rows, args.outdir)
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
