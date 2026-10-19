#!/usr/bin/env bash
### LAI 批处理：短 ID、按样本执行、并发上限、成功校验和定向清理。
### 默认仅生成作业；添加 --submit 才提交。需要 Bash、Python 3 和 Slurm。
set -euo pipefail
if [[ "${1:-}" == "--worker" ]]; then
    # conda 的旧版 activate 不兼容 nounset；激活失败仍立即退出。
    set +u
    source "${LAI_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
    conda activate "${LAI_ENV:-LTR_retriever}"
    set -u
fi
exec "${PYTHON:-python3}" - "$0" "$@" <<'PY'
import argparse
import csv
import datetime
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import uuid

SELF = Path(sys.argv[1]).resolve()
ARGS = sys.argv[2:]
INDEX_EXT = ('suf', 'lcp', 'llv', 'esq', 'des', 'sds', 'ssp', 'prj', 'md5', 'tis', 'ois', 'bck', 'bwt')
DEFAULT_FINDER = os.environ.get('LAI_FINDER', os.path.expanduser('~/software/LTR_FINDER_parallel/LTR_FINDER_parallel'))

def fail(message):
    raise RuntimeError(message)

def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    os.replace(str(tmp), str(path))

def fingerprint(path):
    p = Path(path).resolve(strict=True)
    s = p.stat()
    return {'path': str(p), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}

def normalize(source, dest):
    """流式重命名，保留顺序和序列；原文件从不修改。"""
    before = fingerprint(source)
    opener = gzip.open if str(source).endswith('.gz') else open
    seen, count, bases, length, max_id = set(), 0, 0, 0, 0
    digest = hashlib.sha256()
    out = Path(dest)
    with opener(source, 'rt', encoding='utf-8') as inp, out.open('w') as fasta, \
            out.with_name('contig_id_map.tsv').open('w') as mapping:
        writer = csv.writer(mapping, delimiter='\t', lineterminator='\n')
        writer.writerow(['short_id', 'original_id', 'original_header'])
        for line in inp:
            if line.startswith('>'):
                if count and length == 0:
                    fail('Empty FASTA record')
                header = line[1:].rstrip('\r\n')
                if not header.split():
                    fail('Empty FASTA header')
                old = header.split()[0]
                if old in seen:
                    fail('Duplicate original FASTA ID: ' + old)
                seen.add(old)
                count += 1
                max_id = max(max_id, len(old))
                short = 's{:09d}'.format(count)
                if len(short) > 13:
                    fail('Too many sequences for short IDs')
                writer.writerow([short, old, header])
                fasta.write('>' + short + '\n')
                digest.update(b'>')
                length = 0
            elif line.strip():
                if count == 0:
                    fail('Sequence before first FASTA header')
                seq = ''.join(line.split())
                if not re.fullmatch('[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+', seq):
                    fail('Unsupported FASTA sequence characters')
                fasta.write(seq + '\n')
                digest.update(seq.encode('ascii'))
                bases += len(seq)
                length += len(seq)
    if count == 0 or length == 0:
        fail('Empty FASTA or final record')
    if fingerprint(source) != before:
        fail('Input changed while normalizing')
    result = dict(input=before, contigs=count, bases=bases, max_original_id=max_id,
                  sequence_sha256=digest.hexdigest(), id_format='s%09d')
    atomic_json(out.with_name('input_provenance.json'), result)
    return result

def validate_result(directory, prefix, log='step3.out'):
    d = Path(directory)
    for suffix in ('.fa.out.LAI', '.fa.pass.list', '.fa.LTRlib.fa', '.fa.out'):
        p = d / (prefix + suffix)
        if not p.is_file() or p.stat().st_size == 0:
            fail('Missing/empty final result: ' + str(p))
    logfile = d / log
    if not logfile.is_file() or 'All analyses were finished!' not in logfile.read_text(errors='replace'):
        fail('No LTR_retriever completion message: ' + str(logfile))
    with (d / (prefix + '.fa.out.LAI')).open() as handle:
        header = handle.readline().split()
        if header != ['Chr', 'From', 'To', 'Intact', 'Total', 'raw_LAI', 'LAI']:
            fail('Unexpected LAI columns')
        rows = [line.split() for line in handle if line.startswith('whole_genome\t')]
    if len(rows) != 1 or len(rows[0]) != 7:
        fail('Missing/ambiguous whole_genome LAI row')
    row = rows[0]
    if not all(math.isfinite(float(x)) for x in row[1:]) or float(row[2]) <= 0:
        fail('Invalid whole_genome LAI values')
    return dict(zip(header, row))

def cleanup_candidates(directory, prefix, phase='final', normalized=False):
    d = Path(directory)
    names = {prefix + '.' + ext for ext in INDEX_EXT}
    names.add(prefix + '.mod.ltrTE.fa')
    if phase == 'final':
        g = prefix + '.fa'
        names.update([g + '.masked', g + '.defalse',
                      g + '.pass.list.bed-' + g + '.out.bed',
                      g + '.pass.list.bed-' + g + '.out.bed.homo',
                      g + '.pass.list.bed.cmb',
                      g + '.out.bed-' + g + '.pass.list.bed.cmb'])
        for p in d.iterdir():
            if re.fullmatch(re.escape(g) + r'\.out\.[0-9]+\.LAI\.LTR\.ava\.out', p.name):
                names.add(p.name)
        if normalized:
            names.add(g)
    # 不递归、不跟随软链接；保留全部未列入清单的文件。
    return [d / n for n in sorted(names) if (d / n).is_file() and not (d / n).is_symlink()]

def cleanup(directory, prefix, apply=False, phase='final', normalized=False):
    d = Path(directory)
    paths = cleanup_candidates(d, prefix, phase, normalized)
    rows = [{'file': p.name, 'bytes': p.stat().st_size} for p in paths]
    if apply:
        stamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex[:8]
        manifest = d / ('cleanup.' + phase + '.' + stamp + '.json')
        atomic_json(manifest, {'phase': phase, 'state': 'planned', 'files': rows})
        for p in paths:
            p.unlink()
        atomic_json(manifest, {'phase': phase, 'state': 'completed', 'files': rows})
    return rows

def run_command(command, directory, step, append=False):
    d = Path(directory)
    with (d / (step + '.out')).open('a' if append else 'w') as out, \
            (d / (step + '.err')).open('a' if append else 'w') as err:
        print(datetime.datetime.now().isoformat(), ' '.join(shlex.quote(x) for x in command), flush=True)
        subprocess.run(command, cwd=str(d), stdout=out, stderr=err, check=True)

def valid_done(sample_dir, input_fp):
    done = Path(sample_dir) / 'complete.json'
    if not done.is_file():
        return False
    record = json.loads(done.read_text())
    if record['input'] != input_fp:
        return False
    validate_result(record['directory'], 'genome')
    return True

def worker(manifest_path, index):
    conf = json.loads(Path(manifest_path).read_text())
    sample = conf['samples'][int(index)]
    sample_dir = Path(conf['out']) / sample['id']
    sample_dir.mkdir(parents=True, exist_ok=True)
    with (sample_dir / '.worker.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if valid_done(sample_dir, sample['input']):
            print('Already completed:', sample['id'])
            return
        if fingerprint(sample['input']['path']) != sample['input']:
            fail('Input changed since job generation')
        d = sample_dir / ('run_' + conf['run_id'] + '_' + uuid.uuid4().hex[:8])
        d.mkdir()
        try:
            atomic_json(d / 'run_config.json', conf)
            normalize(sample['input']['path'], d / 'genome.fa')
            run_command(['perl', conf['finder'], '-seq', 'genome.fa', '-threads', str(conf['finder_threads']),
                         '-harvest_out', '-size', '1000000', '-time', '300', '-w', '2', '-C',
                         '-D', '15000', '-d', '1000', '-L', '7000', '-l', '100', '-p', '20', '-M', '0.85'], d, 'step1')
            run_command(['gt', 'suffixerator', '-db', 'genome.fa', '-indexname', 'genome',
                         '-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna'], d, 'step2')
            run_command(['gt', 'ltrharvest', '-index', 'genome', '-similar', '85', '-vic', '10', '-seed', '20',
                         '-seqids', 'yes', '-minlenltr', '100', '-maxlenltr', '7000', '-mintsd', '4',
                         '-maxtsd', '6', '-motif', 'TGCA', '-motifmis', '1', '-out', 'genome.mod.ltrTE.fa'], d, 'harvest')
            os.replace(str(d / 'harvest.out'), str(d / 'genome.harvest.scn'))
            with (d / 'genome.rawLTR.scn').open('wb') as out:
                for name in ('genome.harvest.scn', 'genome.fa.finder.combine.scn'):
                    p = d / name
                    if not p.is_file() or p.stat().st_size == 0:
                        fail('Missing candidate file: ' + name)
                    with p.open('rb') as inp:
                        shutil.copyfileobj(inp, out)
            # LTR_retriever 只读 genome.fa 和 rawLTR.scn；立即释放索引空间。
            if conf['cleanup']:
                cleanup(d, 'genome', apply=True, phase='harvest')
            run_command(['LTR_retriever', '-genome', 'genome.fa', '-inharvest', 'genome.rawLTR.scn',
                         '-threads', str(conf['retriever_threads'])], d, 'step3')
            result = validate_result(d, 'genome')
            atomic_json(d / 'SUCCESS.json', {'input': sample['input'], 'LAI': result})
            if conf['cleanup']:
                cleanup(d, 'genome', apply=True, normalized=True)
            atomic_json(sample_dir / 'complete.json', {'input': sample['input'], 'directory': str(d), 'LAI': result})
            print('SUCCESS', sample['id'], json.dumps(result), flush=True)
        except BaseException as exc:
            atomic_json(d / 'FAILED.json', {'error': str(exc), 'time': time.time()})
            raise

def query_jobs(job_ids=None):
    command = ['squeue', '-h', '-u', os.environ.get('USER', ''), '-o', '%F']
    # 查询用户全部作业，避免已消失的 job ID 让某些 Slurm 版本报错。
    result = subprocess.run(command, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    active = set(result.stdout.split())
    return active if job_ids is None else active.intersection(str(j) for j in job_ids)

def prepare(opts):
    table = Path(opts.samples).resolve(strict=True)
    out = Path(opts.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with (out / '.submit.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state_path = out / 'submitted_jobs.json'
        state = json.loads(state_path.read_text()) if state_path.exists() else []
        if opts.submit and state:
            active = query_jobs([s['job_id'] for s in state])
            if active:
                fail('This output root already has active arrays: ' + ','.join(sorted(active)))
        samples, seen = [], set()
        for n, line in enumerate(table.read_text().splitlines(), 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            row = line.split()
            if len(row) != 2:
                fail('Expected sample_id and FASTA path at line ' + str(n))
            sid, source = row
            if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in seen:
                fail('Unsafe/duplicate sample ID: ' + sid)
            seen.add(sid)
            if opts.only and sid not in opts.only:
                continue
            path = Path(source)
            if not path.is_absolute():
                path = table.parent / path
            fp = fingerprint(path)
            if fp['size'] == 0:
                fail('Empty input: ' + str(path))
            if valid_done(out / sid, fp):
                print('SKIP completed', sid)
                continue
            samples.append({'id': sid, 'input': fp})
        if opts.only and not set(opts.only).issubset(seen):
            fail('Unknown --only sample')
        if not samples:
            print('No samples to run')
            return
        finder = str(Path(opts.finder).resolve(strict=True))
        run_id = datetime.datetime.now().strftime('%Y%m%dT%H%M%S') + '_' + uuid.uuid4().hex[:8]
        batch = out / '_batches' / run_id
        batch.mkdir(parents=True)
        snapshot = batch / SELF.name
        shutil.copy2(str(SELF), str(snapshot))
        shutil.copy2(str(table), str(batch / 'samples.original.txt'))
        conf = {'out': str(out), 'run_id': run_id, 'samples': samples, 'finder': finder,
                'finder_threads': opts.finder_threads, 'retriever_threads': opts.retriever_threads,
                'cleanup': not opts.keep_intermediates}
        manifest = batch / 'manifest.json'
        atomic_json(manifest, conf)
        cpus = max(opts.finder_threads, opts.retriever_threads, 4)
        # 每个数组元素只有一个样本，三个步骤串行，不另行 sbatch。
        header = ['#!/usr/bin/env bash', '#SBATCH -J LAI_t6', '#SBATCH -p ' + opts.partition,
                  '#SBATCH -N 1', '#SBATCH --ntasks=1', '#SBATCH --cpus-per-task=' + str(cpus),
                  '#SBATCH --array=0-{}%{}'.format(len(samples) - 1, opts.max_running),
                  '#SBATCH -o ' + str(batch / '%A_%a.out'), '#SBATCH -e ' + str(batch / '%A_%a.err')]
        if opts.mem:
            header.append('#SBATCH --mem=' + opts.mem)
        header += ['set -euo pipefail', 'export LAI_CONDA_SH=' + shlex.quote(opts.conda_sh),
                   'export LAI_ENV=' + shlex.quote(opts.env),
                   'exec bash {} --worker {} "$SLURM_ARRAY_TASK_ID"'.format(shlex.quote(str(snapshot)), shlex.quote(str(manifest)))]
        job = batch / 'run.slurm'
        job.write_text('\n'.join(header) + '\n')
        subprocess.run(['bash', '-n', str(job)], check=True)
        print('Prepared {} samples; max concurrent {}; CPUs/job {}'.format(len(samples), opts.max_running, cpus))
        print('BATCH:', str(batch))
        if opts.submit:
            response = subprocess.run(['sbatch', '--parsable', str(job)], universal_newlines=True, stdout=subprocess.PIPE, check=True).stdout.strip()
            jid = response.split(';')[0]
            if not jid.isdigit():
                fail('Unexpected sbatch response: ' + response)
            state.append({'job_id': jid, 'batch': str(batch)})
            atomic_json(state_path, state)
            print('SUBMITTED:', jid)
        else:
            print('DRY RUN: no jobs submitted. Add --submit to submit through this script.')

def old_cleanup(opts):
    root = Path(opts.cleanup_existing).resolve(strict=True)
    # 旧流程无锁，执行清理前保守要求该账号没有作业；预览不受限。
    if opts.apply and query_jobs():
        fail('Existing-result cleanup requires no queued/running jobs for this account')
    reports = []
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.is_symlink() or (opts.only and d.name not in opts.only):
            continue
        try:
            result = validate_result(d, d.name)
        except (RuntimeError, ValueError) as exc:
            reports.append({'sample': d.name, 'status': 'SKIP', 'reason': str(exc)})
            continue
        rows = cleanup(d, d.name, apply=opts.apply)
        reports.append({'sample': d.name, 'status': 'APPLIED' if opts.apply else 'PREVIEW',
                        'LAI': result, 'bytes': sum(r['bytes'] for r in rows), 'files': rows})
    print(json.dumps(reports, ensure_ascii=False, indent=2))

def main():
    if ARGS and ARGS[0] == '--worker':
        worker(ARGS[1], ARGS[2])
        return
    if ARGS and ARGS[0] == '--normalize':
        normalize(ARGS[1], ARGS[2])
        return
    parser = argparse.ArgumentParser(description='LAI batch workflow; default: prepare only, no submission')
    parser.add_argument('samples', nargs='?')
    parser.add_argument('--out', default='./results_LTR_t6')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--submit', action='store_true')
    modes.add_argument('--no-submit', action='store_true')
    parser.add_argument('--max-running', type=int, default=50)
    parser.add_argument('--finder-threads', type=int, default=30)
    parser.add_argument('--retriever-threads', type=int, default=5)
    parser.add_argument('--partition', default='x86_64_2cpu')
    parser.add_argument('--mem', help='Optional Slurm memory request, e.g. 64G; no default --mem')
    parser.add_argument('--finder', default=DEFAULT_FINDER)
    parser.add_argument('--conda-sh', default=os.path.expanduser('~/miniconda3/etc/profile.d/conda.sh'))
    parser.add_argument('--env', default='LTR_retriever')
    parser.add_argument('--only', nargs='+', help='Restrict to listed sample IDs')
    parser.add_argument('--keep-intermediates', action='store_true')
    parser.add_argument('--cleanup-existing', metavar='OLD_RESULTS_DIR')
    parser.add_argument('--apply', action='store_true', help='Apply existing-result cleanup (default preview)')
    opts = parser.parse_args(ARGS)
    if not 1 <= opts.max_running <= 50 or min(opts.finder_threads, opts.retriever_threads) < 1:
        parser.error('Concurrency must be 1..50 and threads must be positive')
    for value in (opts.partition, opts.mem or 'none'):
        if not re.fullmatch('[A-Za-z0-9_.-]+', value):
            parser.error('Unsafe partition/memory value')
    if any(c.isspace() for c in str(Path(opts.out).resolve())):
        parser.error('Output path must not contain whitespace (Slurm directive limitation)')
    if opts.cleanup_existing:
        old_cleanup(opts)
    elif opts.samples and not opts.apply:
        prepare(opts)
    else:
        parser.error('Provide sample table, or --cleanup-existing [--apply]')

try:
    main()
except (RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as exc:
    print('ERROR:', exc, file=sys.stderr)
    sys.exit(1)
PY
