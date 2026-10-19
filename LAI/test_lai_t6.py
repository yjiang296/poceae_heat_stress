import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().with_name('submit_ltr_pipeline.sh')
source = SCRIPT.read_text().split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
source = source.rsplit('\ntry:\n    main()', 1)[0]
with patch.object(sys, 'argv', ['-', str(SCRIPT)]):
    api = {}
    exec(compile(source, str(SCRIPT), 'exec'), api)

class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fa = self.root/'input.fa'
        self.fa.write_text('>long_contig_name_one description\nACGTNN\n>long_contig_name_two\nacgt\n')

    def final_files(self, d, prefix='genome'):
        for suffix in ('.fa.pass.list', '.fa.LTRlib.fa', '.fa.out'):
            (d/(prefix+suffix)).write_text('data\n')
        (d/(prefix+'.fa.out.LAI')).write_text('Chr\tFrom\tTo\tIntact\tTotal\traw_LAI\tLAI\nwhole_genome\t1\t10\t0.2\t0.8\t25\t20\n')
        (d/'step3.out').write_text('All analyses were finished!\n')

    def test_rename_preserves_bases_order_mapping_and_input(self):
        before = self.fa.read_bytes()
        out = self.root/'genome.fa'
        report = api['normalize'](self.fa, out)
        self.assertEqual(out.read_text(), '>s000000001\nACGTNN\n>s000000002\nacgt\n')
        self.assertEqual(self.fa.read_bytes(), before)
        self.assertEqual(report['bases'], 10)
        self.assertIn('long_contig_name_one description', (self.root/'contig_id_map.tsv').read_text())

    def test_duplicates_and_empty_records_fail(self):
        for content in ('>x\nAC\n>x\nGT\n', '>x\n>y\nAC\n', '>x\n'):
            self.fa.write_text(content)
            with self.assertRaises(RuntimeError):
                api['normalize'](self.fa, self.root/'genome.fa')

    def test_cleanup_is_exact_and_does_not_follow_symlinks(self):
        self.final_files(self.root)
        (self.root/'genome.suf').write_bytes(b'index')
        (self.root/'genome.fa.masked').symlink_to(self.fa)
        (self.root/'genome.fa.out.20000.LAI.LTR.ava.out').write_bytes(b'temp')
        (self.root/'unknown.large').write_text('retain')
        before = self.fa.read_bytes()
        api['validate_result'](self.root, 'genome')
        preview = api['cleanup'](self.root, 'genome')
        self.assertEqual(sum(x['bytes'] for x in preview), 9)
        self.assertTrue((self.root/'genome.suf').exists())
        api['cleanup'](self.root, 'genome', apply=True)
        self.assertFalse((self.root/'genome.suf').exists())
        self.assertEqual(self.fa.read_bytes(), before)
        self.assertTrue((self.root/'genome.fa.out').exists())
        self.assertTrue((self.root/'unknown.large').exists())

    def test_invalid_lai_or_no_completion_fails(self):
        self.final_files(self.root)
        lai = self.root/'genome.fa.out.LAI'
        lai.write_text(lai.read_text().replace('\t20\n', '\tnan\n'))
        with self.assertRaises(RuntimeError): api['validate_result'](self.root, 'genome')
        self.final_files(self.root)
        (self.root/'step3.out').write_text('Step 3 Runtime: 16s')
        with self.assertRaises(RuntimeError): api['validate_result'](self.root, 'genome')

    def prepare_worker(self):
        cfg = dict(out=str(self.root/'results'), run_id='test', cleanup=True, finder='/fake/finder',
                   finder_threads=30, retriever_threads=5,
                   samples=[{'id':'B73', 'input':api['fingerprint'](self.fa)}])
        manifest = self.root/'manifest.json'
        manifest.write_text(json.dumps(cfg))
        return manifest

    def mock_command(self, cmd, d, step, append=False):
        d = Path(d)
        (d/(step+'.out')).write_text('log\n')
        (d/(step+'.err')).write_text('')
        if step == 'step1': (d/'genome.fa.finder.combine.scn').write_text('finder\n')
        if step == 'step2': (d/'genome.suf').write_text('index')
        if step == 'harvest':
            (d/'harvest.out').write_text('harvest\n')
            (d/'genome.mod.ltrTE.fa').write_text('candidate\n')
        if step == 'step3':
            self.assertFalse((d/'genome.suf').exists())
            self.final_files(d)
            (d/'genome.fa.masked').write_text('mask')

    def test_success_worker_cleanup_skip_and_retained_results(self):
        manifest = self.prepare_worker()
        with patch.dict(api, {'run_command':self.mock_command}):
            api['worker'](manifest, 0)
        root = self.root/'results/B73'
        done = json.loads((root/'complete.json').read_text())
        d = Path(done['directory'])
        self.assertFalse((d/'genome.fa').exists())
        self.assertFalse((d/'genome.fa.masked').exists())
        self.assertTrue((d/'genome.fa.pass.list').exists())
        self.assertTrue((d/'contig_id_map.tsv').exists())
        with patch.dict(api, {'run_command':lambda *a, **k: self.fail('should skip')}):
            api['worker'](manifest, 0)

    def test_failed_command_preserves_evidence_and_no_success(self):
        manifest = self.prepare_worker()
        def broken(cmd, d, step, append=False):
            if step == 'harvest': raise subprocess.CalledProcessError(2, cmd)
            self.mock_command(cmd, d, step)
        with patch.dict(api, {'run_command':broken}):
            with self.assertRaises(subprocess.CalledProcessError): api['worker'](manifest, 0)
        root = self.root/'results/B73'
        self.assertFalse((root/'complete.json').exists())
        d = next(root.glob('run_*'))
        self.assertTrue((d/'FAILED.json').exists())
        self.assertTrue((d/'genome.suf').exists())
        self.assertTrue((d/'genome.fa').exists())

    def test_zero_exit_incomplete_retriever_does_not_clean(self):
        manifest = self.prepare_worker()
        def incomplete(cmd, d, step, append=False):
            if step == 'step3':
                (Path(d)/'step3.out').write_text('failed internally')
            else: self.mock_command(cmd, d, step)
        with patch.dict(api, {'run_command':incomplete}):
            with self.assertRaises(RuntimeError): api['worker'](manifest, 0)
        root = self.root/'results/B73'
        self.assertFalse((root/'complete.json').exists())
        self.assertTrue(next(root.glob('run_*/genome.fa')).exists())

    def test_dry_run_array_limit_and_no_submission(self):
        table = self.root/'samples.txt'
        table.write_text('B73\t'+str(self.fa)+'\n')
        finder = self.root/'finder.pl'
        finder.write_text('# mock\n')
        cmd = ['bash', str(SCRIPT), str(table), '--out', str(self.root/'results'), '--finder', str(finder)]
        result = subprocess.run(cmd, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr)
        job = next((self.root/'results/_batches').glob('*/run.slurm'))
        self.assertIn('#SBATCH --array=0-0%50', job.read_text())
        self.assertNotIn('#SBATCH --mem', job.read_text())
        subprocess.run(['bash', '-n', str(job)], check=True)
        result = subprocess.run(cmd+['--max-running', '51'], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(result.returncode, 0)

    def test_real_subprocess_logging(self):
        api['run_command']([sys.executable, '-c', 'print("smoke")'], self.root, 'smoke')
        self.assertEqual((self.root/'smoke.out').read_text().strip(), 'smoke')

    def test_query_uses_array_parent_ids(self):
        result = subprocess.CompletedProcess([], 0, stdout='123\n123\n999\n', stderr='')
        with patch('subprocess.run', return_value=result) as run:
            self.assertEqual(api['query_jobs'](['123']), {'123'})
            self.assertEqual(run.call_args[0][0][-1], '%F')

    def test_active_array_blocks_second_submission(self):
        table = self.root/'samples.txt'
        table.write_text('B73\t'+str(self.fa)+'\n')
        out = self.root/'results'
        out.mkdir()
        (out/'submitted_jobs.json').write_text('[{"job_id":"123"}]')
        from argparse import Namespace
        opts = Namespace(samples=str(table), out=str(out), submit=True)
        with patch.dict(api, {'query_jobs':lambda ids: {'123'}}):
            with self.assertRaisesRegex(RuntimeError, 'active arrays'): api['prepare'](opts)

if __name__ == '__main__': unittest.main(verbosity=2)
