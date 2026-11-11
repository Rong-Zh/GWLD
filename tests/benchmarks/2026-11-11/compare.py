"""Compare all five engine methods, including complete constant-column cases."""
import argparse
import csv
import hashlib
import json
import os
import pathlib
import statistics
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=pathlib.Path, required=True)
parser.add_argument('--after', type=pathlib.Path, required=True)
parser.add_argument('--output', type=pathlib.Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
binaries = {'before': args.before.resolve(), 'after': args.after.resolve()}
methods = ['MI', 'RMI', 'r2', 'D', "D'"]
records, checks = [], []
environment = dict(os.environ, OMP_DYNAMIC='FALSE')
for samples, variants in [(512, 1024), (4096, 128), (65536, 128)]:
    for constant in [0, 50]:
        for missing in [0, 25]:
            times = {v: {m: [] for m in methods} for v in binaries}
            for repeat in range(3):
                for v in (list(binaries) if repeat % 2 == 0 else list(reversed(binaries))):
                    prefix = str((args.output/v).resolve()) if repeat == 0 else '-'
                    result = subprocess.check_output([str(binaries[v]), str(samples), str(variants),
                        str(missing), '1', prefix, 'uniform', str(constant)], text=True, env=environment)
                    for row in csv.DictReader(result.splitlines(), delimiter='\t'):
                        if int(row['constant_percent']) != constant:
                            raise RuntimeError('Incorrect constant-column workload')
                        times[v][row['method']].append(float(row['median_seconds']))
            for method in methods:
                left = (args.output/('before-'+method+'.bin')).read_bytes()
                right = (args.output/('after-'+method+'.bin')).read_bytes()
                if left != right:
                    raise RuntimeError('Matrix differs: '+str((samples, variants, constant, missing, method)))
                before, after = (statistics.median(times[v][method]) for v in binaries)
                records.append([method, samples, variants, constant, missing, 1, before, after, before/after])
                checks.append({'method': method, 'samples': samples, 'variants': variants,
                    'constant_percent': constant, 'missing_percent': missing,
                    'matrix_sha256': hashlib.sha256(left).hexdigest(), 'byte_identical': True})
            print('Checked:', samples, variants, 'constant', constant, 'missing', missing, flush=True)
with (args.output/'engine.tsv').open('w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['method', 'samples', 'variants', 'constant_percent', 'missing_percent',
                     'threads', 'before_seconds', 'after_seconds', 'speedup'])
    writer.writerows(records)
(args.output/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
print('All', len(checks), 'matrix comparisons passed')
