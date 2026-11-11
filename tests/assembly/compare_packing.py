"""Compare packing and engine binaries compiled against two header revisions."""
import argparse
import csv
import pathlib
import statistics
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=pathlib.Path, required=True,
                    help='Directory containing packing_benchmark and gwld_kernel_benchmark')
parser.add_argument('--after', type=pathlib.Path, required=True)
parser.add_argument('--output', type=pathlib.Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
versions = {'before': args.before.resolve(), 'after': args.after.resolve()}

def write(name, header, rows):
    with (args.output/name).open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)

rows = []
for samples, variants in [(512, 384), (4096, 128), (65536, 128)]:
    for missing in [0, 25, 100]:
        times = {v: [] for v in versions}
        for repeat in range(5):
            for v in (list(versions) if repeat % 2 == 0 else list(reversed(versions))):
                result = subprocess.check_output([str(versions[v]/'packing_benchmark'),
                    str(samples), str(variants), str(missing), '15'], text=True)
                times[v].append(float(result.split()[0]))
        before, after = (statistics.median(times[v]) for v in versions)
        rows.append([samples, variants, missing, before, after, before/after])
write('packing.tsv', ['samples', 'variants', 'missing_percent', 'before_seconds', 'after_seconds', 'speedup'], rows)

rows = []
for samples, variants in [(512, 768), (4096, 128), (65536, 128)]:
    for missing in [0, 25]:
        times = {v: {} for v in versions}
        for repeat in range(3):
            for v in (list(versions) if repeat % 2 == 0 else list(reversed(versions))):
                prefix = str((args.output/v).resolve()) if repeat == 0 else '-'
                result = subprocess.check_output([str(versions[v]/'gwld_kernel_benchmark'),
                    str(samples), str(variants), str(missing), '1', prefix], text=True)
                for row in csv.DictReader(result.splitlines(), delimiter='\t'):
                    times[v].setdefault(row['method'], []).append(float(row['median_seconds']))
        for method in ['MI', 'RMI', 'r2']:
            if (args.output/('before-'+method+'.bin')).read_bytes() != (args.output/('after-'+method+'.bin')).read_bytes():
                raise RuntimeError('Matrix mismatch: '+str((samples, variants, missing, method)))
            before, after = (statistics.median(times[v][method]) for v in versions)
            rows.append([method, samples, variants, missing, 1, before, after, before/after])
write('engine.tsv', ['method', 'samples', 'variants', 'missing_percent', 'threads',
                     'before_seconds', 'after_seconds', 'speedup'], rows)
print('Timing comparisons and byte-identical engine matrices passed:', args.output)
