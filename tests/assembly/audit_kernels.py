"""Compile reproducible probes and check the emitted counting instructions."""
import argparse
import json
import pathlib
import platform
import re
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cxx', default='c++')
parser.add_argument('--arch', choices=['native', 'arm64', 'x86_64'], default='native',
                    help='Cross-architecture selection uses the macOS -arch flag')
parser.add_argument('--output', default='build/assembly')
args = parser.parse_args()
root = pathlib.Path(__file__).resolve().parents[2]
out = pathlib.Path(args.output).resolve()
out.mkdir(parents=True, exist_ok=True)
command = [args.cxx, '-std=c++17', '-O3', '-S', '-I', str(root/'GWLD-R/inst/include'),
           str(root/'tests/assembly/probes.cpp'), '-o', str(out/'kernels.s')]
if platform.system() == 'Darwin':
    sdk = subprocess.check_output(['xcrun', '--sdk', 'macosx', '--show-sdk-path'],
                                  text=True, stderr=subprocess.PIPE).strip()
    command += ['-D_REENTRANT', '-isysroot', sdk]
    if args.arch != 'native':
        command += ['-arch', args.arch]
elif args.arch != 'native':
    parser.error('--arch is supported only by macOS compilers; choose a native/cross compiler with --cxx')
compiled = subprocess.run(command, capture_output=True, text=True)
(out/'compiler.log').write_text(compiled.stdout+compiled.stderr)
if compiled.returncode:
    raise SystemExit(compiled.stderr)
assembly = (out/'kernels.s').read_text()
functions = {}
for match in re.finditer(r'^([^\s;#][^\n]*):[^\n]*\n', assembly, re.M):
    name = match.group(1)
    if not any(key in name for key in ('count_packed_', 'count_four_', 'gwld_neon_four_probe')):
        continue
    body = assembly[match.end():].split('.cfi_endproc', 1)[0].split('.seh_endproc', 1)[0]
    functions[name] = {
        'scalar_popcnt': len(re.findall(r'^\s+popcnt\w*\s', body, re.M)),
        'vector_popcnt': len(re.findall(r'^\s+vpopcnt\w*\s', body, re.M)),
        'avx2_byte_shuffle': len(re.findall(r'^\s+vpshufb\s', body, re.M)),
        'neon_cnt': len(re.findall(r'^\s+cnt(?:\.\w+)?\s', body, re.M)),
        'neon_pair_accumulate': len(re.findall(r'^\s+uadalp(?:\.\w+)?\s', body, re.M)),
    }
# Check instructions only for the architecture that the compiler emitted.
if any('count_packed_avx2' in name for name in functions):
    for key, instruction in [('count_packed_avx2', 'avx2_byte_shuffle'),
                             ('count_packed_avx512', 'vector_popcnt')]:
        selected = [v for k, v in functions.items() if key in k]
        if selected and not any(v[instruction] for v in selected):
            raise SystemExit('Expected vector instruction missing from '+key)
    portable = [v for k, v in functions.items() if 'count_packed_portable' in k]
    if any(v['scalar_popcnt'] or v['vector_popcnt'] for v in portable):
        raise SystemExit('Portable fallback contains POPCNT; check compiler baseline flags')
else:
    neon = [v for k, v in functions.items() if 'gwld_neon_four_probe' in k]
    if neon and not any(v['neon_cnt'] and v['neon_pair_accumulate'] for v in neon):
        raise SystemExit('Expected NEON counting/accumulation instructions missing')
report = {'compiler': subprocess.check_output([args.cxx, '--version'], text=True).splitlines()[0],
          'command': command, 'host': platform.machine(), 'arch': args.arch, 'functions': functions}
(out/'instructions.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(functions, indent=2))
print('Assembly checks passed:', out)
