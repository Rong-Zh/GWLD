# Windows / WSL measurements, 2026-11-09

See [the performance report](../../../docs/performance.md#windows-and-wsl-on-the-same-x86-64-host) for interpreted tables, timing scope and limitations. Production source is `bca4f90`; the original calculation reference is `c8079ab`. Only documentation and the R benchmark's source-loading/batched-timing options were changed while recording this run.

## Raw evidence

The `windows/` and `wsl/` directories contain:

| File | Scope |
| --- | --- |
| `v1-v2-512x384-t1.tsv` | Original versus current full matrix, 25% missing, one thread. |
| `v1-v2-4096x256-t4.tsv` | Original versus current full matrix, 1% missing, four threads. |
| `scalar-vs-packed.tsv` | Current formulas and solver in both execution paths. |
| `scaling-512x2048-t{1,4,12}.tsv` | Current engine, 25% missing, thread scaling. |
| `vcf-read.tsv` | Legacy stream parser versus HTSlib, three alternating executions. |
| `io-end-to-end.csv` | Full `summarize` command, including process startup and file output. |
| `encoding.tsv` | R encoder, three alternating batches of 20 calls; median seconds per call. |
| `r-environment.txt` | R/package versions and RNG settings. |
| `*-test.txt` / `ctest.txt` | Correctness-test output. |

[Input hashes](input-files.csv) identify the exact shared files. [Matrix validation](cross-platform-values.csv) records all-pair Windows/WSL comparisons, and [frequency-output validation](cross-platform-io.txt) records equality across the 16 CLI runs after CRLF/LF normalization. Per-platform output hashes in the IO CSVs retain their original line endings and therefore differ between platforms.

No generated matrices, large VCFs, installed dependencies or executables are committed. The full 2048×2048 binary matrices used for validation were written under the ignored `build-windows-speed/` and `build-wsl-speed/` directories.

## Builds used

Both standalone builds use C++17, `-O3 -DNDEBUG -fopenmp`, without `-ffast-math` or global `-march=native`. Both disable the Armadillo wrapper, BLAS, LAPACK and SuperLU. CPU: Intel Core i5-12400, 12 logical processors. Set `OMP_DYNAMIC=FALSE` and run performance jobs sequentially.

### WSL

Ubuntu 24.04.4 LTS, Linux `6.18.58.1-microsoft-standard-WSL2`, GCC 13.3.0; bundled HTSlib 1.24 and Armadillo headers 12.6.7. The working directory and all benchmark files are on the Windows E: drive mounted at `/mnt/e`.

The actual build used locally extracted Armadillo headers:

```bash
cd /mnt/e/devel/GWLD
export OMP_DYNAMIC=FALSE
cmake -S GWLD-C++ -B build-wsl-speed \
  -DCMAKE_BUILD_TYPE=Release -DGWLD_BUILD_BENCHMARKS=ON \
  -DGWLD_ARMADILLO_INCLUDE_DIR="$PWD/.deps/armadillo/usr/include"
cmake --build build-wsl-speed -j4
ctest --test-dir build-wsl-speed --output-on-failure
```

Supply an installed Armadillo include directory instead when reproducing elsewhere. CMake also requires the bundled HTSlib compression dependencies described in the main build instructions.

### Windows native

Windows 10 Pro 22H2, build 19045.6466; GCC 14.3.0 from `C:/rtools/x86_64-w64-mingw32.static.posix/bin`, with Armadillo 15.0.3 headers in the corresponding `include` directory. The native executables were compiled directly rather than through the Linux HTSlib archive in `deps`.

The following PowerShell commands reproduce the standalone calculation builds used here. Run from the repository root, adjusting installation paths as needed:

```powershell
$env:PATH = 'C:\rtools\x86_64-w64-mingw32.static.posix\bin;' + $env:PATH
$env:OMP_DYNAMIC = 'FALSE'
$compiler = 'C:\rtools\x86_64-w64-mingw32.static.posix\bin\g++.exe'
New-Item -ItemType Directory -Force build-windows-speed | Out-Null
$flags = @('-std=c++17','-O3','-DNDEBUG','-fopenmp',
  '-DARMA_DONT_USE_WRAPPER','-DARMA_DONT_USE_BLAS',
  '-DARMA_DONT_USE_LAPACK','-DARMA_DONT_USE_SUPERLU',
  '-I','C:/rtools/x86_64-w64-mingw32.static.posix/include',
  '-I','GWLD-R/inst/include','-I','GWLD-C++/tests/legacy')
$targets = @{
  gwld_benchmark = @('GWLD-C++/tests/benchmark.cpp',
    'GWLD-C++/tests/legacy/methods.cpp','GWLD-C++/tests/legacy/utils_reference.cpp')
  gwld_engine_benchmark = @('GWLD-C++/tests/engine_benchmark.cpp')
  gwld_kernel_benchmark = @('GWLD-C++/tests/kernel_benchmark.cpp')
  gwld_core_test = @('GWLD-C++/tests/core_test.cpp')
  gwld_packed_test = @('GWLD-C++/tests/packed_test.cpp')
  gwld_simd_test = @('GWLD-C++/tests/simd_test.cpp')
}
foreach ($name in $targets.Keys) {
  & $compiler @flags @($targets[$name]) -o "build-windows-speed/$name.exe"
  if ($LASTEXITCODE -ne 0) { throw "Build failed: $name" }
}
```

For HTSlib targets, the native headers and static library came from `E:/software/R-4.6.1/library/Rhtslib` (Rhtslib 3.8.0). Its `pkgconfig("PKG_LIBS")` output supplied the static dependencies below; these paths/dependencies are specific to this installation, not general Windows package requirements:

```powershell
$htsFlags = $flags + @('-D_FILE_OFFSET_BITS=64',
  '-I','E:/software/R-4.6.1/library/Rhtslib/include',
  '-I','GWLD-C++/include','-I','GWLD-C++/deps/ArgumentsParser')
$htsLibs = @('E:/software/R-4.6.1/library/Rhtslib/usrlib/x64/libhts.a',
  '-lm','-lbz2','-llzma','-lcurl','-lnghttp2','-lpsl',
  '-lbrotlidec','-lbrotlicommon','-lsecur32','-lbcrypt','-lidn2',
  '-lunistring','-liconv','-lssl','-lcrypto','-lz','-lcrypt32',
  '-lwsock32','-lwldap32','-lssh2','-lgcrypt','-lgpg-error',
  '-lws2_32','-lzstd','-lregex')
$htsTargets = @{
  GWLD = 'GWLD-C++/src/main.cpp'
  gwld_vcf_benchmark = 'GWLD-C++/tests/vcf_benchmark.cpp'
  gwld_vcf_fixture = 'GWLD-C++/tests/vcf_fixture.cpp'
}
foreach ($name in $htsTargets.Keys) {
  & $compiler @htsFlags $htsTargets[$name] @htsLibs -o "build-windows-speed/$name.exe"
  if ($LASTEXITCODE -ne 0) { throw "Build failed: $name" }
}
```

## Calculation commands

Use these arguments on **each** platform; substitute `build-windows-speed/<name>.exe` for `build-wsl-speed/<name>` on Windows. The v1 comparison orders its arguments as `samples variants threads missing_percent`; the scaling program instead uses `samples variants missing_percent threads`.

```bash
./build-wsl-speed/gwld_benchmark 512 384 1 25
./build-wsl-speed/gwld_benchmark 4096 256 4 1
./build-wsl-speed/gwld_engine_benchmark 512 384 1 25
./build-wsl-speed/gwld_kernel_benchmark 512 2048 25 1 build-wsl-speed/cross-platform
./build-wsl-speed/gwld_kernel_benchmark 512 2048 25 4
./build-wsl-speed/gwld_kernel_benchmark 512 2048 25 12
```

The first scaling command exports MI/RMI/r² matrices after timing. The Windows run used prefix `build-windows-speed/cross-platform`. Compare every double, accepting paired NaNs and checking missing positions; our cross-platform check required a maximum absolute difference below `1e-12` for all 2,096,128 finite pairs per method. Checksums alone were not used for this check. AVX2/POPCNT dispatch was validated by running the standalone `gwld_simd_test` on each platform.

## Input generation and IO commands

The plain VCF is generated once, with LF line endings, and shared across both platforms. This Python snippet reproduces its recorded SHA256:

```python
from pathlib import Path
path = Path('.deps/vcf-benchmark.vcf')
path.parent.mkdir(parents=True, exist_ok=True)
samples, variants = 512, 12000
with path.open('w', newline='\n', buffering=1024*1024) as out:
    out.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=100000000>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
              + '\t'.join(f's{i}' for i in range(samples)) + '\n')
    suffix = '\t'.join(['0/0:20','0/1:24','1|1:30','./.:.']*(samples//4)) + '\n'
    for i in range(variants):
        out.write(f'1\t{i+1}\trs{i}\tA\tG\t.\t.\t.\tGT:DP\t' + suffix)
```

Create compressed files once on Windows, then let WSL read those same bytes through `/mnt/e`:

```powershell
./build-windows-speed/gwld_vcf_fixture.exe .deps/vcf-benchmark.vcf build-windows-speed/benchmark.vcf.gz wz
./build-windows-speed/gwld_vcf_fixture.exe .deps/vcf-benchmark.vcf build-windows-speed/benchmark.bcf wb
$inputStream = [IO.File]::OpenRead((Resolve-Path .deps/vcf-benchmark.vcf))
$outputStream = [IO.File]::Create((Join-Path $PWD 'build-windows-speed/benchmark-gzip.vcf.gz'))
$gzipStream = New-Object IO.Compression.GZipStream($outputStream,[IO.Compression.CompressionLevel]::Optimal)
try { $inputStream.CopyTo($gzipStream) }
finally { $gzipStream.Dispose(); $outputStream.Dispose(); $inputStream.Dispose() }
```

Reader benchmark, on each platform:

```bash
./build-wsl-speed/gwld_vcf_benchmark .deps/vcf-benchmark.vcf
```

For end-to-end timings, execute each `summarize` combination four times sequentially and discard the first execution. Measure wall time around a direct child-process invocation on its native platform (PowerShell `Stopwatch` on Windows; Python `perf_counter` and `subprocess.run` inside WSL). Do **not** time `wsl.exe` startup. Take the median of the remaining three times. Run with threads 1 and 4, with these four inputs:

```bash
./build-wsl-speed/GWLD summarize --vcf .deps/vcf-benchmark.vcf --threads 1 -o build-wsl-speed/summary-vcf-t1
./build-wsl-speed/GWLD summarize --vcf build-windows-speed/benchmark-gzip.vcf.gz --threads 1 -o build-wsl-speed/summary-gzip-t1
./build-wsl-speed/GWLD summarize --vcf build-windows-speed/benchmark.vcf.gz --threads 1 -o build-wsl-speed/summary-bgzf-t1
./build-wsl-speed/GWLD summarize --vcf build-windows-speed/benchmark.bcf --threads 1 -o build-wsl-speed/summary-bcf-t1
```

Use Windows executable/output paths for the Windows runs. No OS file-cache flush was performed. Compare all output bytes after changing CRLF to LF, retaining original output hashes in the raw CSVs.

## R encoding

Run from the repository root using the native R installation on each platform:

```bash
Rscript --vanilla GWLD-R/tests/benchmarks/codegeno.R --source --iterations=20
```

The measured Windows installation was `E:/software/R-4.6.1`; the WSL run used a locally extracted R under `.deps/r` with its library paths configured. Both need Rcpp, RcppArmadillo and a working native R compilation toolchain. `--source` compiles the current encoder directly, avoiding unrelated plotting-package dependencies; default invocation still uses the installed GWLD package. R encoding uses the native R compilation defaults, not the standalone C++ `-O3` build above.

This benchmark has no discarded warmup. Each version is timed in three alternating batches of 20 calls, with explicit garbage collection before each batch. The median batch time divided by 20 gives seconds per call; output and attributes are compared after each batch. The 20-call option was added after an initial single-call Windows trial returned zero elapsed time for a fast operation. That unusable trial is excluded from the archived results.
