# Exact constant-column counting — 2026-11-11

The baseline is **18c656c**. This is an incremental 2.0-to-2.0 comparison, separate from the 1.x and earlier packing measurements. It measures the complete integer-matrix engine, including preprocessing, caches and allocation, and excludes file I/O, input generation, R conversion/startup and plotting.

## Exact shortcut and eligibility

If locus A is fully observed and every sample has genotype category r, its joint table with locus B is known from B's observed marginals:

- Cell (r, k) equals B's count for category k; all other rows are zero.
- The paired sample count is B's observed total, even if B has missing samples.
- The transposed construction applies when B is the complete constant locus.

The shortcut supports categories 0, 1 and 2. An all-heterozygous column (category 1) is not an allele-fixed locus. All five statistics still use their existing evaluation formulas and LD solver. In particular, the ambiguous unphased LD sign in all-heterozygous cases retains the preceding solver's behavior.

A single-category locus containing missing samples cannot establish the paired overlap by itself, so it is not eligible. All-missing partners yield zero paired observations and the same missing results as before. No variants are removed or reordered.

Classification is prepared once from validated summaries. Matrix dispatch selects a template with constant-column checks enabled only when an eligible column exists; matrices without eligible columns use the ordinary counting template. Existing NEON/POPCNT/AVX2/AVX-512 dispatch and numerical tolerances are unchanged. The new assembly probe reconstructs the table with loads/stores rather than population counts; it still uses a stack buffer in the inspected standalone probe.

## Measurements

Native macOS 15.7.9 arm64, Homebrew Clang 22.1.8, C++17, `-O3 -DNDEBUG`, OpenMP enabled with **one thread**, no fast-math or global native/AVX flags. Before and after executables use the same compiler, Armadillo headers and OpenMP library.

Observed synthetic genotypes use uniform 0/1/2 values, seed 2026. The constant-column fraction is either 0% (control) or 50%. Constant positions are shuffled with a separate seed 2027 and assigned categories 0/1/2 in turn. Missingness of 0% or 25% applies to the other columns; selected constant columns are overwritten with fully observed values. Thus the 50% constant / 25% missing workload has approximately 12.5% missing genotypes overall. This is a controlled synthetic stress case, not an estimate of typical biological data composition or a paper-dataset timing.

Each executable discards one warmup and takes three timed engine evaluations per method. The driver alternates before/after process order across three invocations and reports the median of those invocation medians. Measurements run without concurrent builds or test jobs. Short measurements remain sensitive to host scheduling and CPU frequency.

**65,536 samples × 128 variants, 50% complete constant columns, 25% missing in the other columns:**

| Method | Before (ms) | After (ms) | Speedup |
| --- | ---: | ---: | ---: |
| MI | 9.023 | 4.521 | 2.00× |
| RMI | 9.590 | 4.910 | 1.95× |
| r2 | 8.896 | 4.199 | 2.12× |
| D | 8.873 | 4.222 | 2.10× |
| D' | 8.703 | 4.271 | 2.04× |

The complete [engine table](engine-arm64.tsv) includes 512 × 1024, 4096 × 128 and 65536 × 128, all five methods, both constant fractions and both missingness settings. With 0% deliberately constant columns, speedup ratios range from 0.97 to 1.03. Those controls include small slowdowns as well as small improvements; this is not a universal speedup claim. Ratios from different dated comparisons must not be multiplied.

All **60** method/workload matrices match the baseline byte for byte, including missing positions. Per-matrix SHA-256 values are in [checks.json](checks.json). Compiler information, executable/source hashes and assembly instruction records are in [metadata.json](metadata.json).

## Correctness and platform checks

- Six OpenMP and five serial CTest groups pass, including the independent 37-case LD oracle.
- New count tests check both orientations, all three constant genotype categories, missing partners, completely missing partners and disjoint missing masks. They cover scalar and packed size boundaries and compare against explicit per-sample counting rather than the optimized pair helper.
- Native GCC 16, macOS x86-64 execution through Rosetta and address/undefined-behavior checks pass.
- Seven R scripts pass in OpenMP and serial installations, with CLI/VCF integration. R CMD check reports zero errors, warnings and notes.
- The complete 59-expression example, paper-subset validation and native matrix/stream interruption checks pass. The main local R installation is restored to its OpenMP build and the local CLI is updated.

Native/x86 assembly checks still pass. This session does not execute AVX2/AVX-512 counters on supporting physical x86 hardware, and no new native x86 speed claim is made.

Four exploratory rewrites were discarded because their engine timings had no stable benefit: byte NEON accumulators, expanded Newton divisions, per-method template dispatch and a shorter total-count reduction. Only the exact data-dependent shortcut is retained. Exploratory and verification logs remain under `build/speed-2026-11-11/`.

## Reproduce

Compile the current benchmark source against each header snapshot. From the repository root, on the recorded macOS LLVM setup:

```bash
mkdir -p build/constant-before build/constant-after
git archive 18c656c GWLD-R/inst/include/gwld | tar -x -C build/constant-before
CXX=/opt/homebrew/opt/llvm/bin/clang++
LLVM_PREFIX=/opt/homebrew/opt/llvm
for revision in before after; do
  if [ "$revision" = before ]; then
    SHARED_INCLUDE=build/constant-before/GWLD-R/inst/include
  else
    SHARED_INCLUDE=GWLD-R/inst/include
  fi
  "$CXX" -std=c++17 -O3 -DNDEBUG -D_REENTRANT \
    -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_BLAS \
    -DARMA_DONT_USE_LAPACK -DARMA_DONT_USE_SUPERLU \
    -fopenmp=libomp -I /opt/homebrew/include -I "$SHARED_INCLUDE" \
    GWLD-C++/tests/kernel_benchmark.cpp \
    -L "$LLVM_PREFIX/lib" -Wl,-rpath,"$LLVM_PREFIX/lib" \
    -o "build/constant-$revision/benchmark"
done
python3 tests/benchmarks/2026-11-11/compare.py \
  --before build/constant-before/benchmark \
  --after build/constant-after/benchmark \
  --output build/constant-comparison
```

Use native compiler/library paths on other systems and identical settings in both builds. `std::shuffle` placement may differ across standard-library implementations, so compare within one toolchain rather than treating these data as a cross-toolchain input fixture.

The benchmark now prints all five methods and accepts a final optional `constant_percent` (default 0). A single workload can be run with `benchmark 65536 128 25 1 - uniform 50`. The driver writes timing medians, matrix verification hashes and the last pair of exported matrix sets; any matrix mismatch fails the comparison.
