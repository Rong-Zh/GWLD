# Performance measurements

Results are retained separately by baseline, platform and measurement date. Historical timings describe the recorded revisions; they are not automatically updated when code changes:

- [Latest incremental macOS packing optimization, 2026-11-10](#incremental-packing-optimization--2026-11-10).
- [Windows and WSL comparison, 2026-11-09](#windows-and-wsl-on-the-same-x86-64-host).
- [Original macOS C++ measurements, 2026-11-07](#direct-comparison-with-the-original-implementation).
- [Original macOS R encoding measurements, 2026-11-07](#r-genotype-encoding-versus-1x).

## Incremental packing optimization — 2026-11-10

This comparison uses **68dade2 as the before version** and the subsequent register-based packing implementation as the after version. It measures an additional 2.0-to-2.0 improvement. The original-version and Windows/WSL tables below retain their recorded source revisions and timings; speedups from different tables must not be multiplied.

The new implementation builds three genotype bit planes in registers, writes each word once and reuses allocated storage without clearing words that will be overwritten. ARM64 uses 16-sample NEON comparisons; GCC/Clang x86-64 uses baseline SSE2 comparisons and movemasks. Existing POPCNT, AVX2 and AVX-512 pair counters and their runtime selection are unchanged. Scalar tails preserve the exact sample positions, including missing values. No statistical formula, optimizer tolerance, approximate logarithm or fast-math flag changes are involved.

Native macOS 15.7.9 arm64, Homebrew Clang 22.1.8, C++17, Release `-O3 -DNDEBUG`, OpenMP enabled with **one requested thread**, seed 2026. Genotypes are synthetic uniform 0/1/2 with independently sampled missingness. Complete engine timings include preparation, summaries, packing, caches, matrix allocation and all upper-triangular calculations; input generation, file I/O, R startup and plotting are excluded.

**65,536 samples × 128 variants, 25% missing:**

| Measurement | Before (ms) | After (ms) | Speedup |
| --- | ---: | ---: | ---: |
| Packing only, reused storage | 14.390 | 2.073 | 6.94× |
| Complete MI matrix engine | 21.858 | 9.627 | 2.27× |
| Complete RMI matrix engine | 22.417 | 10.259 | 2.19× |
| Complete r² matrix engine | 21.912 | 9.730 | 2.25× |

**Complete engine speedups at other sizes, 25% missing:**

| Samples × variants | MI | RMI | r² |
| --- | ---: | ---: | ---: |
| 4,096 × 128 | 2.03× | 1.70× | 1.69× |
| 512 × 768 | 1.09× | 1.04× | 1.02× |

Packing-only benchmarks warm up once and time 15 reused builds, then report the median of five alternating before/after invocations. Engine benchmarks warm up once and time three fresh-engine evaluations per method, then report the median of three alternating before/after invocations. Reuse differs deliberately between these two measurements: packing throughput is not a standalone estimate of fresh-engine speed. Millisecond and sub-millisecond measurements are sensitive to host load and CPU frequency.

MI, RMI and r² matrices match the before version byte for byte in all six engine workloads (three dimensions × two missingness settings). D/D′ retain the same counting implementation and numerical solver; independent LD and regression checks pass, but no new D/D′ timings are claimed here. Output equality against this immediate baseline is distinct from the known LD differences against 1.x.

Native Clang/GCC SIMD checks, address/undefined-behavior checks, six OpenMP CTest groups, five serial CTest groups, six R scripts with CLI/VCF integration and R CMD check pass. The earlier full-example, paper and serial-R workflow evidence belongs to the preceding review; those workflows were not rerun for this packing-only follow-up. The reproducible [dated report](../tests/benchmarks/2026-11-10/README.md) includes measured tables, assembly instruction counts, executable hashes and exact reproduction commands.

The x86-64 packing tests also pass under Rosetta. Those supplementary translation timings are not native x86 measurements: mixed data is faster, while all-missing input takes about 33% longer (0.278 to 0.369 ms for 4,096 × 128). AVX2/AVX-512 assembly checks pass, but this follow-up did not execute those counters on supporting physical x86 hardware. Neither a universal input-independent speedup nor a new native Windows/WSL speedup follows from this measurement.

## Windows and WSL on the same x86-64 host

Measured on **2026-11-09**, using production source at `bca4f90` and the same Intel Core i5-12400 host (12 logical processors). Windows runs are native PE x86-64 executables; WSL runs are Linux executables under WSL2. The original-version baseline remains `c8079ab`, compiled into the same executable as the current engine on each platform.

| Environment | Windows native | WSL2 |
| --- | --- | --- |
| OS | Windows 10 Pro 22H2, build 19045.6466 | Ubuntu 24.04; Linux 6.18.58.1-microsoft-standard-WSL2 |
| C++ compiler | Rtools GCC 14.3.0 | Ubuntu GCC 13.3.0 |
| Flags | C++17, `-O3 -DNDEBUG -fopenmp` | C++17, `-O3 -DNDEBUG -fopenmp` |
| Armadillo headers for C++ tests | 15.0.3 | 12.6.7 |
| HTSlib | Windows static library from Rhtslib 3.8.0 | Bundled Linux x86-64 HTSlib 1.24 |
| R / Rcpp / RcppArmadillo | 4.6.1 UCRT / 1.1.2 / 15.6.2.1 | 4.3.3 / 1.0.12 / 0.12.8.1.0 |
| Available x86 paths tested | Portable, automatic dispatch, POPCNT, AVX2 | Portable, automatic dispatch, POPCNT, AVX2 |

Both builds disable Armadillo's BLAS/LAPACK/wrapper, leave runtime CPU dispatch enabled and use neither `-ffast-math` nor global `-march=native`. `OMP_DYNAMIC=FALSE`; requested thread counts are identical. Benchmarks run sequentially, with no other benchmark jobs in parallel. These measurements compare the two installed toolchains, libraries and operating environments together: they do **not** isolate the operating system as the only variable.

### Direct comparison with 1.x

The input generator, seed 2026, workload dimensions and timing protocol are the same on both platforms. Each implementation has one discarded warmup followed by three alternating timed runs; tables report medians. Timings include preprocessing, bit packing, caches, matrix allocation and upper-triangular calculation, but exclude input generation, file I/O and result verification.

**512 samples × 384 variants, 25% missing, one thread** (73,536 pairs):

| Method | Windows 1.x (s) | Windows 2.0 (s) | Speedup | WSL 1.x (s) | WSL 2.0 (s) | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MI | 1.036945 | 0.003889 | 266.7× | 0.946671 | 0.004023 | 235.3× |
| RMI | 1.109446 | 0.014844 | 74.7× | 0.970964 | 0.008552 | 113.5× |
| r² | 0.336274 | 0.010687 | 31.5× | 0.291517 | 0.010978 | 26.6× |
| D | 0.334958 | 0.010404 | 32.2× | 0.286602 | 0.010275 | 27.9× |
| D′ | 0.343568 | 0.011655 | 29.5× | 0.283622 | 0.010474 | 27.1× |

**4096 samples × 256 variants, 1% missing, four threads** (32,640 pairs):

| Method | Windows 1.x (s) | Windows 2.0 (s) | Speedup | WSL 1.x (s) | WSL 2.0 (s) | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MI | 1.592590 | 0.002611 | 609.9× | 1.595028 | 0.002968 | 537.4× |
| RMI | 1.625216 | 0.004440 | 366.0× | 1.740673 | 0.003333 | 522.3× |
| r² | 0.184634 | 0.003144 | 58.7× | 0.166302 | 0.002897 | 57.4× |
| D | 0.185371 | 0.003381 | 54.8× | 0.181178 | 0.002948 | 61.5× |
| D′ | 0.173628 | 0.003103 | 56.0× | 0.173605 | 0.003040 | 57.1× |

MI/RMI maximum absolute differences from 1.x are below `1.25e-14`; all comparisons have zero missing-position mismatches. LD uses the corrected optimizer, so differences from 1.x are expected: across these runs, maximum absolute differences are `3.34e-5` for r², `3.88e-5` for D and `1.70e-4` for D′. These are not claims of bitwise equivalence to the original LD solver.

The separate current-solver scalar-versus-packed benchmark gives MI speedups of 42.15× on Windows and 35.30× on WSL, RMI speedups of 13.52× and 18.42×, and LD speedups of approximately 13.3–13.7×. LD results match exactly between those two paths; MI/RMI differences are below `6.46e-15`. This comparison isolates current execution paths rather than reconstructing a historical version.

### Current-version thread scaling

**512 samples × 2048 variants, 25% missing** (2,096,128 pairs), one discarded warmup followed by three timed runs. Input generation and binary matrix export are outside the timed region.

| Method | Threads | Windows (s) | WSL (s) |
| --- | ---: | ---: | ---: |
| MI | 1 | 0.092852 | 0.100071 |
| MI | 4 | 0.029873 | 0.033718 |
| MI | 12 | 0.017699 | 0.029103 |
| RMI | 1 | 0.403337 | 0.221236 |
| RMI | 4 | 0.124099 | 0.071806 |
| RMI | 12 | 0.094900 | 0.058847 |
| r² | 1 | 0.293537 | 0.304253 |
| r² | 4 | 0.084551 | 0.091440 |
| r² | 12 | 0.048979 | 0.060518 |

At 12 threads, the speedups over each platform's own single-thread run are 5.25× / 3.44× for MI, 4.25× / 3.76× for RMI and 5.99× / 5.03× for r² (Windows / WSL). Windows is faster for MI on this workload, WSL is faster for RMI, and single-thread r² times are close. No universal Windows-versus-WSL speed ranking follows from these results. Scheduling and timer noise matter especially for millisecond-scale runs.

### VCF parsing and complete summary commands

Both platforms read the **same files on the E: drive**, accessed through `/mnt/e` in WSL. The uncompressed input contains 512 samples and 12,000 variants, uses GT:DP fields and is 41,836,454 bytes. It repeats four genotype patterns, making the compressed files unusually small; this is not a representative compression ratio for biological data. Input SHA256 is recorded with the [reproduction notes](../tests/benchmarks/2026-11-09/README.md).

Reader-only comparison, three alternating runs, including opening, parsing, metadata and integer-matrix construction:

| Environment | Legacy stream (s) | HTSlib (s) | Speedup |
| --- | ---: | ---: | ---: |
| Windows | 0.621250 | 0.213558 | 2.91× |
| WSL | 1.201190 | 0.285397 | 4.21× |

Decoded genotypes and metadata agree with the stream reference on both platforms. Complete `GWLD summarize` commands include process startup, reading, frequency calculation and output writing; one warmup is discarded, then three executions are timed:

| Format | Input bytes | Windows, 1 thread (s) | Windows, 4 threads (s) | WSL, 1 thread (s) | WSL, 4 threads (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| VCF | 41,836,454 | 0.253181 | 0.259345 | 0.322848 | 0.324592 |
| Ordinary gzip VCF | 182,685 | 0.274629 | 0.280069 | 0.226802 | 0.226164 |
| BGZF VCF | 252,298 | 0.292762 | 0.272826 | 0.219841 | 0.221570 |
| BCF | 169,452 | 0.086174 | 0.080461 | 0.079201 | 0.077678 |

These are warm-cache measurements, not cold-disk throughput tests. Small compressed inputs show little benefit from extra decompression threads. WSL's mounted-drive access, differing HTSlib/compression libraries and process-startup costs all affect these complete-command timings.

### R encoding and correctness evidence

The R encoding comparison uses the same 2,000,000-genotype workload and seed as the macOS benchmark below. Both the current and original encoder are compiled directly with `Rcpp::sourceCpp`; no complete GWLD package installation is required. Each of the three alternating timed batches executes an operation 20 times, then elapsed time is divided by 20 and the median is reported. Batching avoids zero-duration readings on the Windows elapsed-time clock. Compilation and garbage collection before each batch are excluded; allocations and any garbage collection within a batch are included. Full output and attribute equality is checked after every batch.

| Operation | Windows 1.x (s/call) | Windows 2.0 (s/call) | Speedup | WSL 1.x (s/call) | WSL 2.0 (s/call) | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `recode` | 0.149500 | 0.014000 | 10.68× | 0.143400 | 0.013350 | 10.74× |
| `stats_freq` | 0.088000 | 0.007500 | 11.73× | 0.084800 | 0.007200 | 11.78× |
| `code_Mat` | 0.165000 | 0.015500 | 10.65× | 0.149650 | 0.013700 | 10.92× |

R uses each installation's normal package-compilation flags, rather than the standalone C++ benchmark flags. R, Rcpp and RcppArmadillo versions differ as recorded above. These are encoder measurements, not full R analysis or plotting timings. Full Windows R package installation was not validated in this run because `circlize` was unavailable; the direct-source benchmark does not need it.

Correctness evidence for these runs:

- Windows core, packed and SIMD executables passed. WSL CTest passed all five tests (`utils`, `simd`, `packed`, `core`, `cli`). Native POPCNT and AVX2 paths were executed on both platforms; AVX-512 was not available on this host.
- Every single-thread output cell was compared between the Windows and WSL 2048-variant matrices: 2,096,128 finite pairs per method, with matching missing positions. MI and r² are identical; RMI maximum absolute difference is `8.41e-15`. These checks compare values, not just checksums. Within each platform, thread-count checksums also agree.
- All 16 complete frequency outputs (two platforms × four formats × two thread counts) match after normalizing Windows CRLF and Linux LF line endings.
- The R encoding benchmark checked complete outputs and attributes after every timed batch.

Raw C++ timings, R timings, test logs and cross-platform result checks are archived in [the dated results directory](../tests/benchmarks/2026-11-09/README.md), together with exact commands and environment details. The existing macOS results below remain a separate historical measurement; they were not rerun on this x86-64 host.

## Direct comparison with the original implementation

The baseline is commit `c8079ab`, the pre-rewrite 1.x code retained on `master`. It is not the first commit in the repository's history and is not an intermediate 2.0 implementation. The production kernels in the following macOS measurements are from `c77d8f4`.

`GWLD-C++/tests/legacy/methods.cpp` retains the calculation algorithms from that baseline, with translated comments and references qualified to the frozen `legacy_reference` utilities. Those utilities are kept separately from the optimized compatibility helpers so subsequent utility changes cannot silently accelerate the baseline. The R encoding reference likewise retains `c8079ab`'s implementation with comment-only changes.

### Environment and timing protocol

Measured on 2026-11-07:

- macOS 15.7.9, arm64 Apple Silicon; current packed kernels use NEON.
- Homebrew LLVM/Clang 22.1.8; CMake Release `-O3 -DNDEBUG`, OpenMP available, one requested thread for both implementations.
- Both C++ implementations are compiled into the same benchmark executable with the same compiler and flags, using the installed Armadillo headers (15.4.0). This is a source-algorithm comparison in the current environment, not a comparison of historical release binaries or toolchains.
- Deterministic `std::mt19937` seed 2026. Missingness is sampled independently; other values are uniform over `0/1/2`. Percentages specify generation probabilities, not an exact count of missing cells. These are synthetic datasets, not the paper's complete biological datasets.
- One warmup per implementation and method, followed by three alternating timed executions. Each reported time is the median of those three executions. Runs are sequential, without other benchmark jobs running concurrently.
- The timed region includes locus summaries, bit packing, caches, matrix allocation and all upper-triangular pair calculations. Input generation, result checking, file I/O, plotting and R startup are excluded. Each version performs its own preprocessing; no reusable prepared engine is supplied to the new version.
- Speedup = median original time / median current time. Numbers are measurements from this machine, with no cross-stage multiplication and no extrapolation to other CPUs.

### Full results

| Samples × variants | Missing | Method | 1.x (s) | 2.0 (s) | Speedup | Max. absolute difference |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 512 × 384 | 25% | MI | 0.567922 | 0.002099 | 270.6× | 6.71e-15 |
| 512 × 384 | 25% | RMI | 0.589542 | 0.004193 | 140.6× | 6.57e-15 |
| 512 × 384 | 25% | r² | 0.205812 | 0.005215 | 39.5× | 3.33e-05 |
| 512 × 384 | 25% | D | 0.205867 | 0.005103 | 40.3× | 3.87e-05 |
| 512 × 384 | 25% | D′ | 0.205578 | 0.006204 | 33.1× | 0.00017 |
| 512 × 768 | 0% | MI | 1.930724 | 0.004466 | 432.3× | 7.11e-15 |
| 512 × 768 | 0% | RMI | 1.980312 | 0.005205 | 380.4× | 8.98e-15 |
| 512 × 768 | 0% | r² | 0.731702 | 0.017697 | 41.3× | 2.63e-05 |
| 512 × 768 | 0% | D | 0.730753 | 0.017416 | 42.0× | 3.67e-05 |
| 512 × 768 | 0% | D′ | 0.732281 | 0.022482 | 32.6× | 0.00016 |
| 4096 × 128 | 1% | MI | 0.376551 | 0.001189 | 316.6× | 9.94e-15 |
| 4096 × 128 | 1% | RMI | 0.376939 | 0.001487 | 253.6× | 1.11e-14 |
| 4096 × 128 | 1% | r² | 0.146288 | 0.001584 | 92.4× | 5.34e-06 |
| 4096 × 128 | 1% | D | 0.146111 | 0.001595 | 91.6× | 3.85e-05 |
| 4096 × 128 | 1% | D′ | 0.146224 | 0.001657 | 88.2× | 0.000156 |

Raw precision is preserved in [2026-11-07-v1-v2.tsv](../tests/benchmarks/2026-11-07-v1-v2.tsv). All 15 comparisons had zero missing/nonfinite-position mismatches. The benchmark checks every valid upper-triangular pair; the table reports the maximum absolute difference from the original implementation.

Across these three workloads, MI was 270.6–432.3× faster, RMI 140.6–380.4× faster, r² 39.5–92.4× faster, D 40.3–91.6× faster and D′ 32.6–88.2× faster. These ranges describe these inputs and this run; they are not guaranteed speedups for every dataset. Very short new-version runtimes are sensitive to scheduling and timer noise.

### Why values and speedups differ by method

The original implementation repeatedly copies genotype columns, constructs categories and counts observations for each pair. The new implementation prepares locus statistics and packed category bits once, then counts pairs with bitwise intersections and population counts. Complete columns need only four intersections to reconstruct the 3×3 table. MI/RMI also reuse bounded caches and complete-locus statistics. Those changes remove a large fraction of the original per-pair overhead.

LD additionally solves a likelihood problem for each pair. The new solver uses safeguarded Newton iterations when concavity can be certified and isolates stationary points in potentially multimodal cases. As a result, r²/D/D′ still perform work beyond contingency-table counting; it is expected that MI or RMI can run faster than r² on some inputs.

MI/RMI maximum differences are at floating-point rounding scale (below `1.2e-14` in these runs). For r²/D/D′, 2.0 deliberately changes the optimizer's convergence and solution selection. Nonzero LD differences must not be described as bitwise equivalence. Independent numerical tests, a multimodal counterexample, paper-subset validation and likelihood stress checks provide correctness evidence separately from the timing comparison.

### Reproduce the direct comparison

```bash
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release \
  -DGWLD_BUILD_BENCHMARKS=ON
cmake --build build --target gwld_benchmark -j

# Argument order: samples variants threads missing_percent.
./build/gwld_benchmark 512 384 1 25
./build/gwld_benchmark 512 768 1 0
./build/gwld_benchmark 4096 128 1 1
```

Use the macOS compiler/dependency options from the [installation instructions](../README.md#macos-with-homebrew) when configuring a fresh Apple Silicon build. The benchmark prints both times, speedup, maximum numerical difference and missing-position mismatch count. MI/RMI differences above `1e-10` or missing-position mismatches make it fail. LD differences are reported, rather than rejected merely for differing from the old optimizer; use the numerical regression suite to validate the current solution.

## R genotype encoding versus 1.x

The optional R benchmark uses 2,000,000 character genotypes in a 2,000-row matrix, sampled from `A/A`, `A/G`, `G/A`, `G/G` and `NA` with R seed 2026. It compares `recode`, `stats_freq` and the native matrix encoder against the frozen 1.x source loaded with `Rcpp::sourceCpp`.

The following measurements used R 4.6.1, RcppArmadillo 15.6.2.1 and the locally installed GWLD 2.0.0 on the same macOS arm64 host. Each case uses three alternating executions and median elapsed time, with garbage collection outside the timed call. Unlike the C++ matrix benchmark, this script does not include an explicit discarded warmup. All timed outputs are checked for equality, including attributes.

| Operation | Original 1.x (s) | Current 2.0 (s) | Speedup |
| --- | ---: | ---: | ---: |
| `recode` | 0.130 | 0.006 | 21.7× |
| `stats_freq` | 0.078 | 0.005 | 15.6× |
| `code_Mat` | 0.139 | 0.007 | 19.9× |

The new implementation parses each distinct interned genotype string once per column and reuses open-addressed lookup storage. It supports full allele tokens and missing alleles; this benchmark uses ordinary biallelic inputs where the old and new semantics agree. Raw values are in [2026-11-07-encoding.tsv](../tests/benchmarks/2026-11-07-encoding.tsv). The fastest calls are only a few milliseconds, so the ratios have limited precision.

```bash
Rscript GWLD-R/tests/benchmarks/codegeno.R
```

This script requires an installed GWLD package, Rcpp and a working C++/R toolchain. Compilation and package startup are outside its timed sections.

## Other benchmark tools

| Program | Comparison and scope |
| --- | --- |
| `gwld_benchmark` | Frozen pre-rewrite 1.x calculation source versus the current complete matrix engine; the primary comparison above. |
| `gwld_engine_benchmark` | Scalar versus packed counting/caching using the current numerical solver in both paths. This is an implementation-path comparison, not a frozen historical version. |
| `gwld_kernel_benchmark` | Current engine timing with uniform or variable-MAF HWE synthetic inputs; optional binary matrices support comparisons between builds. |
| `gwld_utils_benchmark` | Frozen compatibility utilities versus optimized helpers, including result and output-byte checks; also compares the production whitespace scanner against a stream-based scanner. |
| `gwld_vcf_benchmark` | Legacy stream parser versus HTSlib reader on the same supported VCF, checking decoded data and timing opening/reading/matrix construction. |

For reference, the separately measured compatibility utility improvements on the same host were approximately 3.24× for splitting, 3.08× for production field parsing, 3.17× for file reading, 3.23× for frequency calculation, 2.25× for allele counting and 4.21× for matrix writing. Those are function-level measurements; they cannot be multiplied into the matrix-engine speedups or interpreted as complete CLI speedups.

Intermediate 2.0-to-2.0 optimization measurements are intentionally not used as the original-version baseline. Re-running the scalar engine benchmark today uses today's solver and summary functions, so it does not reconstruct an earlier commit automatically.

## Limits and correctness checks

Full matrices remain O(variants²) in storage and pair count. Pair-only decay output avoids allocating a full matrix and restricts work to a chromosome/distance window, while input genotypes remain resident. The benefit depends on the number of eligible pairs; there is no universal end-to-end speedup implied by the kernel tables.

Native ARM64 NEON, scalar fallbacks, macOS x86_64 POPCNT, x86 AVX2 under QEMU, and native Windows/WSL x86-64 POPCNT and AVX2 have correctness checks. AVX-512 functions compile with GCC and Clang but have not yet run on supporting hardware. QEMU timings are not used for performance claims. The Windows/WSL timings above exercise native x86 dispatch, but do not isolate the speedup of one SIMD instruction set over another.

For a reproducible explanation of the LD differences, see [independent optimizer verification](../tests/ld/README.md). It includes raw genotype export, a base-R likelihood oracle and recorded per-case evidence.

For numerical and workflow checks, use CTest, R CMD check, [the paper validation script](../tests/paper/README.md) and the optional R interruption test. The preceding dispatch/coordinate review passed OpenMP/serial C++ configurations, serial R regressions, the full R example, native interruption checks and address/undefined-behavior checks. The packing follow-up has its own verification scope listed above. No tests can establish the absence of every possible defect; the benchmark reports the comparisons it actually performs.
