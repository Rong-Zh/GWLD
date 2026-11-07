# Performance measurements

## Direct comparison with the original implementation

The baseline is commit `c8079ab`, the pre-rewrite 1.x code retained on `master`. It is not the first commit in the repository's history and is not an intermediate 2.0 implementation. The current production kernels measured here are from `c77d8f4`. The expanded benchmark harness and these measurements accompany this documentation update.

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

Native ARM64 NEON, scalar fallbacks, macOS x86_64 POPCNT and x86 AVX2 under QEMU have correctness checks. AVX-512 functions compile with GCC and Clang but have not yet run on supporting hardware. QEMU timings are not used for performance claims. These measurements therefore do not quantify Intel/AMD SIMD acceleration.

For a reproducible explanation of the LD differences, see [independent optimizer verification](../tests/ld/README.md). It includes raw genotype export, a base-R likelihood oracle and recorded per-case evidence.

For numerical and workflow checks, use CTest, R CMD check, [the paper validation script](../tests/paper/README.md) and the optional R interruption test. The latest recorded review passed both OpenMP/serial C++ configurations, serial R regressions, the full R example, native interruption checks and address/undefined-behavior checks. No tests can establish the absence of every possible defect; the benchmark reports the comparisons it actually performs.
