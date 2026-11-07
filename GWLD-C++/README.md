# GWLD C++ 2.0

The command-line program shares its numerical engine and HTSlib-based VCF reader with the R package. Shared headers live in [`GWLD-R/inst/include/gwld`](../GWLD-R/inst/include/gwld); CLI orchestration and PLINK input/output live in `src` and `include`.

Build and installation examples are in [the project README](../README.md). All examples below assume the repository root is the working directory.

## Commands

```text
GWLD <command> [options]
```

| Command | Purpose |
| --- | --- |
| `calculate` | Calculate LD matrices, within-chromosome decay pairs and cross-chromosome circos pairs. |
| `convert` | Export standard VCF or ALT dosage TSV without calculating LD. |
| `summarize` | Export REF/ALT frequencies and genotype missing rates without calculating LD. |

Use `GWLD --help`, `GWLD calculate --help`, or `GWLD help calculate` for help. `GWLD -V` and `GWLD --version` print the version.

### Common options

Exactly one of the following input options is required:

| Option | Input |
| --- | --- |
| `-b, --bed-prefix PREFIX` | PLINK SNP-major `.bed`, `.bim` and `.fam` files. |
| `-p, --ped-prefix PREFIX` | PLINK `.ped` and `.map` files. |
| `-v, --vcf PATH` | Biallelic, diploid VCF, gzip/BGZF VCF or BCF. GT may appear anywhere in FORMAT. |

All commands accept:

| Option | Default | Meaning |
| --- | --- | --- |
| `-t, --threads N` | `1` | Positive calculation/decompression thread budget. A serial build accepts only `1`. |
| `-o, --output-prefix PREFIX` | `gwld` | Nonempty output prefix. |
| `-h, --help` | — | Show command help. |

### Calculate options

| Option | Default | Meaning |
| --- | --- | --- |
| `-m, --method NAME` | `rmi` | `rmi`, `mi`, `r2`, `d` or `dprime`; CLI names are lowercase. |
| `-d, --decay` | Off | Write within-chromosome pairs. |
| `--max-distance-kb KB` | `300` | Inclusive decay distance in kb; zero is allowed. An explicit value requires `--decay`. |
| `-c, --circos` | Off | Write cross-chromosome pairs whose value is strictly greater than the threshold. |
| `--circos-threshold VALUE` | `0.1` | Finite, exclusive circos threshold. An explicit value requires `--circos`. |
| `--by-chromosome` | Off | Write a separate matrix for each chromosome. |
| `--heatmap` | Automatic | Request a full matrix in addition to other selected outputs. |

With no pair or chromosome output selected, `calculate` writes a full matrix. Selecting `--decay`, `--circos` or `--by-chromosome` disables that default; add `--heatmap` to request it explicitly. When a full matrix is requested, pair outputs reuse its values.

```bash
./build/GWLD calculate --vcf example/example-3.vcf --method r2 -o ld
./build/GWLD calculate --bed-prefix example/example --decay --circos -t 4 -o pairs
./build/GWLD calculate --ped-prefix example/example-2 --by-chromosome -o chromosomes
./build/GWLD calculate --bed-prefix example/example --decay --heatmap -o combined
```

### Conversion and summary

`convert` accepts `-f, --format vcf|dosage`; the default is `vcf`.

```bash
./build/GWLD convert --ped-prefix example/example-2 --format vcf -o converted
./build/GWLD convert --vcf example/example-3.vcf --format dosage -o genotypes
./build/GWLD summarize --bed-prefix example/example -o frequency
```

## Output schemas

Outputs are tab-separated. Pair and matrix files retain the headerless CLI format; dosage and frequency tables have headers. Numerical output uses enough significant digits to round-trip a double.

| File | Layout |
| --- | --- |
| `PREFIX.heatmap` | Matrix in input variant order. Upper triangle contains results; diagonal and lower triangle are `nan`. |
| `PREFIX.chrNAME.heatmap` | Per-chromosome upper-triangular matrix in the original within-chromosome variant order. |
| `PREFIX.decay` | `ID_1 ID_2 Dist Value`; distance is absolute base-pair distance. |
| `PREFIX.circos` | `CHR_1 POS_1 ID_1 CHR_2 POS_2 ID_2 Value`. |
| `PREFIX.dosage.tsv` | `CHROM POS ID REF ALT`, followed by one column per sample containing ALT dosage `0/1/2`, or `NA`. |
| `PREFIX.freq.tsv` | `CHROM POS ID REF ALT REF_FREQ ALT_FREQ MISSING_RATE`. |
| `PREFIX.vcf` | VCF with diploid GT fields. |

An empty pair selection produces an empty file. Pair output order follows the input variant order and does not change with thread count. Frequency denominators include only called genotypes; all-missing variants have `nan` allele frequencies.

## Genotype and statistical conventions

- Missing genotypes are excluded pairwise. A pair with no complete observations returns a missing result.
- Monomorphic loci yield D = 0 and undefined r²/D′. MI/RMI follow their contingency-table definitions.
- RMI retains the original whole-locus category-count convention; LD retains the original whole-locus marginal-frequency convention, including when the complete-pair subset is smaller.
- VCF preserves REF/ALT dosage orientation. BED retains allele-1 dosage: BIM column 6 is REF and column 5 is ALT. PED uses the most frequent allele as REF, with lexical tie breaking.
- Signed D depends on allele orientation. Align alleles before comparing values from different formats.
- PLINK sample names are `FID:IID`. Within each component, `%` becomes `%25` and `:` becomes `%3A`. Numeric-looking identifiers retain leading zeros. Duplicate family/individual pairs are rejected. VCF sample names are retained.
- The LD optimizer compares feasible stationary points and boundaries for potentially multimodal likelihoods. Certified concave cases use safeguarded Newton iterations with a bracketed fallback. Values can therefore differ from the 1.x optimizer, especially in heterozygote-heavy cases.

Decay analysis indexes variants by chromosome and coordinate, so input need not be sorted. Its window is inclusive and uses absolute distance. Large coordinates are checked before arithmetic. BCF's own coordinate representation can impose additional limits.

A full matrix stores O(variants²) doubles. Pair-only requests avoid allocating that matrix, but still load the input genotype matrix and may produce a quadratic number of pairs. Pair results are written in bounded batches; this does not make the whole input a constant-memory stream.

Malformed input, unsupported method names and invalid thread counts cause a nonzero exit status. Genotype calculations reject multiallelic variants and nonmissing haploid GT values; missing diploid calls and a single missing `.` are supported.

## Migration from 1.x

The old flat invocation and legacy flag aliases have been removed.

| Old option | 2.0 replacement |
| --- | --- |
| `--bfile` | `--bed-prefix` |
| `--file` | `--ped-prefix` |
| `--2vcf` | `convert --format vcf` |
| `--code012` | `convert --format dosage` |
| `--allele-freq` | `summarize` |
| `--by-chrom` | `calculate --by-chromosome` |
| `--max-dist` | `calculate --decay --max-distance-kb` |
| `--threshold` | `calculate --circos --circos-threshold` |
| `--out-prefix` | `--output-prefix` |

Thread selection now uses `-t`; `-p` selects PED input. Update method names to `rmi`, `mi`, `r2`, `d` or `dprime`. The R interface retains its function names and accepts the documented aliases.

## Build options

| CMake option | Default | Purpose |
| --- | --- | --- |
| `GWLD_OPENMP` | `ON` | Enable OpenMP; set `OFF` for a serial build. |
| `GWLD_USE_BUNDLED_HTSLIB` | `ON` on Linux x86-64; otherwise `OFF` | Select the bundled Linux archive or a native HTSlib installation. |
| `GWLD_ARMADILLO_INCLUDE_DIR` | Discovered | Directory containing Armadillo headers. |
| `GWLD_BUILD_BENCHMARKS` | `OFF` | Build the optional performance programs. |
| `BUILD_TESTING` | `ON` | Build and register regression tests. |
| `CMAKE_INSTALL_PREFIX` | CMake default | Destination for the installed `bin/GWLD`. |

Optional benchmark targets require OpenMP for their historical reference even when the production CLI is configured as serial.

C++17 is required by the configured targets. The CLI uses Armadillo matrix storage without BLAS, LAPACK, SuperLU or the Armadillo wrapper library. `pkg-config` is required. Linux x86-64 bundled HTSlib builds require zlib, bzip2, liblzma, libdeflate and thread support. The supplied static archive is HTSlib 1.24 and cannot be used on macOS or another architecture. Native HTSlib must be at least 1.15.

`deps/ArgumentsParser` is a pinned Git submodule from [Rong-Zh/ArgumentsParser](https://github.com/Rong-Zh/ArgumentsParser). Initialize it before configuring CMake. The R package builds the shared reader against Rhtslib and does not use this CLI archive.

## SIMD dispatch

The shared engine keeps multiple counting implementations and selects an available one at runtime:

| Platform | Paths |
| --- | --- |
| ARM64 GCC/Clang | NEON and portable counting. |
| x86 GCC/Clang | AVX512F + VPOPCNTDQ for at least eight 64-bit words; otherwise AVX2 for at least sixteen words; otherwise POPCNT when supported, then SWAR. |
| MSVC x64 | POPCNT when supported, with SWAR fallback. |

No global `-mavx2`, `-mavx512*` or `-march=native` is required. AVX2 byte accumulators are reduced every 31 vectors; NEON accumulators every 4095 vectors. Scalar tails handle partial vectors. Complete columns use four intersections to reconstruct the 3×3 table from known marginals.

Packed counting is enabled for at least 128 samples and four variants; smaller inputs use the scalar path. Packed storage uses three bits per sample/variant, rounded to groups of 64 samples. MI count-log and RMI factorial caches are bounded at 8 MiB each. A separate half-integer log cache, also bounded at 8 MiB, is built only when missing-data workloads can amortize it. These tables use the original mathematical functions, without fast-math approximations.

ARM64 NEON, macOS x86_64 POPCNT and x86 AVX2 under QEMU have been checked against portable counting. AVX-512 has been compiled with GCC and Clang but has not been run on supporting hardware. QEMU checks establish correctness for the exercised paths, not native x86 performance.

## Tests

```bash
ctest --test-dir build --output-on-failure
./build/gwld_simd_test
```

CTest groups cover:

- `core`: independent numerical references, multimodal likelihood cases and distance windows.
- `packed`: scalar/packed equivalence, missing observations, padding and parallel results.
- `simd`: each available forced counting path, vector boundaries and accumulator reduction boundaries.
- `utils`: optimized compatibility helpers versus frozen results and file bytes, plus corrected boundaries.
- `cli`: actual commands, conversions, schemas, input validation and thread consistency; registered when Python is available.

Optional R/CLI integration checks, after installing GWLD:

```bash
GWLD_TEST_CLI="$PWD/build/GWLD" \
GWLD_TEST_VCF_CONVERTER="$PWD/build/gwld_vcf_fixture" \
Rscript GWLD-R/tests/vcf.R
```

The fixture converter creates BGZF VCF and BCF inputs for the same R-reader assertions. A separate serial build can be checked with `-DGWLD_OPENMP=OFF` in a fresh build directory.

## Benchmarks

```bash
cmake -S GWLD-C++ -B build -DGWLD_BUILD_BENCHMARKS=ON
cmake --build build -j

# Original 1.x versus current 2.0: samples, variants, threads, missing percent.
./build/gwld_benchmark 512 384 1 25

# Current scalar versus packed engine: same argument order.
./build/gwld_engine_benchmark 512 384 1 25

# Current engine timing: samples, variants, missing percent, threads,
# optional matrix-output prefix, optional uniform/hwe distribution.
./build/gwld_kernel_benchmark 512 1536 5 1 - hwe

# Compatibility utilities and VCF input parsing.
./build/gwld_utils_benchmark build/utils-perf
./build/gwld_vcf_benchmark example/example-3.vcf
```

The scalar engine benchmark uses the current solver and formulas, so it isolates counting/caching effects; it is not a frozen 1.x or historical 2.0 binary. The original-version benchmark instead links `tests/legacy/methods.cpp` and `utils_reference.cpp`, preserving the algorithms from `c8079ab`.

See [measured results and methodology](../docs/performance.md), [recorded raw results](../tests/benchmarks/2026-11-07-v1-v2.tsv) and [legacy reference notes](tests/legacy/README.md). Historical 1.x manuals remain in the repository; use this reference and `--help` for 2.0 behavior.
