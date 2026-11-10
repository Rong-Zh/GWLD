# GWLD 2.0

GWLD performs genome-wide linkage disequilibrium analysis using mutual information (MI), reduced mutual information (RMI), r², D and D′. It provides an R package for analysis and visualization and a C++ command-line program for batch processing. Both interfaces use the same numerical kernels and HTSlib-based VCF reader.

This README describes the `feature/2.0-rewrite` branch. The command-line interface has changed from 1.x; see [the CLI reference](GWLD-C++/README.md) before migrating existing commands.

## What changed in 2.0

| Area | Changes |
| --- | --- |
| Numerical engine | Shared R/C++ implementation; packed genotype counting; bounded MI/RMI caches; reusable locus statistics; OpenMP parallelism. |
| CPU support | ARM64 NEON and runtime-dispatched x86 AVX-512 VPOPCNTDQ, AVX2, POPCNT and portable SWAR counting. No global AVX requirement. |
| LD estimation | A faster safeguarded solver for certified concave likelihoods; comparison of stationary points and boundaries for potentially multimodal likelihoods. This corrects cases where the old optimizer chose the wrong likelihood mode. |
| CLI | Three verb-based subcommands: `calculate`, `convert` and `summarize`. Lowercase, hyphenated options; legacy options have been removed. |
| Input | Shared streaming reader for VCF, gzip/BGZF VCF and BCF; direct integer genotype matrices in R; preserved sample identities and numeric-looking IDs. |
| Pair analysis | Chromosome-specific distance indexing, unsorted-coordinate support, bounded output batches and optional R file streaming. A full LD matrix is optional. |
| R workflows | Reader results can be passed directly to `HeatMap`; existing matrices can be reused. Numeric genotype matrices are not recoded. Invalid fractional genotypes are rejected before integer conversion. |
| Visualization | An independent `color_legend()` component, automatic label spacing and a compact circos layout with the legend at the upper right. ComplexHeatmap is not required. |
| Encoding and utilities | Cached parsing of distinct genotype strings, reusable hash storage, faster field parsing, genotype summaries and buffered legacy output. |
| Validation | Self-contained R package tests, SIMD and numerical tests, serial/parallel checks, interruption checks and a runnable adaptation of the 2023 paper example. |

See [CHANGELOG.md](CHANGELOG.md) for compatibility details and [performance measurements](docs/performance.md) for comparisons against the pre-rewrite implementation.

## Performance compared with 1.x

The baseline is commit `c8079ab`, the original implementation on `master` immediately before the 2.0 rewrite. These are **direct measurements against 1.x**, rather than products of speedups measured between intermediate 2.0 versions.

The table below measures full upper-triangular matrix calculations on macOS arm64 with LLVM 22.1.8, Release `-O3`, one thread, 512 samples, 384 variants and 25% missing genotypes. Both versions receive the same generated matrix. Each path is warmed up, then timed three times in alternating order; the table reports medians. Timing includes the calculation engine's preparation and matrix allocation, but excludes input generation, file I/O and plotting.

| Method | Original 1.x | Current 2.0 | Speedup |
| --- | ---: | ---: | ---: |
| MI | 0.567922 s | 0.002099 s | 270.6× |
| RMI | 0.589542 s | 0.004193 s | 140.6× |
| r² | 0.205812 s | 0.005215 s | 39.5× |
| D | 0.205867 s | 0.005103 s | 40.3× |
| D′ | 0.205578 s | 0.006204 s | 33.1× |

MI/RMI remain numerically equivalent to the original formulas within floating-point rounding. LD values can differ because 2.0 corrects and tightens the likelihood optimization; the benchmark reports those differences separately. Speedups depend on hardware, sample size, allele distribution, missingness and thread count. They describe the calculation stage, not the wall-clock speed of an entire analysis pipeline.

Additional datasets, numerical differences, raw results and reproduction commands are provided in [docs/performance.md](docs/performance.md). The [Windows/WSL comparison](docs/performance.md#windows-and-wsl-on-the-same-x86-64-host), measured on 2026-11-09 on the same i5-12400 host, adds native x86-64 results, thread scaling, VCF/gzip/BGZF/BCF input timings and R encoder measurements. Both platforms passed native POPCNT/AVX2 correctness checks. AVX-512 functions have been compiled with GCC and Clang, but have not yet been executed on supporting hardware.

## Build the C++ program

Requirements:

- CMake 3.14 or newer and a C++17 compiler.
- Armadillo headers; the CLI does not link the Armadillo wrapper, BLAS or LAPACK.
- OpenMP for parallel calculations, or a serial build with `GWLD_OPENMP=OFF`.
- HTSlib 1.15 or newer and its compression dependencies, discovered through `pkg-config` when using a native installation.
- The pinned [ArgumentsParser](https://github.com/Rong-Zh/ArgumentsParser) Git submodule.

From a checkout of this branch:

```bash
git submodule update --init --recursive
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/GWLD --version
```

The ArgumentsParser submodule currently uses an SSH URL. If you do not have GitHub SSH access, select its public HTTPS URL locally before initialization:

```bash
git config submodule.GWLD-C++/deps/ArgumentsParser.url https://github.com/Rong-Zh/ArgumentsParser.git
git submodule update --init --recursive
```

Linux x86-64 builds default to the bundled HTSlib static archive, which requires development packages for zlib, bzip2, liblzma and libdeflate. Other platforms use a native installation. To select native HTSlib explicitly, pass `-DGWLD_USE_BUNDLED_HTSLIB=OFF`. See [platform and dependency details](GWLD-C++/README.md#build-options).

### macOS with Homebrew

```bash
brew install cmake pkg-config armadillo htslib llvm libomp
git submodule update --init --recursive
cmake -S GWLD-C++ -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DGWLD_USE_BUNDLED_HTSLIB=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/install"
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build
./build/install/bin/GWLD --version
```

Use a fresh build directory when changing compilers. Set `-DGWLD_ARMADILLO_INCLUDE_DIR=/path/to/include` for a nonstandard header location. A serial build uses `-DGWLD_OPENMP=OFF` and must run with one calculation thread. macOS builds define `_REENTRANT` to expose the system's reentrant `lgamma_r` function.

## Command-line quick start

Run the following examples from the repository root. Exactly one input option is required: `--bed-prefix`, `--ped-prefix` or `--vcf`.

```bash
# Full MI matrix; the default method is RMI.
./build/GWLD calculate --bed-prefix example/example --method mi --threads 4 -o result

# Same-chromosome pairs within an inclusive 300 kb window; no full matrix.
./build/GWLD calculate --vcf example/example-3.vcf \
  --decay --max-distance-kb 300 --threads 4 -o decay

# Reuse a full matrix for decay and circos outputs.
./build/GWLD calculate --ped-prefix example/example-2 \
  --method rmi --heatmap --decay --circos -o analysis

# Standard VCF, integer genotype table, and allele-frequency summary.
./build/GWLD convert --bed-prefix example/example --format vcf -o converted
./build/GWLD convert --vcf example/example-3.vcf --format genotype -o genotypes
./build/GWLD summarize --vcf example/example-3.vcf -o frequencies

./build/GWLD --help
./build/GWLD calculate --help
```

| Output | Contents |
| --- | --- |
| `.heatmap` | Full matrix storage with values in the upper triangle; diagonal and lower triangle are `nan`. |
| `.decay` | Within-chromosome pairs: variant IDs, absolute distance in bp and value; no header. |
| `.circos` | Cross-chromosome pairs exceeding the threshold: two chromosome/position/ID triples and value; no header. |
| `.genotype.tsv` | Header, variant metadata and per-sample integer genotype codes `0/1/2/NA`; `0/1/2` count ALT alleles. This is not a VCF. |
| `.freq.tsv` | Header, variant metadata, REF/ALT frequencies and genotype missing rate. |
| `.vcf` | VCF with diploid GT fields. |

CLI options, defaults, allele orientation and migration notes are documented in [GWLD-C++/README.md](GWLD-C++/README.md). Full matrices require quadratic memory; use pair-only outputs when a matrix is unnecessary.

## Install the R package

The package requires R 4.0 or newer, a C++17 toolchain and the dependencies declared in [DESCRIPTION](GWLD-R/DESCRIPTION). It builds its shared VCF reader against Rhtslib, independently of the CLI's bundled Linux archive.

Install dependencies in R:

```r
install.packages(c("Rcpp", "RcppArmadillo", "viridisLite", "circlize", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Rhtslib")
```

Install the local checkout from a terminal:

```bash
R CMD INSTALL GWLD-R
```

Or install this branch from GitHub after installing `remotes`:

```r
remotes::install_github("Rong-Zh/GWLD", ref = "feature/2.0-rewrite", subdir = "GWLD-R")
```

OpenMP availability depends on the R toolchain. Selecting Homebrew LLVM for the CLI does not automatically configure R: R uses its own compiler settings and user Makevars. For a serial package build:

```bash
GWLD_OPENMP=OFF R CMD INSTALL GWLD-R
```

Use `cores = 1` with that build. Parallel requests are rejected explicitly when OpenMP is unavailable.

## R analysis workflow

Rows represent samples and columns represent variants. The calculation interfaces accept exact `0/1/2/NA` codes; `-1` also denotes missing input. Fractional values, infinities and other codes are rejected. Missing observations are excluded pairwise.

```r
library(GWLD)

input <- vcf2matrix("example/example-3.vcf", threads = 1)
g <- input$Genotype
info <- input$Info

mi <- MI(g)
rmi <- GWLD(g, method = "RMI", cores = 1)
r2 <- LD(g, method = "r2", cores = 1)

# Draw from a reader result, or reuse an existing matrix.
HeatMap(input, method = "r2")
HeatMap(input, LDmatrix = r2)

within <- calc_decay(g, info, method = "RMI", maxdist = 300)
between <- calc_circos(g, info, method = "RMI", threshold = 0.1)

# Reuse values rather than calculating them again.
within_cached <- decay(rmi, info, maxdist = 300)
between_cached <- circos(rmi, info, threshold = 0.1)

# Stream large pair outputs to TSV files with headers.
calc_decay(g, info, method = "RMI", file = "decay.tsv")
calc_circos(g, info, method = "RMI", threshold = 0.1, file = "circos.tsv")
```

`bed2matrix()` and `ped2matrix()` also return `list(Info, Genotype)` with integer genotypes. The manual VCF transposition and numeric conversion steps used in older examples are unnecessary. `ped2matrix(prefix, genotype = "char")` requests allele strings explicitly; `read.vcf()` retains the character-table interface. `codegeno()` converts diploid allele strings, supports `/` and `|` separators, and preserves already encoded numeric matrices.

C++ uses `int64_t` for coordinates; R reader POS columns remain numeric. Negative, fractional or out-of-range MAP/BIM coordinates are rejected.

VCF genotype codes count ALT alleles. BED retains the original allele-1 counting convention; PED uses the major allele as genotype code zero, with lexical tie breaking. Align allele orientation before comparing signed D between formats. PLINK sample names use escaped `FID:IID`, preserve leading zeros and reject duplicate family/individual pairs. VCF sample names are retained.

Named cached heatmaps must match the supplied variant IDs and order. Matrix and pair computations check for R interrupts between batches. An interrupted streamed file retains completed rows. R pair outputs use double-precision text and include headers when written with `file=`; CLI pair files retain their headerless format.

The complete example remains in [example/example.R](example/example.R). After installation, run it from the repository root:

```r
source("example/example.R", chdir = TRUE)
```

### Circos plots and legends

`circos.ideogram(info)` accepts chromosome/position metadata, and `circos.linksnp(between)` accepts the pair table directly. The independent `color_legend()` grid component supports vertical and horizontal color bars, titles, physical sizing, exact tick positions and label spacing. ComplexHeatmap is not a dependency.

The default circos legend sits at the upper right. `legend.position = c(x, y)` sets its center in normalized device coordinates; `legend.space` on `circos.ideogram()` reserves additional space to the right. Square graphics devices give a compact circular plot; rectangular devices preserve the circle's aspect ratio and can leave extra space. See `?color_legend`, `?circos.linksnp` and `?circos.ideogram` for controls.

## Tests and paper validation

```bash
# C++ numerical, SIMD, packed-count, utility and CLI checks.
ctest --test-dir build --output-on-failure

# R package installation, namespace, documentation, examples and six test scripts.
R CMD build GWLD-R
R CMD check --no-manual GWLD_2.0.0.tar.gz

# Optional native interruption checks after installing the package.
python3 GWLD-R/tests/interrupt_test.py
```

R test fixtures are included in the source package, so package checks do not depend on the repository's working directory. Individual scripts in `GWLD-R/tests` can also be run from the repository root. Set `GWLD_TEST_OPENMP=OFF` when testing a serial build. Cross-interface VCF checks can additionally use `GWLD_TEST_CLI` and `GWLD_TEST_VCF_CONVERTER`; see [the test instructions](GWLD-C++/README.md#tests).

An additional [independent LD accuracy check](tests/ld/README.md) reconstructs genotype counts in base R, compares every stationary/boundary candidate and checks a separate grid-refined search. It includes the maximum-difference pairs, a wrong-mode counterexample and tied solutions.

The latest macOS verification on 2026-11-07 passed R CMD check with **zero errors, warnings or notes**, all five CTest groups in both OpenMP and serial builds, all six R scripts in a serial build, the full example, native interruption checks and address/undefined-behavior checks. An additional stress run covered 100,000 likelihood cases.

For a runnable adaptation of Box 1 and numerical cross-checks from the [2023 GWLD paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10468308/):

```r
install.packages(c("genetics", "infotheo"))
```

```bash
Rscript tests/paper/paper_2023.R build/paper-validation 2
```

It checks 1,770 pairs from the packaged duck subset against independent calculations. This validates the subset and example workflow; it does not reproduce all species, full datasets or performance results in the paper. Sources, thresholds and output descriptions are in [tests/paper/README.md](tests/paper/README.md).

The [2026-11-10 assembly-assisted review](tests/assembly/README.md) records emitted SIMD instructions, corrected POPCNT dispatch requirements, coordinate validation and reproduction commands.

## Repository layout

| Path | Purpose |
| --- | --- |
| `GWLD-R/inst/include/gwld/` | Shared numerical, packed-count and VCF-reader implementation. |
| `GWLD-R/src/`, `GWLD-R/R/` | Native R adapters, genotype encoding, R workflows and plotting. |
| `GWLD-C++/src/`, `GWLD-C++/include/` | CLI orchestration and PLINK input/output. |
| `GWLD-C++/tests/`, `GWLD-R/tests/` | Regression tests and optional benchmarks. |
| `tests/paper/` | Paper example adaptation and archived source material. |
| `tests/benchmarks/` | Recorded benchmark measurements. |
| `example/` | Example data and the complete R demonstration. |

Historical implementations under `tests/legacy` are benchmark references; they are not production calculation paths. The optimized compatibility utilities and their frozen references are distinguished in [the legacy helper notes](GWLD-C++/tests/legacy/README.md).

## Citation and license

If you use GWLD in research, cite *GWLD: an R package for genome-wide linkage disequilibrium analysis*, [DOI: 10.1093/g3journal/jkad154](https://doi.org/10.1093/g3journal/jkad154). The R package declares GPL (>= 2); see [its license information](GWLD-R/DESCRIPTION) and the licenses supplied with bundled dependencies.
