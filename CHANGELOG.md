# Changelog

## Unreleased — R interface review, 2026-11-11

- Preserve literal `NA` metadata in `read.vcf()`; restrict missing-value conversion to generated genotype placeholders.
- Make `savefig()` apply scale and physical/pixel units, support custom device functions, validate sizes and restore the previous device after drawing errors.
- Add VCF metadata and figure-export regressions; document verification in `docs/review-2026-11-11.md`.

## Unreleased — assembly review, 2026-11-10

- Build packed genotype planes in registers with 16-sample NEON/SSE2 blocks, replacing per-sample scatter updates. Preserve exact sample positions, missing/invalid-code handling and scalar tails; reuse packed allocations without redundant clearing.
- Record incremental packing and engine benchmarks against `68dade2`, with byte-identical MI/RMI/r² matrices and explicit x86/Rosetta limits in `tests/benchmarks/2026-11-10/README.md`.

- Match x86 SIMD dispatch to the instructions emitted by Clang: vector kernels now explicitly require scalar POPCNT and AVX2; the AVX-512 path additionally requires AVX512F/VPOPCNTDQ. CPUs with POPCNT masked fall back safely.
- Keep C++ coordinates as `int64_t`, with direct integer parsing, scientific-notation compatibility and ordinary genomic range validation.
- Keep R POS columns consistently numeric; remove the oversized-coordinate text conversion and arbitrary decimal parser. Validate ordinary MAP/BIM positions before numeric coercion.
- Add reproducible native/x86 assembly probes, feature-selection regression checks, ordinary coordinate validation tests.
- Record the review evidence and hardware-validation limits in `tests/assembly/README.md`.

## 2.0.0 — rewrite branch, updated 2026-11-07

### Computation and performance

- Share the C++ numerical engine between the CLI and R package.
- Replace repeated per-pair genotype scans with packed category bitsets when the dataset is large enough.
- Add ARM64 NEON, x86 AVX2 and AVX-512 VPOPCNTDQ implementations with runtime selection and POPCNT/SWAR fallbacks.
- Reuse complete-locus marginals and bounded count-log/factorial caches for MI/RMI; cache half-integer logarithms for suitable missing-data workloads.
- Use four intersections for complete pairs, reuse packed diagonals for summaries and write each matrix element once.
- Accelerate concave LD likelihoods with safeguarded Newton iterations. Compare stationary points and endpoints for potentially multimodal cases, correcting wrong-mode results from the original optimizer.
- Cache distinct R genotype strings using reusable open-addressed storage. Optimize compatibility parsing, allele counting, frequency calculation and buffered numeric output.

### Input and workflows

- Use one HTSlib reader for VCF, compressed VCF and BCF in R and C++.
- Return integer genotype matrices directly from R readers, avoiding manual transposition and conversion.
- Preserve PLINK family/individual identity, separator escaping and leading zeros; reject duplicate sample identities.
- Reject fractional or invalid numerical genotypes before native integer conversion.
- Index decay windows by chromosome and coordinate, including unsorted input and large positions.
- Add bounded R pair-file output, double-precision pair text and cancellation checkpoints.
- Reuse cached LD matrices in heatmaps and pair extraction, checking named variant order.

### Visualization

- Add an independent continuous `color_legend()` grid component without ComplexHeatmap.
- Support horizontal/vertical bars, titles, exact tick placement and automatic label spacing.
- Compact the circos canvas, shift the circle slightly left and place the default legend at the upper right.
- Accept native metadata and character pair tables directly in circos plotting.

### Compatibility changes

- Replace the flat CLI with `calculate`, `convert` and `summarize` subcommands. Remove old flag aliases.
- Export integer genotype tables with `convert --format genotype` as `.genotype.tsv`, distinct from valid VCF output; write frequency summaries as `.freq.tsv`.
- Require C++17 in CMake and R build configurations.
- Keep existing R function names, but return direct integer data from `vcf2matrix()` and by default from `ped2matrix()`; request `genotype = "char"` for PED allele strings.
- Preserve numeric genotype allele orientation in `codegeno()`; do not recode existing numeric input against its major allele.
- LD values may change because the optimizer now selects the appropriate likelihood solution. Allele orientation still matters for signed D.

### Documentation and validation

- Translate maintained project READMEs to English and document installation, migration, schemas, numerical conventions and current limitations.
- Record direct comparisons against the pre-rewrite 1.x implementation, rather than multiplying intermediate speedups.
- Include R workflow fixtures in the source package for independent package checks.
- Add numerical, SIMD, utility, streaming, identity, legend and interruption regressions.
- Add an independent base-R LD optimizer oracle, raw genotype export, all stationary/boundary candidate checks and grid-refined checks, with 37 recorded cases including ambiguous signs.
- Add a documented runnable adaptation of the 2023 paper example with independent MI/RMI and LD checks.

Benchmark conditions and measurements are in [docs/performance.md](docs/performance.md). The AVX-512 path still requires execution on supporting hardware.
