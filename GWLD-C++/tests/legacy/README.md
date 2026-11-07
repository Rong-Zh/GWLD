# Historical references and utility compatibility helpers

`methods.cpp` retains the v1 calculation code from commit `c8079ab`, apart
from translated comments and qualified references to its historical utilities.
`utils_reference.cpp` and `utils_reference.h` freeze the original utility
algorithms in `legacy_reference`. The optional `gwld_benchmark` links these
frozen utilities so changes to the optimized helpers do not alter its baseline.
The benchmark accepts samples, variants, threads and missing percentage, covers
all five measures, warms both paths and reports three-run alternating medians.
MI/RMI outputs must agree within 1e-10; LD differences are reported because
the production optimizer intentionally corrects the original solution search.
Measured 1.x comparisons are recorded in [the performance report](../../../docs/performance.md).

`utils.cpp` and `utils.h` now contain the optimized compatibility helpers:
const-reference inputs, delimiter lookup, direct file-field parsing, buffered
numeric writing, shared vectorizable genotype counts, and checked metadata
selection. Valid rectangular output retains the original trailing delimiters
and six significant digits. Recent C++17 libraries use `to_chars`; older
language/library configurations retain the classic-locale stream formatter.
Empty output is supported, invalid input throws, and the final chromosome
is no longer lost when it starts on the last row. SNP metadata selection
returns exactly the requested columns, including selections other than three.

These compatibility helpers are not linked into the production CLI. Production
calculations use `GWLD-R/inst/include/gwld`; production PLINK field parsing
uses `include/dataset.hpp`. The optimized genotype summary is shared by the
helpers, the CLI frequency command and the R/C++ calculation interfaces.

`gwld_utils_test` compares valid results and file bytes against the frozen
reference and checks corrected boundary cases. With benchmarks enabled,
`gwld_utils_benchmark [output_prefix]` alternates three timed runs after warmup,
reports medians, and verifies each result. It measures splitting, production
whitespace parsing, reading, frequencies, allele counting and matrix writing.
