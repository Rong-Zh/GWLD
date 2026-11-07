# Independent verification of the LD optimizer

This directory makes the LD accuracy claim reproducible. It separates a mathematical characterization of the optimum from finite-precision evidence for the implementation. The tests validate the fixed-marginal likelihood used by GWLD; they do not establish that this statistical model is appropriate for every biological dataset.

## Run

From the repository root, with the C++ build dependencies and base R installed:

```bash
cmake -S GWLD-C++ -B build -DGWLD_BUILD_BENCHMARKS=ON
cmake --build build --target gwld_ld_validation_export -j
Rscript --vanilla tests/ld/validate_ld.R \
  build/gwld_ld_validation_export build/ld-validation
```

Use the platform/compiler options from the [installation guide](../../README.md) when configuring a fresh build. The exporter requires OpenMP to compile the frozen original implementation, even though it runs that implementation with one thread. It does not require an installed GWLD R package. The independent oracle requires only base R.

When benchmarks, testing and `Rscript` are available, CMake also registers this check:

```bash
ctest --test-dir build -R ld_independent --output-on-failure
```

## What constitutes evidence

For fixed allele frequencies `p` and `q`, let the four haplotype frequencies be

```text
h = (x, p-x, q-x, 1-p-q+x)
max(0, p+q-1) <= x <= min(p, q)
D = x - p*q
r² = D² / (p*q*(1-p)*(1-q))
```

Under the unphased diploid model, the 3×3 genotype probabilities are

```text
          h1²                 2*h1*h2              h2²
          2*h1*h3             2*(h1*h4+h2*h3)      2*h2*h4
          h3²                 2*h3*h4              h4²
```

For observed genotype counts `n[i,j]`, the objective is

```text
log L(x) = sum(n[i,j] * log(P[i,j](x)))
```

Zero-count terms contribute zero. A zero probability with positive count gives negative-infinite log likelihood. A multinomial coefficient, constant in `x`, may be omitted without changing the optimum.

The objective is differentiable in the interior of the feasible interval. Its derivative is a sum of rational terms from four linear factors and one quadratic factor. Clearing the active positive denominators gives a polynomial of degree at most five. Thus, for a nonconstant objective on a nondegenerate interval, a global maximum must be at a feasible boundary or at an interior real root of that derivative polynomial. Comparing every such candidate characterizes the global solution in exact arithmetic. Tied global maxima must all be considered valid solutions.

This argument explains what a correct solver must find. It is not, by itself, a formal proof that every floating-point branch in the C++ implementation is correct. The independent numerical checks below exercise the implementation against that characterization.

## Independent implementation

`export_cases.cpp` runs the actual production engine and the frozen pre-rewrite `c8079ab` implementation. It exports the raw genotype observations as well as the old/new LD results. It selects the maximum absolute r², D and D′ differences from each of the three recorded benchmark workloads.

`validate_ld.R` does **not** load GWLD or call its likelihood, counting or root-solving code. It:

1. Rebuilds allele frequencies and all nine genotype counts from the raw exported observations using base R, then checks them against the C++ summaries.
2. Constructs the full nine-genotype probability formula above and evaluates its log likelihood directly.
3. Forms the stationary polynomial independently in R and uses `polyroot`, a different numerical root algorithm from the production root-isolation implementation. It evaluates all feasible real-root candidates and both boundaries.
4. Performs a second search using a 10,001-point grid and `optimize()` around every grid-local minimum of negative log likelihood. This search does not use the polynomial roots.
5. Checks the production result against both searches and independently derives r² and D′ from D.
6. Writes the inputs, candidate values, old/new log likelihoods, numerical gaps and pass/fail results. Failures produce a nonzero exit status.

GWLD retains the original whole-locus allele-frequency convention even when missing observations reduce the complete-pair contingency table. The oracle uses the same stated convention but reconstructs it independently from raw observations. The claim concerns maximizing over `x` conditional on those fixed marginals, not jointly estimating every frequency under a different missing-data model.

## Cases and acceptance criteria

There are 37 cases:

- Nine maximum-difference selections: three methods across the three original-version benchmark datasets. Some selections refer to the same pair; these nine selections are not nine distinct pairs.
- One known multimodal counterexample.
- Perfect coupling, perfect repulsion and tied boundary maxima.
- Twenty-four additional deterministic pairs with missing observations, seed 654321.

The maximum-difference workload generator uses seed 2026 and matches [the recorded performance comparison](../../docs/performance.md).

Each case must satisfy:

- Independent allele frequencies within `1e-14` and exactly matching genotype counts.
- A finite production likelihood, within `1e-8` absolute negative log likelihood of the independent stationary/boundary optimum and the grid-refined optimum.
- Production likelihood no worse than the original result, allowing `1e-8` numerical tolerance.
- Production haplotype frequency within `1e-7` of a feasible independently computed stationary/boundary candidate.
- Independently derived r² and D′ within `1e-12` of the production values.

The tolerances are assertions on explicit quantities, not a blanket assertion that different LD values are interchangeable. Near-zero D values can have large relative differences despite small absolute changes.

## Recorded evidence

On 2026-11-07, macOS arm64, all 37 cases passed:

| Quantity | Observed maximum |
| --- | ---: |
| Distance from production haplotype frequency to an independent candidate | `4.0e-14` |
| Absolute negative-log-likelihood gap to the independent optimum | `1.82e-12` |

For the multimodal counterexample, the observed table is

```text
0    2    0
2  262    5
2    5    2
```

The original result was `D = +0.2333543571`; the production result was `D = -0.2336136616`. The independently evaluated log likelihood increased from approximately `-275.6128404` to `-273.2742052`, a gain of `2.3386352`. The oracle found three interior stationary candidates and selected the same solution as production. Here, a large D difference exposes the original wrong-mode choice.

The all-double-heterozygote case illustrates a different issue: both `D = -0.25` and `D = +0.25` are exact global maxima. The sign is not identifiable from those unphased observations under this model. A sign difference in that case does not establish that either exact boundary solution is wrong; compare likelihood and the set of valid solutions.

Recorded snapshots are in [`recorded/cases.tsv`](recorded/cases.tsv), [`recorded/independent_checks.tsv`](recorded/independent_checks.tsv) and [`recorded/stationary_points.tsv`](recorded/stationary_points.tsv). They contain synthetic data only. Re-running the command regenerates raw genotypes, these tables, grid-refinement candidates and `sessionInfo.txt` in the selected output directory. The snapshots are review evidence, not fixed expected values that replace the independent oracle.

These checks provide reproducible evidence for these inputs. They do not formally certify every possible floating-point input. The broader core tests also cover missing/monomorphic cases, allele reversal, numerical references and additional likelihood cases.
