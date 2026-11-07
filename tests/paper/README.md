# 2023 paper example and numerical validation

`paper_2023.R` adapts Box 1 and the methods described in Zhang et al. (2023) to GWLD 2.0. It uses the packaged `duck$SNP` subset: 500 samples and 60 variants on chromosome 26, giving 1,770 unique variant pairs. It writes numerical checks, pair tables and figures.

## Sources and provenance

- Paper: [GWLD: an R package for genome-wide linkage disequilibrium analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10468308/), [DOI: 10.1093/g3journal/jkad154](https://doi.org/10.1093/g3journal/jkad154).
- `reference/box1_original.txt`: Box 1 extracted from the paper XML. It retains the original quotation marks and input-path placeholders for provenance and is not directly executable R code.
- `reference/jkad154_supplementary_data.zip`: the original supplementary files, including the C++ manual `File_S1_G3-2023-404081.docx` and `Supplemental_Tables_and_Figures_G3-2023-404081.docx`. These files do not contain an independent validation script. SHA-256: `988f8aee2d823d1f1a886ccd663354207373527d59d4e473ca9a758db11e3339`.
- Retrieval endpoints: [Europe PMC full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10468308/fullTextXML) and [supplementary files](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10468308/supplementaryFiles). Retrieved 2026-11-05. The paper excerpt is retained with attribution under CC BY 4.0.

The runnable script is a newly written adaptation, not original validation code extracted from the supplement. It removes automatic remote package installation and supplies chromosome/position metadata to the current heatmap interface instead of the original `SnpPosition = Info$POS` call.

## Run

Install the local `GWLD-R` package and the independent reference packages:

```r
install.packages(c("genetics", "infotheo"))
```

From the repository root:

```bash
Rscript tests/paper/paper_2023.R
```

The default output directory is `build/paper-validation`, with two threads. Supply an output directory and thread count explicitly if desired:

```bash
Rscript tests/paper/paper_2023.R build/paper-validation-4 4
```

For a project-local R library:

```bash
R_LIBS="$PWD/build/R-library" Rscript tests/paper/paper_2023.R
```

Use one thread with a serial GWLD build. The script prints the loaded GWLD version and library path; it does not install or update packages. Missing reference packages or failed assertions produce a nonzero exit status.

## Checks and tolerances

| Check | Reference | Absolute tolerance |
| --- | --- | --- |
| MI | `infotheo::mutinformation`, converted from natural-log units to bits | `1e-10` |
| RMI | Independent R expression for the approximation used by GWLD | `1e-10` |
| r², D, D′ | `genetics::LD`, with D aligned to the major-allele orientation | `1e-3` |
| Parallel execution | One thread versus the requested thread count, all five methods | Exact equality |
| R entry points | `GWLD` versus `MI`, `RMI` and `LD` | Default R numerical comparison |
| Pair extraction | Direct calculations versus reuse of the RMI matrix | Exact equality |

The LD tolerance accommodates optimizer convergence differences between `genetics` and GWLD's explicit comparison of feasible solutions. It does not imply bitwise equality. MI uses complete observations for each pair; RMI retains GWLD's whole-locus category convention.

The script also writes a Pearson correlation matrix for MI/RMI/r² and three heatmaps. `circos-stored-duck.png` uses the packaged historical genome-wide circos data. It is not a genome-wide calculation from the 60 chromosome-26 variants; that subset has no cross-chromosome pairs.

## Outputs

| File | Contents |
| --- | --- |
| `checks.tsv` | Reference, number of pairs, maximum error, tolerance and pass/fail result. |
| `pairs-1.tsv` through `pairs-5.tsv` | Pair-level results and reference errors for MI, RMI, r², D and D′, respectively. |
| `correlations.tsv` | Correlations among the three summary measures. |
| `decay.tsv`, `circos.tsv` | RMI pair outputs. |
| `results.rds` | Genotypes, metadata, five result matrices and reference values. |
| `heatmap-*.png`, `circos-stored-duck.png` | Figures. |
| `sessionInfo.txt` | R and package versions. |

## Scope and observed results

The complete duck dataset in the paper has 542 samples and 40,601 SNPs. This script uses a packaged subset and does not attempt to reproduce the exact correlations in supplementary Table S1. It does not run PLINK, TASSEL or the full human, pig and maize datasets, and it is not a reproduction of all paper results or paper benchmarks.

On 2026-11-07, macOS arm64, R 4.6.1 and the locally installed GWLD 2.0.0 passed all five numerical checks, parallel checks and both pair-extraction paths. Maximum absolute errors were approximately `7.33e-15` for MI, `7.88e-15` for RMI, `3.70e-4` for r², `7.64e-5` for D and `6.48e-4` for D′.
