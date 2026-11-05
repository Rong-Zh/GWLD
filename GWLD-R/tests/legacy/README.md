# Genotype-encoding benchmark reference

`codeGenotype.cpp` is the unchanged implementation before the 2.0 encoding
optimization. It is used only by the optional `tests/benchmarks/codegeno.R`
benchmark through `Rcpp::sourceCpp`; it is not compiled into the GWLD package.
