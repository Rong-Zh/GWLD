# Genotype-encoding benchmark reference

`codeGenotype.cpp` retains the 1.x implementation from commit `c8079ab` before the 2.0 encoding
optimization; Chinese comments have been translated into English without
changing the code. It is used only by the optional `tests/benchmarks/codegeno.R`
benchmark through `Rcpp::sourceCpp`; it is not compiled into the GWLD package.
