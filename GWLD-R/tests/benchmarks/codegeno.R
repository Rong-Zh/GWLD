# Optional benchmark: run from repository root after installing GWLD.
# Three alternating runs, medians; all timed calls also compare full outputs.
library(GWLD)
legacy <- new.env()
Rcpp::sourceCpp("GWLD-R/tests/legacy/codeGenotype.cpp",env=legacy)
set.seed(2026)
g <- matrix(sample(c("A/A","A/G","G/A","G/G",NA_character_),2000000,TRUE),nrow=2000)
v <- as.vector(g)
cases <- list(
  recode=list(old=function() legacy$recode(v,"/"),new=function() recode(v,"/")),
  stats_freq=list(old=function() legacy$stats_freq(v,"/"),new=function() stats_freq(v,"/")),
  code_Mat=list(old=function() legacy$code_Mat(g,"/"),new=function() GWLD:::code_Mat(g,"/")))
cat("method\tlegacy_seconds\toptimized_seconds\tspeedup\n")
for (name in names(cases)) {
  times <- list(old=numeric(3),new=numeric(3))
  for (repeat_id in 1:3) {
    answers <- list()
    for (kind in if(repeat_id%%2) c("old","new") else c("new","old")) {
      gc()
      times[[kind]][repeat_id] <- system.time(answers[[kind]] <- cases[[name]][[kind]]())[["elapsed"]]
    }
    stopifnot(isTRUE(all.equal(answers$old,answers$new,check.attributes=TRUE)))
  }
  old <- median(times$old); fresh <- median(times$new)
  cat(name,old,fresh,old/fresh,sep="\t"); cat("\n")
}
