# Optional benchmark: run from repository root after installing GWLD.
# Pass --source to compile the current encoder directly without installing GWLD.
# Use --iterations=20 to resolve short calls on coarse elapsed-time clocks.
# Three alternating batches, medians per call; compare outputs after each batch.
arguments <- commandArgs(trailingOnly=TRUE)
iteration_option <- grep("^--iterations=",arguments,value=TRUE)
iterations <- if (length(iteration_option)) as.integer(sub("^--iterations=","",iteration_option)) else 1L
stopifnot(length(iterations)==1L,!is.na(iterations),iterations>=1L)
if ("--source" %in% arguments) {
  current <- new.env()
  Rcpp::sourceCpp("GWLD-R/src/codeGenotype.cpp",env=current)
  recode <- current$recode
  stats_freq <- current$stats_freq
  optimized_matrix <- current$code_Mat
} else {
  library(GWLD)
  optimized_matrix <- getFromNamespace("code_Mat","GWLD")
}
legacy <- new.env()
Rcpp::sourceCpp("GWLD-R/tests/legacy/codeGenotype.cpp",env=legacy)
set.seed(2026)
g <- matrix(sample(c("A/A","A/G","G/A","G/G",NA_character_),2000000,TRUE),nrow=2000)
v <- as.vector(g)
cases <- list(
  recode=list(old=function() legacy$recode(v,"/"),new=function() recode(v,"/")),
  stats_freq=list(old=function() legacy$stats_freq(v,"/"),new=function() stats_freq(v,"/")),
  code_Mat=list(old=function() legacy$code_Mat(g,"/"),new=function() optimized_matrix(g,"/")))
cat("method\tlegacy_seconds\toptimized_seconds\tspeedup\n")
for (name in names(cases)) {
  times <- list(old=numeric(3),new=numeric(3))
  for (repeat_id in 1:3) {
    answers <- list()
    for (kind in if(repeat_id%%2) c("old","new") else c("new","old")) {
      gc()
      times[[kind]][repeat_id] <- system.time({
        for (iteration in seq_len(iterations)) answers[[kind]] <- cases[[name]][[kind]]()
      })[["elapsed"]]/iterations
    }
    stopifnot(isTRUE(all.equal(answers$old,answers$new,check.attributes=TRUE)))
  }
  old <- median(times$old); fresh <- median(times$new)
  cat(name,old,fresh,old/fresh,sep="\t"); cat("\n")
}
