# Independent base-R validation of the fixed-marginal haplotype likelihood.
# No GWLD library or production likelihood/solver is used by this oracle.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript tests/ld/validate_ld.R [exporter] [output_directory]")
exporter <- if (length(args)) args[1L] else "build/gwld_ld_validation_export"
out <- if (length(args) > 1L) args[2L] else "build/ld-validation"
exporter <- normalizePath(exporter, mustWork = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
status <- system2(exporter, shQuote(c(file.path(out, "cases.tsv"), file.path(out, "genotypes.tsv"))))
if (status != 0L) stop("Case exporter failed")
cases <- read.delim(file.path(out, "cases.tsv"), check.names = FALSE)
raw <- read.delim(file.path(out, "genotypes.tsv"))

# Polynomial coefficients are in increasing powers. This implementation uses
# base R polyroot (a different root algorithm from production root isolation).
multiply <- function(a, b) {
  result <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) for (j in seq_along(b)) result[i+j-1L] <- result[i+j-1L] + a[i]*b[j]
  result
}
add <- function(a, b) {
  length(a) <- length(b) <- max(length(a), length(b))
  a[is.na(a)] <- 0; b[is.na(b)] <- 0
  a+b
}
derivative <- function(a) if (length(a) < 2L) 0 else a[-1L]*seq_len(length(a)-1L)

checks <- roots_output <- refinements <- vector("list", nrow(cases))
for (row in seq_len(nrow(cases))) {
  case <- cases[row, ]
  observations <- raw[raw$case_id == case$case_id, ]
  a <- observations$g1; b <- observations$g2
  stopifnot(length(a) > 0L, all(a %in% -1:2), all(b %in% -1:2))
  # Rebuild both marginals and contingency counts from the raw observations.
  p <- mean((2-a[a >= 0])/2); q <- mean((2-b[b >= 0])/2)
  complete <- a >= 0 & b >= 0
  table <- table(factor(a[complete], levels = 0:2), factor(b[complete], levels = 0:2))
  counts <- as.numeric(t(table))
  stopifnot(abs(p-case$p) < 1e-14, abs(q-case$q) < 1e-14,
            identical(counts, as.numeric(case[paste0("c", 0:8)])))
  lo <- max(0, p+q-1); hi <- min(p, q); span <- hi-lo
  stopifnot(span > 0, sum(counts) > 0)
  probabilities <- function(x) {
    h <- c(x, p-x, q-x, 1-p-q+x)
    c(h[1]^2, 2*h[1]*h[2], h[2]^2,
      2*h[1]*h[3], 2*(h[1]*h[4]+h[2]*h[3]), 2*h[2]*h[4],
      h[3]^2, 2*h[3]*h[4], h[4]^2)
  }
  # Full nine-genotype likelihood, including the heterozygote factors of two.
  nll <- function(x) {
    if (x < lo || x > hi) return(Inf)
    probs <- probabilities(x); used <- counts > 0
    if (any(probs[used] <= 0)) return(Inf)
    -sum(counts[used]*log(probs[used]))
  }
  # Clearing positive denominator factors from d log L / dx gives degree <=5.
  # Scale to u=(x-lo)/span before forming the polynomial.
  weights <- c(2*counts[1]+counts[2]+counts[4], 2*counts[3]+counts[2]+counts[6],
               2*counts[7]+counts[4]+counts[8], 2*counts[9]+counts[6]+counts[8], counts[5])
  factors <- list(c(lo,span), c(p-lo,-span), c(q-lo,-span), c(1-p-q+lo,span),
    c(lo*(1-p-q+lo)+(p-lo)*(q-lo), span*(4*lo+1-2*p-2*q), 2*span^2))
  factors <- lapply(factors, function(x) x/max(abs(x)))
  active <- which(weights > 0)
  polynomial <- 0
  for (i in active) {
    term <- derivative(factors[[i]])
    for (j in setdiff(active, i)) term <- multiply(term, factors[[j]])
    polynomial <- add(polynomial, weights[i]*term)
  }
  while (length(polynomial) > 1L && tail(polynomial, 1L) == 0) polynomial <- head(polynomial, -1L)
  polynomial <- polynomial/max(abs(polynomial))
  roots <- if (length(polynomial) > 1L) polyroot(polynomial) else complex()
  real <- Re(roots[abs(Im(roots)) <= 1e-7*pmax(1, abs(Re(roots)))])
  real <- real[real >= -1e-9 & real <= 1+1e-9]
  x <- sort(unique(c(lo, hi, lo+span*pmax(0,pmin(1,real)))))
  score <- vapply(x, nll, numeric(1))
  best_score <- min(score)
  roots_output[[row]] <- data.frame(case_id = case$case_id, candidate_x = x, negative_log_likelihood = score,
                                   boundary = x == lo | x == hi)

  # A second numerical check uses no stationary-polynomial roots: dense grid
  # followed by independent optimize() calls around every grid-local minimum.
  grid <- seq(lo, hi, length.out = 10001L)
  grid_score <- vapply(grid, nll, numeric(1))
  centers <- which(grid_score[2:10000] <= grid_score[1:9999] & grid_score[2:10000] <= grid_score[3:10001])+1L
  refined <- c(lo, hi)
  for (k in centers) refined <- c(refined, optimize(nll, c(grid[k-1L], grid[k+1L]), tol = 1e-12)$minimum)
  refined_score <- vapply(refined, nll, numeric(1))
  refinements[[row]] <- data.frame(case_id = case$case_id, candidate_x = refined, negative_log_likelihood = refined_score)

  new_x <- case$new_D+p*q; old_x <- case$old_D+p*q
  new_score <- nll(new_x); old_score <- nll(old_x)
  nearest <- min(abs(x-new_x))
  # Check reported LD measures independently from the selected haplotype frequency.
  expected_r2 <- case$new_D^2/(p*q*(1-p)*(1-q))
  bound <- if (case$new_D > 0) min(p*(1-q),q*(1-p)) else max(-p*q,-(1-p)*(1-q))
  expected_Dprime <- case$new_D/bound
  oracle_gap <- new_score-best_score
  grid_gap <- new_score-min(refined_score)
  passed <- is.finite(new_score) && abs(oracle_gap) <= 1e-8 && abs(grid_gap) <= 1e-8 &&
    new_score <= old_score+1e-8 && nearest <= 1e-7 &&
    abs(case$new_r2-expected_r2) <= 1e-12 && abs(case$new_Dprime-expected_Dprime) <= 1e-12
  checks[[row]] <- data.frame(case_id=case$case_id, dataset=case$dataset, selection=case$selection,
    old_D=case$old_D, new_D=case$new_D, abs_D_difference=abs(case$new_D-case$old_D),
    relative_D_difference=if (case$old_D == 0) NA_real_ else abs((case$new_D-case$old_D)/case$old_D),
    old_log_likelihood=-old_score, new_log_likelihood=-new_score,
    log_likelihood_gain=old_score-new_score, oracle_nll_gap=oracle_gap,
    grid_refined_nll_gap=grid_gap, nearest_candidate_distance=nearest,
    stationary_candidates=sum(x > lo & x < hi), passed=passed)
}
checks <- do.call(rbind, checks)
write_tsv(checks, "independent_checks.tsv")
write_tsv(do.call(rbind, roots_output), "stationary_points.tsv")
write_tsv(do.call(rbind, refinements), "grid_refinement.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(checks[, c("case_id", "dataset", "log_likelihood_gain", "oracle_nll_gap", "passed")], row.names=FALSE)
cat("Max distance to independent stationary/boundary candidate:", max(checks$nearest_candidate_distance), "\n")
cat("Max absolute gap to independent best likelihood:", max(abs(checks$oracle_nll_gap)), "\n")
if (!all(checks$passed)) stop("Independent LD validation failed; inspect ", out)
cat("All", nrow(checks), "LD cases passed. Evidence:", out, "\n")
