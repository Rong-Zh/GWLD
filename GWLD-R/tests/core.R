# Runs both after sourceCpp(src/method.cpp) and against an installed package.
if (!exists("RMIC", mode = "function")) {
  library(GWLD)
  for (name in c("RMIC", "RMIC_Mat", "MIC", "MIC_Mat", "LDC", "LDC_Mat", "gwld"))
    assign(name, getFromNamespace(name, "GWLD"))
}
near <- function(a, b, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)))
}
fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
reference <- function(a, b, reduced = FALSE) {
  a[a == -1L] <- NA_integer_
  b[b == -1L] <- NA_integer_
  m <- table(factor(a, levels = sort(unique(na.omit(a)))),
             factor(b, levels = sort(unique(na.omit(b)))))
  n <- sum(m)
  if (!n) return(NA_real_)
  if (!reduced) {
    h <- function(x) { p <- x[x > 0] / n; -sum(p * log2(p)) }
    return(h(rowSums(m)) + h(colSums(m)) - h(m))
  }
  R <- nrow(m); S <- ncol(m)
  I <- lgamma(n + 1) + sum(lgamma(m + 1)) - sum(lgamma(rowSums(m) + 1)) - sum(lgamma(colSums(m) + 1))
  w <- n / (n + 0.5 * R * S)
  x <- (1-w) / R + w * rowSums(m) / n
  y <- (1-w) / S + w * colSums(m) / n
  nu <- (S+1) / (S*sum(x*x)) - 1/S
  mu <- (R+1) / (R*sum(y*y)) - 1/R
  omega <- (R-1)*(S-1)*log(n+0.5*R*S) + 0.5*(R+nu-2)*sum(log(y)) +
    0.5*(S+mu-2)*sum(log(x)) +
    0.5*(lgamma(mu*R)+lgamma(nu*S)-R*(lgamma(S)+lgamma(mu))-S*(lgamma(R)+lgamma(nu)))
  (I-omega) / (n*log(2))
}
set.seed(2026)
g <- matrix(sample(c(0L, 1L, 2L, NA_integer_), 1200, replace = TRUE), nrow = 100)
for (i in seq_len(ncol(g)-1L)) for (j in seq.int(i+1L, ncol(g))) {
  near(MIC(g[,i], g[,j]), reference(g[,i], g[,j]))
  near(RMIC(g[,i], g[,j]), reference(g[,i], g[,j], TRUE))
}
for (method in c("RMI", "MI", "r^2", "D", "D'")) {
  one <- gwld(g, method, 1L)
  two <- gwld(g, method, 2L)
  stopifnot(identical(one, two), all(is.na(one[lower.tri(one, diag = TRUE)])))
}
near(MIC_Mat(g, 2L), gwld(g, "MI", 2L))
near(RMIC_Mat(g, 2L), gwld(g, "RMI", 2L))
near(LDC_Mat(g, "r^2", 2L), gwld(g, "r^2", 2L))
# Exercise the packed backend and RMI lookup path above the small-cohort cutoff.
large <- rbind(g,g,g)
for (method in c("MI","RMI","r^2","D","D'")) {
  packed <- gwld(large,method,1L)
  stopifnot(identical(packed,gwld(large,method,2L)))
  for (i in seq_len(ncol(large)-1L)) for(j in seq.int(i+1L,ncol(large))) {
    scalar <- switch(method,MI=MIC(large[,i],large[,j]),RMI=RMIC(large[,i],large[,j]),
                     LDC(large[,i],large[,j],method))
    near(packed[i,j],scalar)
  }
}
a <- c(0L, 0L, 1L, 1L, 2L, 2L)
near(LDC(a, a, "r^2"), 1)
near(LDC(a, 2L-a, "D'"), 1)
stopifnot(is.na(MIC(a, rep(NA_integer_, 6))), is.na(LDC(a, rep(0L, 6), "r^2")))
fails(MIC(a, a[-1]))
fails(MIC(c(3L, 1L), c(1L, 1L)))
fails(gwld(g, "unknown", 1L))
fails(gwld(g, "MI", 0L))
stopifnot(identical(dim(gwld(matrix(integer(),0,4),"RMI")),c(4L,4L)),
          all(is.na(gwld(matrix(integer(),0,4),"MI"))))
info <- cbind(rep(c("1", "2"), 6), c(2000, 10, 1000, 900000, 3000, 5000000000, 0, 11, 1000, 12, 4000, 13), paste0("s", 1:12))
for (method in c("RMI", "MI", "r^2", "D", "D'")) {
  m <- gwld(g, method)
  near(calc_decay(g, info, method, 1L, 1L), calc_decay(g, info, method, 2L, 1L))
  near(calc_decay(g, info, method, 2L, 1L), decay(m, info, 1L))
  near(calc_circos(g, info, method, 2L, 0.1), circos(m, info, 0.1))
}
empty <- info
empty[,1] <- as.character(seq_len(nrow(empty)))
stopifnot(identical(dim(calc_decay(g, empty)), c(0L, 4L)))
stopifnot(identical(dim(calc_circos(g, info, threshold = 100)), c(0L, 7L)))
fails(calc_decay(g, info[-1,,drop=FALSE]))
fails(calc_decay(g, info, maxdist = -1L))
fails(decay(matrix(0, 2, 3), info))
if ("package:GWLD" %in% search()) {
  near(GWLD(a, a, "MI"), MI(a, a))
  near(GWLD(a, a, "RMI"), RMI(a, a))
  near(GWLD(as.data.frame(g)), GWLD(g))
  bed <- tempfile(fileext = ".bed")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, 0xe4)), bed)
  near(GWLD:::read_bed(bed, 4L)[,1], c(2L, NA_integer_, 1L, 0L))
  fails(GWLD:::read_bed(bed, 0L))
  writeBin(as.raw(c(0x6c, 0x1b, 0x00, 0xe4)), bed)
  fails(GWLD:::read_bed(bed, 4L))
  unlink(bed)
}
cat("R numerical, parallel, metadata and empty-result tests passed\n")
