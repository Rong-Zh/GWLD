# Reproduce Box 1 and validate the measures discussed in Zhang et al. (2023).
# https://doi.org/10.1093/g3journal/jkad154
# This is a GWLD 2.0 adaptation, not an unpublished original validation script.
# Usage: Rscript tests/paper/paper_2023.R [output_directory] [cores]

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  root <- dirname(dirname(dirname(normalizePath(script, mustWork = TRUE))))
  out <- if (length(args)) args[1L] else file.path(root, "build", "paper-validation")
  cores <- if (length(args) > 1L) suppressWarnings(as.integer(args[2L])) else 2L
  if (length(args) > 2L || is.na(cores) || cores < 1L)
    stop("Usage: Rscript tests/paper/paper_2023.R [output_directory] [positive_cores]")
  required <- c("GWLD", "genetics", "infotheo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
  if (utils::packageVersion("GWLD") < "2.0.0") stop("This adaptation requires GWLD >= 2.0.0")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)
  write_tsv <- function(x, name) utils::write.table(x, file.path(out, name),
    sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

  env <- new.env()
  utils::data("duck", package = "GWLD", envir = env)
  duck <- env$duck
  g <- GWLD::codegeno(duck$SNP$genotype, sep = "/")
  info <- duck$SNP$info
  cat("GWLD", as.character(utils::packageVersion("GWLD")), "from", find.package("GWLD"), "\n")
  cat(nrow(g), "samples;", ncol(g), "variants; chromosome(s)",
      paste(unique(info[, 1L]), collapse = ","), "\n")
  methods <- c("MI", "RMI", "r^2", "D", "D'")
  matrices <- setNames(lapply(methods, function(m) GWLD::GWLD(g, method = m, cores = 1L)), methods)
  for (m in methods) {
    stopifnot(identical(matrices[[m]], GWLD::GWLD(g, method = m, cores = cores)))
  }
  stopifnot(isTRUE(all.equal(matrices$MI, GWLD::MI(g), check.attributes = FALSE)),
            isTRUE(all.equal(matrices$RMI, GWLD::RMI(g), check.attributes = FALSE)),
            isTRUE(all.equal(matrices[["r^2"]], GWLD::LD(g, method = "r^2"), check.attributes = FALSE)))

  # Independent R expression for the RMI approximation, equation 6 and the
  # Newman et al. approximation used by GWLD. Preserve whole-locus categories.
  reference_rmi <- function(a, b) {
    tab <- table(factor(a, levels = sort(unique(stats::na.omit(a)))),
                 factor(b, levels = sort(unique(stats::na.omit(b)))))
    n <- sum(tab)
    if (!n) return(NA_real_)
    r <- nrow(tab); s <- ncol(tab)
    aa <- rowSums(tab); bb <- colSums(tab)
    w <- n / (n + 0.5 * r * s)
    x <- (1 - w) / r + w * aa / n
    y <- (1 - w) / s + w * bb / n
    nu <- (s + 1) / (s * sum(x^2)) - 1 / s
    mu <- (r + 1) / (r * sum(y^2)) - 1 / r
    log_omega <- (r - 1) * (s - 1) * log(n + 0.5 * r * s) +
      0.5 * (r + nu - 2) * sum(log(y)) + 0.5 * (s + mu - 2) * sum(log(x)) +
      0.5 * (lgamma(mu * r) + lgamma(nu * s) -
             r * (lgamma(s) + lgamma(mu)) - s * (lgamma(r) + lgamma(nu)))
    (lgamma(n + 1) + sum(lgamma(tab + 1)) - sum(lgamma(aa + 1)) -
       sum(lgamma(bb + 1)) - log_omega) / (n * log(2))
  }
  pairs <- which(upper.tri(matrices$MI), arr.ind = TRUE)
  references <- matrix(NA_real_, nrow(pairs), length(methods), dimnames = list(NULL, methods))
  # A is the major allele after codegeno(); this aligns signed D with genetics.
  genotypes <- lapply(seq_len(ncol(g)), function(j)
    genetics::genotype(c("A/A", "A/G", "G/G")[g[, j] + 1L], reorder = "no"))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]; j <- pairs[k, 2L]
    complete <- stats::complete.cases(g[, c(i, j)])
    # infotheo reports nats; GWLD reports bits.
    references[k, "MI"] <- infotheo::mutinformation(g[complete, i], g[complete, j]) / log(2)
    references[k, "RMI"] <- reference_rmi(g[, i], g[, j])
    ref <- genetics::LD(genotypes[[i]], genotypes[[j]])
    direction <- if (mean(g[, i], na.rm = TRUE) > 1) -1 else 1
    direction <- direction * if (mean(g[, j], na.rm = TRUE) > 1) -1 else 1
    references[k, c("r^2", "D", "D'")] <- c(ref[["R^2"]], direction * ref$D, ref[["D'"]])
  }
  checks <- lapply(methods, function(m) {
    actual <- matrices[[m]][pairs]
    expected <- references[, m]
    errors <- abs(actual - expected)
    tolerance <- if (m %in% c("MI", "RMI")) 1e-10 else 1e-3
    passed <- identical(is.na(actual), is.na(expected)) &&
      all(is.finite(errors[!is.na(errors)])) && all(errors <= tolerance, na.rm = TRUE)
    write_tsv(data.frame(ID_1 = info[pairs[, 1L], 3L], ID_2 = info[pairs[, 2L], 3L],
      GWLD = actual, reference = expected, abs_error = errors),
      paste0("pairs-", match(m, methods), ".tsv"))
    data.frame(method = m, reference = if (m == "MI") "infotheo (bits)" else
      if (m == "RMI") "independent R expression" else "genetics",
      pairs = length(actual), max_abs_error = max(c(0, errors), na.rm = TRUE),
      tolerance = tolerance, passed = passed)
  })
  checks <- do.call(rbind, checks)
  write_tsv(checks, "checks.tsv")
  print(checks, row.names = FALSE)

  # Reproduce Box 1's two ways of calculating decay and circos.
  decay <- GWLD::calc_decay(g, info, method = "RMI", cores = cores)
  circos <- GWLD::calc_circos(g, info, method = "RMI", cores = cores, threshold = 0.2)
  stopifnot(identical(decay, GWLD::decay(matrices$RMI, info)),
            identical(circos, GWLD::circos(matrices$RMI, info, threshold = 0.2)))
  write_tsv(decay, "decay.tsv")
  write_tsv(circos, "circos.tsv")
  values <- vapply(matrices[c("MI", "RMI", "r^2")], function(m) m[pairs], numeric(nrow(pairs)))
  correlations <- stats::cor(values, use = "pairwise.complete.obs")
  write_tsv(data.frame(method = rownames(correlations), correlations, check.names = FALSE), "correlations.tsv")
  print(correlations)
  saveRDS(list(Info = info, Genotype = g, matrices = matrices, reference = references,
               pairs = pairs, decay = decay, circos = circos), file.path(out, "results.rds"))
  capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))

  plot_png <- function(name, draw) {
    grDevices::png(file.path(out, name), width = 2400, height = 1600, res = 160)
    on.exit(grDevices::dev.off())
    draw()
  }
  for (m in c("MI", "RMI", "r^2")) {
    plot_png(paste0("heatmap-", if (m == "r^2") "r2" else m, ".png"), function()
      GWLD::HeatMap(LDmatrix = matrices[[m]], method = m, SnpPosition = info[, 1:2], SnpName = info[, 3L],
                    cores = cores, color = "YellowTored"))
  }
  # Box 1 plots the stored genome-wide results, not the chromosome-26 subset.
  plot_png("circos-stored-duck.png", function() {
    GWLD::circos.ideogram(duck$Circos$chr)
    GWLD::circos.linksnp(duck$Circos$linkdata, bg.col = grDevices::rainbow(29))
  })
  cat("Decay pairs:", nrow(decay), "; circos pairs:", nrow(circos), "\n")
  if (!all(checks$passed)) stop("Validation failed; see ", file.path(out, "checks.tsv"))
  cat("All paper validation checks passed. Results:", out, "\n")
}

main()
