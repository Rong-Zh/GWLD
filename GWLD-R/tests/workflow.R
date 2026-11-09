# Run from the repository root or through R CMD check after installing the package.
library(GWLD)
near <- function(a, b) stopifnot(isTRUE(all.equal(a, b, check.attributes = FALSE)))
fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))

fixtures <- if (dir.exists("GWLD-R/tests/fixtures")) "GWLD-R/tests/fixtures" else "fixtures"

# Reader outputs are interchangeable, including explicit PED character output.
ped <- ped2matrix(file.path(fixtures, "example-2"))
chars <- ped2matrix(file.path(fixtures, "example-2"), genotype = "char")
stopifnot(is.integer(ped$Genotype), identical(ped$Info, chars$Info))
stopifnot(identical(ped$Genotype, codegeno(chars$Genotype)))
vcf <- vcf2matrix(file.path(fixtures, "example-3.vcf"))
bed <- bed2matrix(file.path(fixtures, "example"))
for (input in list(ped, vcf, bed)) {
  stopifnot(is.integer(input$Genotype), identical(codegeno(input$Genotype), input$Genotype))
}
# An encoded ALT allele count must never be encoded again relative to its major allele.
fixed <- matrix(c(2L, 2L, 1L, NA_integer_), ncol = 1, dimnames = list(NULL, "alt-major"))
stopifnot(identical(codegeno(fixed), fixed), identical(codegeno(fixed * 1), fixed))
for (bad in c(3, -1, 0.5, Inf)) fails(codegeno(matrix(bad, 2, 2)))

# Missing and multi-character PED alleles, including a frequency tie.
root <- tempfile("gwld-workflow-")
dir.create(root)
prefix <- file.path(root, "alleles")
writeLines(c("1 a 0 1", "1 b 0 2", "1 missing 0 3"), paste0(prefix, ".map"))
writeLines(c("f s1 0 0 1 -9 A A A AT 0 0", "f s2 0 0 1 -9 A G AT AT 0 0",
             "f s3 0 0 1 -9 G G A A 0 0", "f s4 0 0 1 -9 0 A 0 0 0 0"), paste0(prefix, ".ped"))
fixture <- ped2matrix(prefix)
near(fixture$Genotype, cbind(c(0L, 1L, 2L, NA_integer_), c(1L, 2L, 0L, NA_integer_), rep(NA_integer_, 4)))
stopifnot(identical(fixture$Genotype, codegeno(ped2matrix(prefix, "char")$Genotype)))

grDevices::png(file.path(root, "plots.png"), width = 2400, height = 1600, res = 160)
# Passing only a cached matrix makes an accidental recalculation impossible.
r2 <- LD(bed$Genotype)
cached <- HeatMap(LDmatrix = r2, SnpPosition = bed$Info[, 1:2], SnpName = bed$Info$ID)
inferred <- HeatMap(bed)
reader_cached <- HeatMap(bed, LDmatrix = r2)
near(cached$LDmatrix, inferred$LDmatrix)
near(reader_cached$LDmatrix, r2)
stopifnot(identical(inferred$SnpName, bed$Info$ID))
near(HeatMap(chars)$LDmatrix, LD(ped$Genotype))
fails(HeatMap(LDmatrix = matrix(0, 2, 3)))
fails(HeatMap(bed, LDmatrix = matrix(0, 2, 2)))
fails(HeatMap(LDmatrix = r2, SnpName = "wrong"))

# Circos receives the native character pair table without manual coercion.
info <- matrix(c("chr1", "9", "a", "chr1", "100", "b", "chrX", "20", "c"), ncol = 3, byrow = TRUE)
circos.ideogram(info)
stopifnot(circlize::get.cell.meta.data("xlim", sector.index = "chr1")[2] == 100)
links <- matrix(c("chr1", "9", "a", "chrX", "20", "c", "0.2"), nrow = 1)
circos.linksnp(links)
circos.ideogram(info)
circos.linksnp(matrix(character(), 0, 7))
fails(circos.linksnp(matrix("bad", 1, 7)))
grDevices::dev.off()
unlink(root, recursive = TRUE)
cat("PED numeric reading, encoding reuse, cached heatmaps and native circos inputs passed\n")
