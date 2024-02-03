# install GWLD
# Install Rtools first for Windows. Visit https://cran.r-project.org/bin/.
# install LLVM and libomp first for macOS M chip.
install.packages("devtools")
devtools::install_github("Rong-Zh/GWLD/GWLD-R")

# loading package
library(GWLD)

# loading example data with base-type genotype
data("duck")
dta <- duck$SNP
SNP <- dta$genotype
Info <- dta$info
geno012 <- codegeno(dta$genotype, sep = "/")

## VCF format
# Approach 1
SNP <- read.vcf(filename = "example-3.vcf", genotype = "int")
SNP <- data.frame(SNP)
Info <- SNP[, 1:3]
SNP <- data.frame(t(SNP[, 10:ncol(SNP)]))
geno012 <- apply(SNP, 2, as.numeric)
row.names(geno012) <- row.names(SNP)
colnames(geno012) <- Info$ID

# Approach 2
SNP <- vcf2matrix("example-3.vcf")
Info <- SNP$Info
head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

## Plink format
# read plink format(ped, map)
SNP <- ped2matrix(FilePrefix = "example-2")
Info <- SNP$Info
head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

# read Plink format data(bed, bim, fam)
# Approach 1
SNP <- read.plink("example")
Info <- SNP$bim
Info <- Info[, c(1, 4, 2)]
colnames(Info) <- c("CHROM", "POS", "ID")
geno012 <- SNP$bed

# Approach 2
SNP <- bed2matrix("example")
Info <- SNP$Info
geno012 <- SNP$Genotype

# Using the plink data above as an example
# Calculating values with different methods.
r2 <- LD(geno012, cores = 1)
mi <- MI(geno012, cores = 1)
RMI <- RMI(geno012, cores = 1)
value <- GWLD(geno012, method = "RMI", cores = 1)

# Calculating the value between two loci.
r2 <- LD(geno012[, 1], geno012[, 2], method = "r^2", cores = 1)
value <- GWLD(geno012[, 1], geno012[, 2], method = "RMI", cores = 1)

# Calculate the distance between different chromosomes
rmi_circos <- calc_circos(geno012, Info, method = "RMI", threshold = 0.1, cores = 1)

# For more details, see https://doi.org/10.1093/g3journal/jkad154 Box 1
