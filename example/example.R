# install GWLD
install.packages("devtools")

devtools::install_github("Rong-Zh/GWLD/GWLD-R")

# loading package
library(GWLD)

# loading example data with base-type genotype
# data("duck")
# data <- duck$SNP
# SNP <- data$genotype
# Info <- data$info
# geno012 <- codegeno(data$genotype, sep = "/")

## VCF format
# 1、read vcf
SNP <- read.vcf(filename = "example-3.vcf", genotype = "int")
SNP <- data.frame(SNP)
Info <- SNP[, 1:5]
SNP <- data.frame(t(SNP[, 10:ncol(SNP)]))
geno012 <- apply(SNP, 2, as.numeric)
row.names(geno012) <- row.names(SNP)
colnames(geno012) <- Info$ID

# 2 read vcf
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
SNP <- read.plink("example")
Info <- SNP$bim
geno012 <- SNP$bed

# Using the plink data above as an example
# Calculate LD,
r2 <- LD(geno012, cores = 2)
mi <- LD(geno012, cores = 2)
RMI <- RMI(geno012, cores = 2)

# Calculate the distance between different chromosomes
rmi_circos <- calc_circos(geno012, Info, method = "RMI", threshold = 0.1, cores = 1)

# For more details, see https://doi.org/10.1093/g3journal/jkad154 Box 1




