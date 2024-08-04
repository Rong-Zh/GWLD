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

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

## VCF format
# Approach 1
vcf <- read.vcf(filename = "example-3.vcf", genotype = "int")
Info <- data.frame(CHROM = vcf[, 1], POS = as.numeric(vcf[, 2]), ID = vcf[, 3])
SNP <- data.frame(t(vcf[, 10:ncol(vcf)]))
geno012 <- apply(SNP, 2, as.numeric)
row.names(geno012) <- row.names(SNP)
colnames(geno012) <- Info$ID

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

# Approach 2
SNP <- vcf2matrix("example-3.vcf")
Info <- SNP$Info

head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

## Plink format
# read plink format(ped, map)
SNP <- ped2matrix(FilePrefix = "example-2")
Info <- SNP$Info
head(SNP$Genotype)
geno012 <- codegeno(SNP$Genotype, sep = "/")

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

# read Plink format data(bed, bim, fam)
# Approach 1
SNP <- read.plink("example")
bim <- SNP$bim
Info <- data.frame(CHROM = bim[, 1], POS = as.numeric(bim[, 4]), ID = bim[, 2])
geno012 <- SNP$bed

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

# Approach 2
SNP <- bed2matrix("example")
Info <- SNP$Info
geno012 <- SNP$Genotype

HeatMap(geno012, method = "r^2", SnpPosition = Info[, 1:2], SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)

# Using the plink data above as an example
# Calculating values with different methods.
r2 <- LD(geno012, cores = 1)
mi <- MI(geno012, cores = 1)
RMI <- RMI(geno012, cores = 1)
value <- GWLD(geno012, method = "RMI", cores = 1)

# Calculating the value between two loci.
r2 <- LD(geno012[, 1], geno012[, 2], method = "r^2", cores = 1)
value <- GWLD(geno012[, 1], geno012[, 2], method = "RMI", cores = 1)

# Calculate the decay from result
RMI <- RMI(geno012, cores = 1)
rmi_decay <- decay(RMI, Info)

# or calculate circos result from data
rmi_decay <- calc_decay(geno012, Info, method = "RMI")

# Calculate the decay from result
RMI <- RMI(geno012, cores = 1)
rmi_circos <- circos(RMI, Info, threshold = 0.1)

# Calculate the distance between different chromosomes
rmi_circos <- calc_circos(geno012, Info, method = "RMI", threshold = 0.1, cores = 1)

#_________________________________Plot decay____________________________________
library(ggplot2)
duck_decay <- duck$Decay |> dplyr::filter(method == "RMI")

ggplot(data = duck_decay, aes(x = dist / 1000, y = mean)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "Physical distance/ Kb", y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 15, hjust = 0),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = c(.90, .85),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 300), expand = c(0, 5), breaks = seq(0, 300, 20))

# The raw result calculated by the calc_decay function
tbl <- readRDS(file = "decay.RDS")
ld_avg <- function(ld, bin) {
  ld_out <- c()
  for (i in 1:(length(dist) - 1)) {
    index <- (ld[, 1] < dist[i + 1]) & (ld[, 1] > dist[i])
    ld_out[i] <- mean(ld[index, 2])
  }
  ld_decay <- cbind(dist = dist[-1], mean = ld_out)
  return(ld_decay)
}

# Calculate the average of the bins
bin <- c(seq(0, 1000, 100), seq(1000, 10000, 1000), seq(10000, 300000, 1000))
duck_decay <- data.frame(na.omit(ld_avg(tbl[, c(3, 4)], bin)))

# plot
ggplot(data = duck_decay, aes(x = dist / 1000, y = mean)) +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  labs(x = "Physical distance/ Kb", y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 15, hjust = 0),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = c(.90, .85),
    legend.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 300), expand = c(0, 5), breaks = seq(0, 300, 20))

#_________________________________plot Ciros____________________________________
circosdata <- duck$Circos
circos.ideogram(circosdata$chr)
circos.linksnp(circosdata$linkdata)

# plot using the example data,
# It is possible to construct all or partial chromosome length data
chr_info <- Info |>
  data.frame() |>
  dplyr::mutate_at(c("CHROM", "POS"), .funs = as.numeric)

rmi_res <- rmi_circos |>
  data.frame() |>
  dplyr::mutate_at(c("CHR_1", "POS_1", "CHR_2", "POS_2", "Value"), .funs = as.numeric)

# show all chromosomes of the duck data
circos.ideogram(circosdata$chr)
circos.linksnp(rmi_res)

# show partial chromosomes of the duck data
circos.ideogram(chr_info)
circos.linksnp(rmi_res)

# Construct chromosome length information
chr_info <- structure(
  list(
    chr = c(1, 2, 3, 4, 5, 15),
    start = c(0, 0, 0, 0, 0, 0),
    end = c(164695780, 124299693, 49759293, 10956466, 60265033, 7563935)
  ),
  class = c("list", "custom")
)

circos.ideogram(chr_info)
circos.linksnp(rmi_res)

# For more details, see https://doi.org/10.1093/g3journal/jkad154 Box 1
