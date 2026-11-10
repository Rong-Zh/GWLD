library(GWLD)
near <- function(a,b,tolerance=1e-10) stopifnot(isTRUE(all.equal(a,b,tolerance=tolerance,check.attributes=FALSE)))
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
cores <- if (Sys.getenv("GWLD_TEST_OPENMP","ON")=="OFF") 1L else 2L

# Genotype validation must happen before any Rcpp integer conversion.
a <- c(0,1,2,0,1,2)
for (bad in c(0.9,1.9,2.9,Inf,-Inf,3,-2)) {
  invalid <- a; invalid[1] <- bad
  for (method in c("MI","RMI","r2","D","Dprime")) fails(GWLD(invalid,a,method))
  fails(MI(invalid,a)); fails(RMI(a,invalid)); fails(LD(invalid,a))
  fails(GWLD(matrix(invalid,3,2)))
  fails(calc_decay(matrix(invalid,3,2),cbind("1",1:2,c("a","b"))))
}
near(GWLD(a,a,"r2"),1)
for (invalid in c(0,1.5,Inf,NA_real_)) fails(GWLD(matrix(a,3,2),cores=invalid))

# A heterozygote-heavy table has two likelihood modes. The better one is negative.
table <- matrix(c(0L,2L,0L,2L,262L,5L,2L,5L,2L),3,3,byrow=TRUE)
x <- y <- integer()
for (i in 1:3) for (j in 1:3) {
  x <- c(x,rep(i-1L,table[i,j])); y <- c(y,rep(j-1L,table[i,j]))
}
near(LD(x,y,method="D"),-0.23361366206657,1e-8)
near(LD(2L-x,y,method="D"),-LD(x,y,method="D"),1e-10)
near(LD(x,y),LD(2L-x,y),1e-10)

# Exercise several parallel batches and lossless in-memory/file output.
set.seed(2026)
g <- matrix(sample(c(0L,1L,2L,NA_integer_),256*80,TRUE),256,80)
colnames(g) <- paste0("v",1:80)
info <- data.frame(CHROM=rep(c("chr1","chrX"),40),POS=seq_len(80)*10,ID=colnames(g))
root <- tempfile("gwld-regressions-"); dir.create(root)
for (method in c("MI","RMI","r2","D","Dprime")) {
  m <- GWLD(g,method=method,cores=cores)
  stopifnot(identical(dimnames(m),list(colnames(g),colnames(g))))
  d <- calc_decay(g,info,method=method,cores=cores)
  c <- calc_circos(g,info,method=method,cores=cores,threshold=-1)
  near(d,decay(m,info)); near(c,circos(m,info,threshold=-1))
  for (kind in c("decay","circos")) for (from_matrix in c(FALSE,TRUE)) {
    path <- file.path(root,paste(method,kind,from_matrix,"tsv",sep="."))
    result <- if (kind=="decay") {
      if (from_matrix) decay(m,info,file=path) else calc_decay(g,info,method,cores,file=path)
    } else {
      if (from_matrix) circos(m,info,threshold=-1,file=path) else calc_circos(g,info,method,cores,threshold=-1,file=path)
    }
    stopifnot(identical(result,path))
    streamed <- as.matrix(read.table(path,header=TRUE,sep="\t",colClasses="character",na.strings=character(),check.names=FALSE))
    near(streamed,if (kind=="decay") d else c)
  }
}
tiny <- matrix(NA_real_,2,2); tiny[1,2] <- 1e-8
tiny_info <- cbind(c("chr1","chrX"),1:2,c("a","b"))
stopifnot(as.numeric(circos(tiny,tiny_info,threshold=1e-9)[1,7])==1e-8)
empty <- file.path(root,"empty.tsv")
calc_circos(g,info,threshold=100,file=empty)
stopifnot(length(readLines(empty))==1L)
fails(calc_decay(g,info,file=""))
fails(calc_decay(g,info,maxdist=0.5))

# Cached matrices are never silently relabeled. Direct plotting supports palettes.
grDevices::png(file.path(root,"plots.png"),width=2400,height=1600,res=160)
input <- list(Genotype=g,Info=info)
m <- LD(g)
permuted <- m[c(2,1,3:80),c(2,1,3:80)]
fails(HeatMap(input,LDmatrix=permuted))
fails(HeatMap(LDmatrix=m,SnpName=rev(colnames(g))))
near(HeatMap(input,LDmatrix=m,color=c("yellow","orange","red"))$LDmatrix,m)
HeatMap(LDmatrix=matrix(NA_real_,1,1),SnpName="single",SnpPosition=100)
unnamed <- g; colnames(unnamed) <- rep(".",ncol(unnamed))
HeatMap(unnamed,color="YellowTored")
grDevices::dev.off()

# Family/individual identity is preserved, including separator escaping.
prefix <- file.path(root,"samples")
writeLines("1 a 0 10",paste0(prefix,".map"))
writeLines(c("f1 sample 0 0 1 -9 A A","f2 sample 0 0 1 -9 G G",
             "f:3 sample% 0 0 1 -9 A G"),paste0(prefix,".ped"))
stopifnot(identical(rownames(ped2matrix(prefix)$Genotype),c("f1:sample","f2:sample","f%3A3:sample%25")))
writeLines(c("f1 sample 0 0 1 -9 A A","f1 sample 0 0 1 -9 G G"),paste0(prefix,".ped"))
fails(ped2matrix(prefix))
# Numeric-looking IDs and alleles are identifiers, so leading zeros matter.
writeLines("01 001 0 10",paste0(prefix,".map"))
writeLines(c("01 001 0 0 1 -9 01 01","1 1 0 0 1 -9 1 1"),paste0(prefix,".ped"))
ped <- ped2matrix(prefix)
stopifnot(identical(rownames(ped$Genotype),c("01:001","1:1")),
          identical(colnames(ped$Genotype),"001"), identical(unname(ped$Genotype[,1]),c(0L,2L)))
writeLines("01 001 0 10 G A",paste0(prefix,".bim"))
writeLines(c("01 001 0 0 1 -9","1 1 0 0 1 -9"),paste0(prefix,".fam"))
writeBin(as.raw(c(0x6c,0x1b,0x01,0x03)),paste0(prefix,".bed"))
bed <- bed2matrix(prefix)
stopifnot(identical(rownames(bed$Genotype),c("01:001","1:1")),identical(colnames(bed$Genotype),"001"))
# Reader metadata must not lose coordinate bits before native validation.
writeLines("01 001 0 9007199254740993",paste0(prefix,".map"))
stopifnot(identical(ped2matrix(prefix)$Info$POS,"9007199254740993"))
writeLines("01 001 0 4.611686018427387903e18 G A",paste0(prefix,".bim"))
stopifnot(identical(bed2matrix(prefix)$Info$POS,"4611686018427387903"))
writeLines("01 001 0 1.00000000000000001",paste0(prefix,".map"))
fails(ped2matrix(prefix))
writeLines("01 001 0 1.00000000000000001 G A",paste0(prefix,".bim"))
fails(bed2matrix(prefix))
if (cores==1L) fails(GWLD(g,cores=2L))
unlink(root,recursive=TRUE)
cat("Multimodal LD, exact genotype validation, pair streaming/precision, matrix IDs and sample IDs passed\n")
