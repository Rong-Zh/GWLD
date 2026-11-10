library(GWLD)
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,check.attributes=FALSE)))
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
root <- tempfile("gwld-vcf-")
dir.create(root)
path <- file.path(root,"input.vcf")
header <- paste0("##fileformat=VCFv4.2\n",
  "##contig=<ID=1,length=1000000>\n",
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts_1\ts_2\ts_3\ts_4\n")
text <- paste0(header,
  "1\t2000\ta\tA\tG\t.\t.\t.\tDP:GT\t12:0/0\t10:1/0\t.:0|.\t.\n",
  "1\t1000\tb\tC\tT\t.\t.\t.\tGT\t1|1\t0/0\t0/1\t1/.")
writeChar(text,path,eos=NULL,useBytes=TRUE)
expected <- matrix(c(0L,1L,NA_integer_,NA_integer_,2L,0L,1L,NA_integer_),nrow=4)
data <- vcf2matrix(path)
near(data$Genotype,expected)
stopifnot(is.integer(data$Genotype),identical(rownames(data$Genotype),paste0("s_",1:4)),
          identical(colnames(data$Genotype),c("a","b")))
near(data$Info$POS,c(2000,1000))
gzip <- paste0(path,".gz")
connection <- gzfile(gzip,"wb")
writeChar(text,connection,eos=NULL,useBytes=TRUE)
close(connection)
near(vcf2matrix(gzip,threads=2L)$Genotype,expected)
near(GWLD(vcf2matrix(gzip)$Genotype,method="MI"),GWLD(expected,method="MI"))
table <- read.vcf(path,genotype="int")
near(matrix(as.integer(table[,10:13]),nrow=2),t(expected))
stopifnot(identical(unname(read.vcf(path,"char")[1,11]),"G/A"),
          identical(unname(read.vcf(path,"None")[1,11]),"10:1/0"))
alleles <- read.vcf(path,"allele")
stopifnot(identical(dim(alleles),c(2L,17L)),
          identical(colnames(alleles)[10:11],c("s_1_A","s_1_a")),
          identical(alleles[1,12:13],c(s_2_A="1",s_2_a="0")))
near(read.vcf(gzip,"int"),table)
fails(vcf2matrix(path,threads=0L))
fails(read.vcf(path,"bad-mode"))
converter <- Sys.getenv("GWLD_TEST_VCF_CONVERTER")
if (nzchar(converter)) {
  for (mode in c("wz","wb")) {
    converted <- file.path(root,if (mode=="wz") "bgzf.vcf.gz" else "input.bcf")
    stopifnot(system2(converter,c(shQuote(path),shQuote(converted),mode))==0L)
    near(vcf2matrix(converted,threads=2L)$Genotype,expected)
  }
}
cli <- Sys.getenv("GWLD_TEST_CLI")
if (nzchar(cli)) {
  prefix <- file.path(root,"cli")
  stopifnot(system2(cli,c("convert","--vcf",shQuote(gzip),"--format","genotype","-o",shQuote(prefix)))==0L)
  codes <- as.matrix(read.table(paste0(prefix,".genotype.tsv"),header=TRUE,comment.char="",check.names=FALSE)[,6:9])
  near(t(codes),expected)
  stopifnot(system2(cli,c("calculate","--vcf",shQuote(gzip),"--method","mi","--heatmap","--decay","-o",shQuote(prefix)))==0L)
  near(as.matrix(read.table(paste0(prefix,".heatmap"))),GWLD(expected,method="MI"))
  near(as.matrix(read.table(paste0(prefix,".decay"),colClasses="character")),calc_decay(expected,data$Info,method="MI"))
}
# POS remains exact when VCF coordinates exceed the precision of an R double.
large_path <- file.path(root,"large-position.vcf")
writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
             "1\t9007199254740993\ta\tA\tG\t.\t.\t.\tGT\t0/1"),large_path)
stopifnot(identical(vcf2matrix(large_path)$Info$POS,"9007199254740993"))
writeChar(header,path,eos=NULL,useBytes=TRUE)
stopifnot(identical(dim(read.vcf(path)),c(0L,13L)))
fails(vcf2matrix(path))
unlink(root,recursive=TRUE)
cat("Shared VCF reader, gzip/BGZF/BCF and CLI/R consistency tests passed\n")
