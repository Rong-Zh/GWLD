library(GWLD)
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,check.attributes=FALSE)))
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
reference <- function(g,sep) {
  parsed <- lapply(g,function(x) {
    if (is.na(x) || x %in% c("", ".", "NA")) return(c(NA_character_,NA_character_))
    p <- strsplit(x,if (sep=="") "" else "[/|]",fixed=sep=="")[[1]]
    stopifnot(length(p)==2L)
    p[p %in% c(".","NA","")] <- NA_character_
    p
  })
  called <- table(unlist(parsed))
  if (length(called)) called <- called[order(-as.numeric(called),names(called),method="radix")]
  major <- if(length(called)) names(called)[1] else NA_character_
  code <- vapply(parsed,function(p) if(anyNA(p)) NA_integer_ else as.integer(2L-sum(p==major)),integer(1))
  list(code=code,alleles=as.character(names(called)),counts=as.numeric(called),missing=sum(is.na(code)))
}
check <- function(g,sep) {
  ref <- reference(g,sep)
  near(recode(g,sep),ref$code)
  stats <- stats_freq(g,sep)
  stopifnot(identical(stats$allele,ref$alleles))
  near(stats$`allele counts`,ref$counts)
  near(stats$`allele freq`,ref$counts/sum(ref$counts))
  near(stats$`missing counts`,ref$missing)
  near(stats$`missing rate`,if(length(g)) ref$missing/length(g) else NA_real_)
}
set.seed(2026)
for (pool in list(c("A/A","A/G","G/A","G/G","A|G","G|A",NA,"./.","A/.",""),
                  c("A/AT","AT/AT","A/A","AT|A",NA,"NA/AT","."),
                  c("C/C","A/G","C/A","G/C","G/G",NA))) {
  for (trial in 1:20) check(sample(pool,200,replace=TRUE),"/|")
}
check(c("AA","AG","GA","GG",NA,"..","A."),"")
check(rep(NA_character_,8),"/|")
check(character(),"/")
check(c("AT/AT","AT/AT",NA),"/")
near(recode(c("A/AT","AT/AT","A/A"),"/"),c(1L,2L,0L))
near(recode(c("G/G","A/A"),"/"),c(2L,0L)) # lexical tie-break
g <- matrix(sample(c("A/A","A/G","G/G",NA,"./."),6000,TRUE),nrow=200)
dimnames(g) <- list(paste0("s",1:nrow(g)),paste0("v",1:ncol(g)))
result <- codegeno(g,"/")
expected <- vapply(seq_len(ncol(g)),function(j) reference(g[,j],"/")$code,integer(nrow(g)))
near(result,expected)
stopifnot(identical(dimnames(result),dimnames(g)))
near(codegeno(as.data.frame(g),"/"),expected)
stopifnot(identical(dim(codegeno(matrix(character(),0,3),"/")),c(0L,3L)),
          identical(dim(codegeno(matrix(character(),3,0),"/")),c(3L,0L)))
near(codegeno(matrix(NA_character_,5,2),"/"),matrix(NA_integer_,5,2))
for (bad in c("A","A/G/T","A//","/A","A/","AAA")) fails(recode(bad,"/"))
fails(recode("AAA",""))
cat("Genotype encoding, allele-frequency, missing-value and empty-matrix tests passed\n")

# Force hash-table growth and then reuse it across unrelated matrix columns.
many <- c(paste0("allele", 1:96, "/allele", 96:1), NA_character_, "allele1/.")
check(many, "/|")
expanded <- matrix(rep(many, 3), ncol = 3)
expanded[, 2] <- rep(c("A/A", "A/G", "G/G", NA_character_), length.out = nrow(expanded))
expected <- vapply(seq_len(ncol(expanded)), function(j) reference(expanded[, j], "/|")$code, integer(nrow(expanded)))
near(GWLD:::code_Mat(expanded, "/|"), expected)
cat("Many-allele hash growth and matrix-column reuse passed\n")
