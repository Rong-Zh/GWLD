library(GWLD)
check_savefig <- function() {
  root <- tempfile("gwld-savefig-")
  dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  fail <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
  image_size <- function(path) {
    file <- file(path,"rb")
    on.exit(close(file))
    readBin(file,"raw",16L)
    readBin(file,"integer",2L,size=4L,endian="big")
  }
  picture <- grid::rectGrob()
  path <- file.path(root,"plot.png")
  savefig(path,picture,width=1,height=0.5,dpi=80,scale=2)
  stopifnot(identical(image_size(path),c(160L,80L)))
  savefig(path,picture,width=25.4,height=12.7,units="mm",dpi=80)
  stopifnot(identical(image_size(path),c(80L,40L)))
  savefig(path,picture,width=80,height=40,units="px",dpi=80)
  stopifnot(identical(image_size(path),c(80L,40L)))

  savefig(path,picture,width=1,height=0.5,dpi=80,device=grDevices::png)
  stopifnot(identical(image_size(path),c(80L,40L)))
  pdf <- file.path(root,"dimensions.pdf")
  savefig(pdf,picture,width=2.54,height=5.08,units="cm",compress=FALSE)
  stopifnot(any(grepl("/MediaBox [0 0 72 144]",readLines(pdf,warn=FALSE),fixed=TRUE,useBytes=TRUE)))

  # Custom devices receive inches, just like PDF; preserve the prior device.
  grDevices::pdf(file.path(root,"prior.pdf"))
  prior <- grDevices::dev.cur()
  on.exit(if (prior %in% grDevices::dev.list()) grDevices::dev.off(prior),add=TRUE)
  seen <- NULL
  custom <- function(filename,width,height,...) {
    seen <<- c(width,height)
    grDevices::pdf(file=filename,width=width,height=height,...)
  }
  savefig(file.path(root,"custom.pdf"),picture,width=2.54,height=5.08,units="cm",scale=2,device=custom)
  stopifnot(isTRUE(all.equal(seen,c(2,4))),identical(grDevices::dev.cur(),prior))
  # A drawing failure must close the output device and restore the prior one.
  devices <- grDevices::dev.list()
  fail(savefig(file.path(root,"failed.pdf"),new.env(),device=custom))
  stopifnot(identical(grDevices::dev.list(),devices),identical(grDevices::dev.cur(),prior))
  for (value in list(0,-1,NA_real_,Inf,c(1,2))) {
    fail(savefig(path,picture,width=value))
    fail(savefig(path,picture,scale=value))
    fail(savefig(path,picture,dpi=value))
  }
  fail(savefig(path,picture,units="invalid"))
}
check_savefig()
cat("Figure scaling, units, custom devices and cleanup passed\n")
