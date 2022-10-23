#' linkage disequilibrium heatmap
#'
#' @import grid
#'
#' @param data SNP genotype data.frame or matrix
#' @param method measure linkage disequilibrium method
#' @param SnpPosition vector
#' @param SnpName vector
#' @param color heatmap corlor
#' @param showLDvalues logical
#' @param title heatmap's title
#' @param family Font for heatmap
#' @param ... NULL
#' @param cores thread number
#' 
#' @return HeatMap
#' @export
#' @examples NULL

HeatMap <- function(data, method="r^2", SnpPosition=NULL, SnpName=NULL, cores=1,
                    color=NULL, showLDvalues = FALSE, title="Pairwise LD", family=NULL, ...) {
  
  requireNamespace("grid", quietly = TRUE)
  if (method=="r^2") {
    LDmatrix <- LD(data, method="r^2", cores=cores)
    lower <- 0
    upper <- 1
  } else if (method=="D'") {
    LDmatrix <- LD(data, method="D'", cores=cores)
    lower <- 0
    upper <- 1
  } else if (method=="D") {
    LDmatrix <- LD(data, method="D", cores=cores)
    lower <- -1
    upper <- 1
  } else if (method=="RMI") {
    LDmatrix <- RMI(data, cores=cores)
    lower <- -0.1
    upper <- 1.6
  } else if (method=="MI") {
    LDmatrix <- MI(data, cores=cores)
    lower <- 0
    upper <- 1.6
  } else {
    stop("Invalid LD measurement")
  }
  #____________________________创建新的视图框_________________________________#
  grid::grid.newpage()
  mainvp <- viewport(width = unit(0.8, "snpc"), height = unit(0.8, "snpc"),
                     name="mainvp")
  flipVP <- viewport(width = unit(0.85, "snpc"), height = unit(0.85, "snpc"),
                     y=0.6, angle=-45, name="flipVP")
  
  #______________________________颜色选择______________________________________#
  if (is.null(color)) {
    color <- grey.colors(30)
  } else if (color=="blueTored") {
    colors <- colorRampPalette(rev(c("blue", "purple", "red")), space = "rgb")
    color <- colors(30)
  } else if (color=="YellowTored") {
    colors <- colorRampPalette(rev(c("white", "yellow", "orange", "red")), space = "rgb")
    color <- colors(30)
  } else {
    color <- color
  }
  #颜色选择
  mybreak <- -10:20/10
  color <- color[length(color):1]
  col_key <- color[cut(seq(lower, upper, 0.1), mybreak, include.lowest = T, right = F)]
  ###修改
  colcut <- as.character(cut(LDmatrix, mybreak, labels=as.character(color), right = T))
  
  #______________________________创建方形框____________________________________#
  rectcol <- ifelse(nrow(LDmatrix)<= 60, "white", NULL)
  
  ImageRect <- MakeImageRect(dim(LDmatrix)[1], dim(LDmatrix)[2], colcut,
                             name="heatmapRect", byrow=FALSE, rectlwd=1.1,
                             rectcol=rectcol)
  #___________________________设置绘图中的字体_________________________________#
  family <- ifelse(is.null(family), "serif", family)
  
  #________________________在方形框上是否显示LD的值____________________________#
  if(showLDvalues) {
    ImageText <- MakeImageText(dim(LDmatrix)[1], dim(LDmatrix)[2],
                               round(LDmatrix, digits = 2), name="heatmaptext",
                               flip = TRUE, family=family)
  } else {
    ImageText <- NULL
  }
  #____________________把方形框和LD值的图片放入flipVP视图中____________________#
  heatMap <- gTree(children=gList(ImageRect, ImageText), vp=flipVP, name="heatMap")
  #___________________________绘图中的主标题设置_______________________________#
  ##主标题放在heatmapVP
  main_title <- textGrob(title, 0.5, 1.05,
                         gp=gpar(cex=1.2, fontfamily=family),
                         name="main_title")
  
  ##标签标题放在heatmapVP
  if (is.null(SnpPosition)){
    mapLabel <- NULL
  } else {
    mapLabel <- paste("Physical Length: ", round(((max(SnpPosition)-min(SnpPosition))/1000),1), " Kb", sep="")
  }
  
  label_title <- textGrob(mapLabel, 0.5, 0.95,
                          gp=gpar(cex=1, fontfamily=family),
                          just="centre", name= "label_title")
  title <- gTree(children=gList(main_title, label_title), name = "title")
  #___________________________绘图中SNP的标签设置______________________________#
  label <- HeatmapLabels(LDmatrix, SnpPosition, SnpName, vp=flipVP, family=family)
  #________________________________设置图例____________________________________#
  legend <- HeatmapLegend(col_key, method, family=family)
  
  LDheatmap <- gTree(children=gList(heatMap, legend, label, title),
                     vp=mainvp, name= "LDheatmap")
  
  grid.draw(LDheatmap)
  
  invisible(structure(list(LDmatrix=LDmatrix, LDheatMap=LDheatmap,
                           SnpName=SnpName, SnpPosition=SnpPosition,
                           Color=color), class = "HeatMap"))
}

##绘图函数的类为HeatMap
#' @importFrom grid grid.draw
#' @export
grid.draw.HeatMap <- function(x, recording = TRUE) {
   grid.draw(x$LDheatMap)
}

#' @export
print.HeatMap <- function(x, ...) {
   grid.draw(x)
}