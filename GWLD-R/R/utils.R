#' @import grid
MakeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE, rectlwd=NULL, rectcol=NULL) {
  xx <- (1:ncol)/ncol
  yy <- (1:nrow)/nrow
  # 创建坐标对，重复列号(如果 byrow = TRUE)或行号(如果 byrow = FALSE)以强制该类型的填充
  if(byrow) {
      right <- rep(xx, nrow)
      top <- rep(yy, each=ncol)
  } else {
      right <- rep(xx, each=nrow)
      top <- rep(yy, ncol)
  }
  if (is.null(rectcol)) {
      rectcol = NA ###设置矩形框的颜色为NA
  }
  rectGrob(x=right, y=top,
           width=1/ncol, height=1/nrow,
           just=c("right", "top"),
           gp=gpar(lwd=rectlwd, col=rectcol, fill=cols),
           name=name)
}

MakeImageText <- function(nrow, ncol, cols, name, flip = FALSE, family=family) {
  cols <- as.character(cols)
  cols[is.na(cols)] <- ""
  xx <- (1:ncol)/ncol
  yy <- (1:nrow)/nrow

  # 需要以不同的顺序填充单元格，就像生成图像一样
  if(flip) {
      right <- rep(xx, each = nrow)
      top <- rep(yy, ncol)
  } else {
      right <- rep(xx, nrow)
      top <- rep(yy, each=ncol)
  }
  textGrob(cols, x=right-0.5/nrow, y=top-0.5/ncol, rot= 45,
           gp=gpar(cex=0.45, fontfamily=family),
           name=name)
}

HeatmapLegend <- function(x, y, legend.colors, legend.title, breaks, n.breaks, 
                          labels, legend.family) {
  #归一化到0,1
  labels_pos <- seq(0, 1, length.out=n.breaks)
  legend_breaks <- seq(breaks[1], breaks[2], length.out=n.breaks+1)
  pos <- labels_pos[cut(labels, legend_breaks, labels=FALSE, right=TRUE, include.lowest=T)]
  right  <- rep((1:n.breaks)/n.breaks, each=1)
  top<- rep(1, n.breaks)
  legend_title <- textGrob(legend.title, x=0.5, y=1.3, name="legend_title", just = "bottom", 
                           gp=gpar(cex=0.85, fontfamily=legend.family))
  legend_rect <- rectGrob(x=right, y=top, width=1/n.breaks, height=1,
                          just=c("right", "top"), name="legend_rect", gp = gpar(fill=legend.colors, col=NA))
  legend_box <- linesGrob(x=c(0,1,1,0,0), y=c(0,0,1,1,0), name = "legend_box", gp=gpar(col="black", lwd=0.4))
  legend_ticks <- segmentsGrob(x0=pos, 
                               y0=rep(0,length(pos)), 
                               x1=pos, 
                               y1=rep(-0.15,length(pos)),
                               gp=gpar(col="black", lwd=0.7),
                               name = "legend_ticks")
  legend_text <- textGrob(labels, x=pos, y=rep(-0.2, length(pos)), just = "top",
                          gp=gpar(cex=0.75, fontfamily=legend.family),
                          name="legend_text")
  legend_vp <- viewport(x=x, y=y, width=0.4, height=0.035, name = "legend_vp")
  legend <- gTree(children = gList(legend_title, legend_rect, legend_ticks, legend_box, legend_text), vp=legend_vp, name="legend")
  legend
}

HeatmapLabels <- function(LDmatrix, pos, name, vp, label.size=label.size, family=family) {
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  snp <- ((1:nsnps-1) + 0.5) / nsnps
  geneMapLocation <- ifelse(is.null(pos), -0.004, -0.15)
  seq.x <- c(0.5*geneMapLocation + 1/(nsnps*2),
             1+0.5*geneMapLocation - 1/(nsnps*2))
  seq.y <- c(-0.5*geneMapLocation + 1/(nsnps*2),
             1-0.5*geneMapLocation - 1/(nsnps*2)) 
  ##设置距离
  if(is.null(pos)) {
    pos <- 1:nsnps ##为数值
    total.dist <- max(pos)-min(pos)
    min.dist <- min(pos)
    regionx <- seq.x[1] + ((pos-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
    regiony <- seq.y[1] + ((pos-min.dist)/total.dist)*(seq.y[2]-seq.y[1])
    diagonal <- NULL
    segments <- NULL
  } else {
    ###直线段
    diagonal <- linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal", vp=vp)
    ###斜线段
    total.dist <- max(pos)-min(pos)
    min.dist <- min(pos)
    regionx <- seq.x[1] + ((pos-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
    regiony <- seq.y[1] + ((pos-min.dist)/total.dist)*(seq.y[2]-seq.y[1])
    segments <- segmentsGrob(snp, snp, regionx, regiony, name="segments", vp=vp)
  }
  ##设置SNP名
  if (!is.null(segments) & is.null(name)) {
    ###点的位置和符号,符号大小得微调
    symbols_pch <- pointsGrob(snp[1:nsnps], snp[1:nsnps], pch="*",
                              gp=gpar(cex=-0.02*nsnps+2, bg="blue", col="blue"),
                              name="symbols_pch", vp=vp)
    SNPnames <- NULL
  } else if(!is.null(name)) {
    ##宽度转换
    name_gap <- convertWidth(grobWidth(textGrob("r")), "npc", valueOnly=TRUE)/sqrt(2)
    ind <- match(name, row.names(LDmatrix), nomatch=0)
    #点的位置和符号
    symbols_pch <- pointsGrob(snp[ind], snp[ind], pch="*",
                              gp=gpar(cex=-0.02*nsnps+2, bg="blue", col="blue"),
                              name="symbols_pch", vp=vp)
    #距离检测(避免SNP字体重叠并更新ind, 让重叠的名字不显示)
    min_dist <- convertHeight(grobHeight(textGrob("rs", gp=gpar(fontsize=label.size))), "npc", valueOnly=TRUE)*1.1
    if(length(ind)>=2) {
      start <- ind[1]
      res <- c(start)
      for(i in 2:length(ind)) {
        if(abs(regionx[ind[i]] - regionx[start]) >= min_dist) {
          res <- append(res, ind[i])
          start <- ind[i]
        }     
      }
      ind <- res
    }
    ## 如果绘图窗口太小,字体会出界
    SNPnames <- textGrob(name[ind], just=c("left", "bottom"), rot=90,
                         regionx[ind]-sqrt(2+0.5)*name_gap,
                         regiony[ind]+sqrt(2+0.5)*name_gap,
                         gp=gpar(fontsize=label.size, col="blue", fontfamily=family),
                         name="SNPnames", vp=vp)
  } else {
    symbols_pch <- NULL
    SNPnames <- NULL
  }
  label <- gTree(children=gList(symbols_pch, segments, diagonal, SNPnames), name ="label")
  label
}


