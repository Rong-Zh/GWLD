#' @import grid
MakeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE,rectlwd=NULL,rectcol=NULL) {
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

MakeImageText <- function(nrow, ncol, cols, name, flip = FALSE, family=NULL) {
  cols <- as.character(cols)
  cols[is.na(cols)] <- ""
  cols <- paste("   ", cols)
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

  #设置字体
  family <- ifelse(is.null(family), "serif", family)

  textGrob(cols, x=right-0.035, y=top-0.015, rot= 45,
           gp=gpar(cex=0.45, fontfamliy=family),
           just=c("right", "top"),
           name=name)
}

##legend 需要修改
HeatmapLegend <- function(color, LDmeasure, family=NULL) {
  legend_ImageRect <- MakeImageRect(2, length(color), cols=c(rep(NA,length(color)), color[1:length(color)]),
                                    "legend_colorkey")
  legendvp <- viewport(x=1.1, y=-0.10, height=0.10, width=0.5, just=c("right","bottom"), name="legend")
  #添加标签 'Color key'
  if(LDmeasure=="r^2") {
    ttt <- expression(paste(r^2," Color Key"))
    #在颜色调上添加标签
    legend_labels <- textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="legend_labels")
    #在颜色调底部轴上绘制刻度
    legend_ticks <- segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="legend_ticks")
  } else if(LDmeasure=="D'") {
    ttt <- "D' Color Key"
    #在颜色调上添加标签
    legend_labels <- textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="legend_labels")
    #在颜色调底部轴上绘制刻度
    legend_ticks <- segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6), name="legend_ticks")
  } else if(LDmeasure=="D") {
    ttt <- "D Color Key"
    #在颜色调上添加标签
    legend_labels <- textGrob(paste(seq(-1,1,0.5)), x=0.25*0:4, y=0.25, gp=gpar(cex=0.6), name="legend_labels")
    #在颜色调底部轴上绘制刻度
    legend_ticks <- segmentsGrob(x0=c(0:4)*0.25, y0=rep(0.4,6), x1=c(0:4)*0.25, y1=rep(0.5,6), name="legend_ticks")
  } else if(LDmeasure=="RMI") {
    ttt <- "RMI Color Key"
    #在颜色调上添加标签
    legend_labels <- textGrob(paste(c(-0.1, 0.24, 0.58, 0.92, 1.26, 1.6)), x=0.2*0:5, y=0.25, gp=gpar(cex=0.6), name="legend_labels")
    #在颜色调底部轴上绘制刻度
    legend_ticks <- segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="legend_ticks")
  } else if(LDmeasure=="MI") {
    ttt <- "MI Color Key"
    #在颜色调上添加标签
    legend_labels <- textGrob(paste(0.32*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="legend_labels")
    #在颜色调底部轴上绘制刻度
    legend_ticks <- segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="legend_ticks")
  }
  # 设置字体
  family <- ifelse(is.null(family), "serif", family)
  legend_title <- textGrob(ttt, x=0.5, y=1.25, name="legend_title", gp=gpar(cex=0.8, fontfamliy=family))
  #在颜色条周围加一个框
  legend_box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="legend_box")
  legend <- gTree(children=gList(legend_ImageRect, legend_title, legend_labels, legend_ticks, legend_box), name = "legend", vp=legendvp)
  legend
}

HeatmapLabels <- function(LDmatrix, pos, name, vp, family=NULL) {
  
  ###设置字体
  family <- ifelse(is.null(family), "serif", family)
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
    pos <- 1:nsnps ##有数值
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
                              gp=gpar(cex=-0.01*nsnps+1.6, bg="blue", col="blue"),
                              name="symbols_pch", vp=vp)
    SNPnames <- NULL
  } else if(!is.null(name)) {
    name_gap <- convertWidth(grobWidth(textGrob("rs")), "npc", valueOnly=TRUE)/sqrt(2)##宽度转换
    ind <- match(name, row.names(LDmatrix), nomatch=0)##整数
    #点的位置和符号
    symbols_pch <- pointsGrob(snp[ind], snp[ind], pch="*",
                          gp=gpar(cex=-0.01*nsnps+1.6, bg="blue", col="blue"),
                          name="symbols_pch", vp=vp)
    #距离检测,避免SNP字体重叠并更新ind, 让重叠的名字不显示
    dist <- (pos-min.dist)/total.dist
    if(length(ind)>=2) {
      start = ind[1]
      res <- c(start)
      for(i in 2:length(ind)) {
        if(abs(dist[ind[i]] - dist[start]) >= 0.015) {
          res <- append(res, ind[i])
          start <-ind[i]
        }     
      }
      ind <- res
    }
    SNPnames <- textGrob(name, just="left", rot=90,
                         regionx[ind]-sqrt(2+0.5)*name_gap,
                         regiony[ind]+sqrt(2+0.5)*name_gap,
                         gp=gpar(cex=-0.006*nsnps+1.1, col="blue", fontfamily=family),
                         name="SNPnames", vp=vp)
  } else {
    symbols_pch <- NULL
    SNPnames <- NULL
  }
  label <- gTree(children=gList(symbols_pch, segments, diagonal, SNPnames), name ="label")
  label
}


