#' circos snp
#'
#' @param data data.frame
#' @param ylim Range of data on y-axis
#' @param cell.padding only affect current track
#' @param track.height Height of the track. It is the percentage to the radius of the unit circles
#' @param cex.lab a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the defaul
#' @param bg.border Color for the border of the plotting regions
#' @param family the name of a font family for drawing text
#' @param ... other graphical parameters
#' 
#' @import circlize
#' @importFrom graphics par rect text
#' @rdname snpcircos
#' @examples NULL
#' @export

circos.ideogram <- function(data, ylim=c(0,1), cell.padding=c(0.02, 0, 0.02, 0),
                            track.height=0.1, cex.lab=0.6, bg.border=NA, family=NULL, ...) {
  ###输入的数据应该有3列,CHROM,POS,ID, 重要的是前两列可以直接拿vcf文件中前5列
  if(!all(grepl("chr", data[,1]))) {
    df <- data.frame()
    i = 1
    for (i in unique(data[,1])) {
      tmp_dat <- data[data[,1]==i,]
      chr <- paste0("chr", i)
      end <- tmp_dat[nrow(tmp_dat), 2]
      df <- rbind(df, c(chr,0,end))
    }
    colnames(df) <- c("chr", "start", "end")
    rownames(df) <- df[,1]
  } else {
    data <- data
  }
  #设置字体
  family <- ifelse(is.null(family), "serif", family)
  #清除circos参数
  circos.clear()
  circos.par(cell.padding = cell.padding)
  ##必须使用因子，不然顺序会改变
  circos.initialize(factors = factor(df[,1], levels=df[,1]), xlim =df[,2:3])
  ##坐标轴和标签
  circos.track(ylim = ylim, track.height=track.height, bg.border = bg.border,
               panel.fun = function(x, y, ...) {
                 sector.index <- get.current.sector.index()
                 #显示坐标轴
                 circos.axis(h="bottom", labels = F)
                 chr.index <- which(get.all.sector.index()==sector.index)
                 #显示染色体号
                 xlim <- get.cell.meta.data("xlim")
                 circos.text(mean(xlim), 0.6, labels = chr.index, facing= "downward",
                             cex=cex.lab, niceFacing = TRUE, family=family)
               })
}

#' link snp
#'
#' @param data  A data frame with 7 columns;Chr_1, Pos_1, ID_1, Chr_2, Pos_2, ID_2 value
#' @param ylim  Range of data on y-axis
#' @param color.theme The color theme in the legend
#' @param track.height Height of the track
#' @param legend.show logical whether to show the legend
#' @param legend.position legend position
#' @param legend.fontsize the font size of legend
#' @param legend.limits A numeric vector of length two providing limits of the scale
#' @param legend.breaks A numeric vector of positions
#' @param legend.labels A vector giving labels (must be same length as legend.breaks)
#' @param legend.n.breaks An integer guiding the number of major breaks
#' @param bg.border Color for the border of the plotting regions
#' @param bg.col Background color for the plotting regions
#' @param family The font family of legend
#' @param lwd The line width
#' @param lty The line type
#' @param ... other graphical parameters
#' 
#' @import viridisLite
#' @importFrom grDevices colorRampPalette grey.colors rainbow col2rgb rgb 
#' @export
#'
#' @examples NULL

#数据有7列，CHROM_1, POS_1, ID_1 , CHROM_2, POS_2, ID_2, Value,
#根据value值的大小设置连线的颜色
circos.linksnp <- function(data, ylim=c(0,1), color.theme="viridis", track.height=0.05, legend.show=TRUE,
                           legend.position=c(0.88, 0.75), legend.fontsize=12, legend.limits=NULL, 
                           legend.breaks=NULL, legend.labels=NULL, legend.n.breaks=NULL, family=NULL, 
                           bg.border="black", bg.col=NULL, lwd = par("lwd"), lty = par("lty"), ... ) {
  #第一列和第三列的染色体位置必需是数值
  if(!all(grepl("chr", data[,1])) && !all(grepl("chr", data[,4]))) {
    data[,1] <- paste0("chr", data[,1])
    data[,4] <- paste0("chr", data[,4])
  }
  #图例使用多少种颜色
  n.breaks <- ifelse(is.null(legend.n.breaks), 100,  legend.n.breaks)
  #使用viridisLite包中的颜色
  if(color.theme %in% c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako")) {
    link.colors <- do.call(color.theme, list(n.breaks, alpha = 1, begin = 0, end = 1, direction = -1))
  } 
  if (is.null(legend.limits)) {
    #获取数据中的最大值和最小值
    limits <- c(floor(min(data[, 7])*10)/10, ceiling(max(data[, 7])*10)/10)
  } else {
    limits <- legend.limits
  }
  color_breaks <- seq(limits[1], limits[2], length.out=n.breaks+1)
  colorkey <- cut(data[, 7], color_breaks, labels=link.colors, include.lowest=T)
  # #转换成RGB,归一化至[0,1]
  # rgb_col <- data.frame(col2rgb(colorkey)/255)
  # #设置透明度,归一化至[0,1]
  # alpha <- as.numeric(colorkey)/n.breaks
  # #连接SNP, 并根据值的大小设置颜色
  # for (i in 1:nrow(data)) {
  #   circos.link(data[i, 1], data[i, 2], data[i, 4], data[i, 5],
  #               col = rgb(rgb_col[1, i], rgb_col[2, i], rgb_col[3, i], alpha[i]),
  #               lwd = lwd,
  #               lty = lty)
  # }
  for (i in 1:nrow(data)) {
    circos.link(data[i, 1], data[i, 2], data[i, 4], data[i, 5],
                col = as.character(colorkey[i]),
                lwd = lwd,
                lty = lty)
  }
  ##显示SNP块颜色
  if(is.null(bg.col)) {
    bg.col <- rainbow(length(get.all.sector.index()))
  } else {
    bg.col <- bg.col
  }
  #显示SNP块
  circos.track(ylim=ylim, track.height=track.height,
               bg.col=bg.col, bg.border=bg.border,
               bg.lty=par("lty"), bg.lwd=par("lwd"))
  
  #设置字体
  family <- ifelse(is.null(family), "serif", family)
  #是否显示图例
  if(legend.show) {
    if(is.null(legend.breaks)) {
      legend.breaks <- seq(limits[1], limits[2], length.out=5)
    } else {
      legend.breaks <- legend.breaks
    }
    if(is.null(legend.labels)) {
      labels <- legend.breaks
    } else {
      labels <- legend.labels
    }
    legend <- circos.legendlink(x=legend.position[1], y=legend.position[2], legend.limits=limits,
                                legend.breaks=legend.breaks,legend.n.breaks=n.breaks, legend.col=link.colors, 
                                legend.labels=labels, legend.family=family, legend.fontsize=legend.fontsize, ...)
    grid.draw(legend)
  }
  circos.clear()##清除参数
}

#' circos legend
#'
#' @param x numeric
#' @param y numeric
#' @param legend.limits A numeric vector of length two providing limits of the scale
#' @param legend.n.breaks An integer guiding the number of major breaks.
#' @param legend.breaks A numeric vector of positions
#' @param legend.labels A vector giving labels (must be same length as legend.breaks)
#' @param legend.col A specification for the default plotting color
#' @param legend.fontsize numeric; the font size of legend
#' @param legend.family The font family of legend
#' @param ... 	other graphical parameters
#'
#' @return NULL
#' @export
#'

circos.legendlink <- function(x=x, y=y, legend.limits=c(min, max), legend.n.breaks=legend.n.breaks,
                              legend.breaks=legend.breaks, legend.labels=legend.labels, legend.col=legend.col, 
                              legend.fontsize=legend.fontsize, legend.family=legend.family, ...) {
  #图例标签位置
  legend_label_pos <- seq(0, 1, length.out=legend.n.breaks)
  breaks <- seq(legend.limits[1], legend.limits[2], length.out=legend.n.breaks+1)
  pos <- cut(legend.breaks, breaks, labels=FALSE, right=TRUE, include.lowest=T)
  ##移除可能没有匹配到的标签
  ind <- pos[!is.na(pos)]
  legend.labels <- legend.labels[!is.na(pos)]
  
  right <- rep(1, legend.n.breaks)
  top <- rep((1:legend.n.breaks)/legend.n.breaks, each=1)
  ##图例视图
  legend_vp <- viewport(x=x, y=y, width=unit(0.035, "snpc"), height=unit(0.25, "snpc"), name = "legend_vp")
  legend_rect <- rectGrob(x=right, y=top, width=1, height=1/legend.n.breaks,
                          just=c("right", "top"), name="legend_rect", gp = gpar(fill=legend.col, col=NA))
  legend_ticks <- segmentsGrob(x0=rep(c(0,1), each=length(legend.labels)), 
                               y0=rep(legend_label_pos[ind], length(legend.labels)), 
                               x1=rep(c(0.15, 0.85), each=length(legend.labels)), 
                               y1=rep(legend_label_pos[ind], length(legend.labels)),
                               gp=gpar(col="white", lwd=0.7),
                               name = "legend_ticks")
  legend_box <- linesGrob(x=c(0,1,1,0,0), y=c(0,0,1,1,0), name = "legend_box", gp=gpar(col="white"))
  ##调整宽度
  width <- convertUnit(unit(1.5, "lines"), "npc", valueOnly = T) 
  legend_text <- textGrob(legend.labels, x=rep(1, length(legend.labels))+width, y=legend_label_pos[ind], just = c("left"),
                          gp=gpar(fontsize=legend.fontsize, fontfamily=legend.family),
                          name="legend_text")
  
  legend <- gTree(children = gList(legend_rect, legend_ticks, legend_box, legend_text), vp=legend_vp, name="legend")
  legend
}