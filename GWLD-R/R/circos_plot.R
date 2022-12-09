#' circos snp
#'
#' @param data data.frame
#' @param ylim Range of data on y-axis
#' @param cell.padding only affect current track
#' @param track.height Height of the track. It is the percentage to the radius of the unit circles
#' @param cex.lab a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the defaul
#' @param bg.border Color for the border of the plotting regions
#' @param ... other graphical parameters
#' 
#' @import circlize
#' @importFrom graphics par rect text
#' @rdname snpcircos
#' @examples NULL
#' @export

circos.main <- function(data, ylim=c(0,1), cell.padding=c(0.02, 0, 0.02, 0),
                        track.height=0.1, cex.lab=0.6, bg.border=NA, ...) {
  ###输入的数据应该有3列,CHROM,POS,ID, 重要的是前两列
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
  #清楚circos参数
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
                             cex=cex.lab, niceFacing = TRUE)
               })

}

#' link snp
#'
#' @param data  data.frame
#' @param ylim  Range of data on y-axis
#' @param link.col link lines' color
#' @param track.height Height of the track
#' @param legend logical whether to show the legend
#' @param legend.pos legend position
#' @param legend.border Color for the border of the legend
#' @param legend.font font for the legend labels
#' @param legend.labels.cex Color for the border of the plotting regions
#' @param bg.border Color for the border of the plotting regions
#' @param bg.col Background color for the plotting regions
#' @param family font of legend label
#' @param lwd The line width
#' @param lty The line type
#' @param ... other graphical parameters
#'
#' @importFrom grDevices colorRampPalette grey.colors rainbow
#' @export
#'
#' @examples NULL
###将会构造成7列数据，CHROM_1, POS_1, ID_1 , CHROM_2, POS_2, ID_2, Value,
##根据value值的大小设置连线的颜色
circos.snplink <- function(data, ylim=c(0,1), link.col=NULL,
                            track.height=0.05, legend = TRUE, legend.pos=c(1, 0.5),
                            legend.border=NA, legend.font=1, legend.labels.cex = 0.5,
                            bg.border = "black", bg.col=NA, family=NULL,
                            lwd = par("lwd"), lty = par("lty"),...) {

  #第一列和第三列的染色体位位置必需是数值
  if (!all(grepl("chr",data[,1])) & !all(grepl("chr", data[,4]))) {
    data[,1] <- paste0("chr", data[,1])
    data[,4] <- paste0("chr", data[,4])
  }
  ##图例的显示的数据大小
  legend.min <- floor(round(min(data[,7]),digits = 2)*10)/10
  legend.max <- ceiling(round(max(data[,7]),digits = 2)*10)/10
  n_break <- seq(legend.min, legend.max, 0.05)
  n_color <- length(n_break)-1
  if(is.null(link.col)){
    link.col <- rainbow(4+n_color)[4+1:n_color]
  } else if(length(link.col)!=n_color) {
    stop(paste0("link.col parameter must input ", n_color, " colors"))
  } else {
    link.col <- link.col
  }
  #连接SNP，并根据值的大小设置颜色
  for (i in seq_len(nrow(data))) {
    colorkey <- cut(data[i, 7], n_break, labels=link.col)
    circos.link(data[i, 1], data[i, 2], data[i, 4], data[i, 5],
                col = as.character(colorkey),
                lwd = lwd,
                lty = lty)
  }
  ##显示SNP块
  circos.track(ylim=ylim, track.height=track.height,
               bg.col=bg.col, bg.border = bg.border,
               bg.lty = par("lty"), bg.lwd = par("lwd"))
  

  ##是否显示图例
  if (legend) {
    circos.legendlink(legend.pos[1], legend.pos[2], legend.limits=c(legend.min,legend.max),
                      legend.bin=1, legend.col=link.col, legend.border = legend.border,
                      legend.family=family, legend.cex=legend.labels.cex,
                      legend.lty = par("lty"), legend.lwd = par("lwd"),
                      legend.font=legend.font)
  }
  circos.clear() ##清除
}


#' circos legend
#'
#' @param x numeric
#' @param y numeric
#' @param legend.limits legend's limits;c(min, max)
#' @param legend.bin int; legend's bin
#' @param legend.col A specification for the default plotting color
#' @param legend.cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param legend.border Color for the border of the legend
#' @param legend.family The name of a font family for drawing text
#' @param legend.lty line width for borders and shading
#' @param legend.lwd line width for borders and shading
#' @param legend.font An integer which specifies which font to use for text
#' @param ... 	other graphical parameters
#'
#' @return NULL
#' @export
#'
circos.legendlink <- function(x=0.8, y=0.5, legend.limits=c(min,max),legend.bin=1, legend.col=legend.col,
                               legend.cex=0.5, legend.border=NA, legend.family=NULL,
                               legend.lty = par("lty"), legend.lwd = par("lwd"),
                               legend.font=1,...) {

  legend.min <- min(legend.limits)
  legend.max <- max(legend.limits)
  legend.rect <- (legend.max - legend.min)*10*2
  ##图例框
  rect(xleft = rep(x, legend.rect),
       ybottom = seq(y, y+0.025*(legend.rect-1), 0.025),
       xright = rep(x+0.06, legend.rect),
       ytop = seq(y+0.025,y+0.025*legend.rect, 0.025),
       col=legend.col, border = legend.border, lty = legend.lwd, lwd = legend.lwd)
  ##图例标签
  text(x=rep(x+0.085, legend.rect/2+1), y=seq(y, y+0.05*(legend.rect/2), 0.05*legend.bin),
       adj=c(0.5,0.5), labels = c(seq(legend.min, legend.max, 0.1*legend.bin)),
       cex=legend.cex, font = legend.font, family=legend.family)
}

