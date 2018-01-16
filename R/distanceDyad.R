#' Distance of potential nucleosome dyad
#' @description Calculate the distance of potential nucleosome dyad and the linear model for V.
#' @param vPlotOut The output of \link{vPlot}.
#' @param fragLenRanges  A numeric vector (length=3) for fragment size of nucleosome free and mono-nucleosome.
#' Default c(60, 180, 250).
#' @param draw Plot the results or not. Default TRUE.
#' @param ... Prameters could be passed to plot.
#' @return an invisible list with distance of nucleosome and the linear model.
#' @importFrom stats loess.smooth lm
#' @importFrom graphics abline
#' @export
#' @author Jianhong Ou
#' @seealso \link{vPlot}
#' @examples 
#'bamfile <- system.file("extdata", "GL1.bam",
#'                        package="ATACseqQC")
#'library(MotifDb)
#'CTCF <- query(MotifDb, c("CTCF"))
#'CTCF <- as.list(CTCF)
#'library(BSgenome.Hsapiens.UCSC.hg19)
#'vp <- vPlot(bamfile, pfm=CTCF[[1]],
#'      genome=Hsapiens,
#'      min.score="95%", seqlev="chr1",
#'      draw=FALSE)
#'distanceDyad(vp)
#'      
distanceDyad <- function(vPlotOut, fragLenRanges=c(60, 180, 250), draw=TRUE, ...){
  stopifnot(length(fragLenRanges)==3)
  stopifnot(is(fragLenRanges, "numeric"))
  if(!identical(colnames(vPlotOut), c("distanceToBindingSite", "FragmentLength"))){
    stop("vPlotOut must be the output of vPlot.")
  }
  vPlotOut <- vPlotOut[vPlotOut$FragmentLength>=fragLenRanges[1] & 
                         vPlotOut$FragmentLength<=fragLenRanges[3], , drop=FALSE]
  if(length(unique(vPlotOut$distanceToBindingSite))<300){
    return(c(distanceNucl=NA, bindingWidth=NA))
  }
  data4nucl <- vPlotOut[vPlotOut$FragmentLength>=fragLenRanges[2] & 
                          vPlotOut$FragmentLength<=fragLenRanges[3], , drop=FALSE]
  data4nucl$FragmentLength <- 1
  data4nucl$distanceToBindingSite <- abs(data4nucl$distanceToBindingSite)
  data4nucl <- rowsum(data4nucl$FragmentLength, data4nucl$distanceToBindingSite)
  delta <- diff(range(data4nucl[, 1]))/4
  data4nucl <- tryCatch(loess.smooth(as.numeric(rownames(data4nucl)), data4nucl[, 1], 
                                     span=1/10, evaluation = length(data4nucl)),
                        error=function(e) NULL)
  if(length(data4nucl)==2){
    data4nucl.peak <- peakdet(data4nucl$y, delta = delta)
    if(length(data4nucl.peak$peakpos)>0){
      data4nucl.peak$x <- data4nucl$x[data4nucl.peak$peakpos]
      data4nucl.peak$y <- data4nucl$x[data4nucl.peak$peakpos]
      data4ncul.peak <- data4nucl.peak$x[which.max(data4nucl.peak$y)]
    }else{
      data4ncul.peak <- max(data4nucl$x)
    }
  }else{
    data4ncul.peak <- max(data4nucl$x)
  }
  
  distanceNucl <- c(-1 * data4ncul.peak, data4ncul.peak)
  
  data4bw.bk.1 <- data4bw <- vPlotOut[vPlotOut$distanceToBindingSite > -250 & 
                                        vPlotOut$distanceToBindingSite < 250,, drop=FALSE]
  data4bw$FragmentLength <- floor(data4bw$FragmentLength/10)*10 + 5
  data4bw <- split(data4bw$distanceToBindingSite, data4bw$FragmentLength)
  data4bw <- lapply(data4bw, table)
  data4bw.names <- unique(unlist(sapply(data4bw, names), use.names = FALSE))
  data4bw <- do.call(cbind, lapply(data4bw, function(.ele) .ele[match(data4bw.names, names(.ele))]))
  data4bw[is.na(data4bw)] <- 0
  data4bw <- as.data.frame(data4bw)
  data4bw.names <- as.numeric(data4bw.names)
  suppressWarnings({
    data4bw <- lapply(data4bw, function(.ele){
      .ele <- tryCatch(loess.smooth(data4bw.names, .ele, span=1/10, evaluation = length(data4bw.names)),
                       error=function(e) NULL)
      if(length(.ele)==0){
        return(NULL)
      }
      .peak <- peakdet(.ele$y)
      .valley <- .ele$x[.peak$valleypos]
      vl <- .valley[.valley < 0]
      vr <- .valley[.valley > 0]
      if(length(vl)) {
        vl <- vl[which.max(vl)]
        #if(vl < distanceNucl[1]) return(NULL)
      }else{
        return(NULL)
      }
      if(length(vr)){
        vr <- vr[which.min(vr)]
        #if(vr > distanceNucl[2]) return(NULL)
      }else{
        return(NULL)
      }
      c(vl=vl, vr=vr)
    })
  })
  data4bw <- do.call(cbind, data4bw)
  data4bw <- data4bw[, apply(data4bw, 2, function(.ele) !any(is.na(.ele)))]
  data4bw.colsum <- colSums(data4bw)
  data4bw.mean <- mean(abs(data4bw))
  data4bw <- data4bw[, abs(data4bw.colsum) < data4bw.mean/5, drop=FALSE]
  data4bw <- colMeans(abs(data4bw))
  data4bw <- data4bw[seq.int(which.max(data4bw))]
  
  if(length(data4bw)>=3){
    data4bw <- data.frame(x=data4bw, y=as.numeric(names(data4bw)))
    lm <- lm(y ~ 0+x, data4bw) ## cross zero
    if(summary(lm)$r.squared>.9){
      coefficient <- summary(lm)$coefficients[1]
      ## rotate matrix
      coefficients <- coefficient + seq(-0.5, .5, by = 0.02)
      angels <- -atan(coefficients)
      data4bw.bk.rot <- lapply(angels, function(angel){
        R <- matrix(c(cos(angel), sin(angel), -sin(angel), cos(angel)), nrow = 2)
        data4bw.bk <- data4bw.bk.1
        data4bw.bk$distanceToBindingSite <- abs(data4bw.bk$distanceToBindingSite)
        data4bw.bk <- t(R %*% t(as.matrix(data4bw.bk)))
        data4bw.bk <- round(data4bw.bk, digits = 0)
        list(a=data4bw.bk, b=sum(data4bw.bk[, 2]==0))
      })
      data4bw.bk.rot.min <- which.min(sapply(data4bw.bk.rot, function(.ele) .ele$b))
      data4bw.bk.rot.min <- data4bw.bk.rot.min[which.min(abs(data4bw.bk.rot.min-25))]
      coefficient <- coefficients[data4bw.bk.rot.min]
      data4bw.bk <- data4bw.bk.rot[[data4bw.bk.rot.min]]$a
      data4bw.bk <- sapply(split(data4bw.bk[, 1], data4bw.bk[, 2]), mean, na.rm=FALSE)
      data4bw.bk <- tryCatch(loess.smooth(x=as.numeric(names(data4bw.bk)), y=data4bw.bk, 
                                 span=1/10, evaluation = length(data4bw.bk)),
                             error=function(e) NULL)
      if(length(data4bw.bk)==2){
        data4bw.bk <- as.data.frame(data4bw.bk)
        data4bw.bk.peak <- peakdet(data4bw.bk$y)
        if(length(data4bw.bk.peak$valleypos)>0){
          data4bw.bk.bottom <- min(data4bw.bk$x[data4bw.bk.peak$valleypos])
        }else{
          data4bw.bk.bottom <- 0
        }
        if(length(data4bw.bk.peak$peakpos)>0){
          zero <- which(data4bw.bk$x==data4bw.bk.bottom)
          data4bw.bk.peak <- data4bw.bk.peak$peakpos[which.min(abs(data4bw.bk.peak$peakpos - zero))]
          if(data4bw.bk.peak<zero){
            Apoint <- data4bw.bk.peak
            Bpoint <- data4bw.bk.peak + zero
          }else{
            Apoint <- zero - (data4bw.bk.peak - zero) 
            Bpoint <- data4bw.bk.peak
          }
          Apoint <- max(c(1, Apoint))
          Bpoint <- min(c(Bpoint, nrow(data4bw.bk)))
          data4bw.bk <- data4bw.bk[Apoint:Bpoint, ]
        }
        q50 <- quantile(data4bw.bk$y, probs = .25)
        data4bw.bk.left <- data4bw.bk[data4bw.bk$x<data4bw.bk.bottom & data4bw.bk$y<q50, "x"]
        data4bw.bk.right <- data4bw.bk[data4bw.bk$x>data4bw.bk.bottom & data4bw.bk$y<q50, "x"]
        if(length(data4bw.bk.left)>0 && length(data4bw.bk.right)>0){
          bindingWidth <- max(data4bw.bk.right) - min(data4bw.bk.left)
        }else{
          bindingWidth <- NA
        }
      }else{
        bindingWidth <- NA
      }
    }else{
      coefficient <- NA
      bindingWidth <- NA
    }
  }else{
    lm <- NA
    coefficient <- NA
    bindingWidth <- NA
  }
  
  if(draw){
    plot(vPlotOut, ...)
    #lines(data4nucl, col="red")
    abline(v=distanceNucl, col="red", lty="dotdash")
    if(!is.na(coefficient)){
      abline(a=0, b=-1*coefficient, col="red", lty="dotdash")
      abline(a=0, b=coefficient, col="red", lty="dotdash")
      if(!is.na(bindingWidth)){
        abline(a=bindingWidth, b=coefficient, col="red", lty="dotdash")
        abline(a=-bindingWidth, b=coefficient, col="red", lty="dotdash")
        abline(a=bindingWidth, b=-coefficient, col="red", lty="dotdash")
        abline(a=-bindingWidth, b=-coefficient, col="red", lty="dotdash")
      }
    }
  }
  return(invisible(list(distanceNucl=diff(distanceNucl), lm=lm, bindingWidth=bindingWidth)))
}
