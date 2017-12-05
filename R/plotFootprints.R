#' Plots a footprint estimated by Centipede
#'
#' Visualizing the footprint profile
#'
#' @param Profile A vector with the profile estimated by CENTIPEDE
#' @param Mlen Length of the motif for drawing vertical lines delimiting it
#' @param xlab Label of the x axis
#' @param ylab Label for the y axis
#' @param legTitle Title for one of the plot corners
#' @param xlim xlim
#' @param ylim ylim
#' @param newpage Plot the figure in a new page?
#' @param motif a pfm object.
#' @param segmentation the segmentation position and abundance
#' @importFrom grid grid.newpage viewport plotViewport pushViewport upViewport
#' gpar grid.xaxis grid.yaxis convertX convertY 
#' unit grid.legend grid.text grid.lines grid.segments
#' @importFrom motifStack plotMotifLogoA
#' @importClassesFrom motifStack pfm
#' @return Null.
#' @author Jianhong Ou
#' @examples
#' library(MotifDb)
#' CTCF <- query(MotifDb, c("CTCF"))
#' CTCF <- as.list(CTCF)
#' motif <- new("pfm", mat=CTCF[[1]], name="CTCF")
#' ATACseqQC:::plotFootprints(Profile=sample.int(500), 
#'                         Mlen=ncol(CTCF[[1]]), motif=motif)
plotFootprints <- function (Profile, Mlen = 0,
                            xlab = "Dist. to motif (bp)",
                            ylab = "Cut-site probability",
                            legTitle, xlim, ylim,
                            newpage = TRUE, motif, segmentation)
{
  stopifnot(is(motif, "pfm"))
  if(newpage) grid.newpage()

  S <- length(Profile)
  W <- ((S/2) - Mlen)/2
  vp <- plotViewport(margins=c(5.1, 5.1, 4.1, 2.1), name="plotRegion")
  pushViewport(vp)
  if(missing(xlim)){
    xlim <- c(0, S/2+1)
  }
  if(missing(ylim)){
    ylim <- c(0, max(Profile) * 1.12)
  }
  vp1 <- viewport(y=.4, height=.8,
                  xscale=xlim,
                  yscale=ylim,
                  name="footprints")
  pushViewport(vp1)
  grid.lines(x=1:(S/2),
             y=Profile[1:(S/2)],
             default.units="native",
             gp=gpar(lwd = 2, col = "darkblue"))
  grid.lines(x=1:(S/2),
             y=Profile[(S/2 + 1):S],
             default.units="native",
             gp=gpar(lwd = 2, col = "darkred"))
  if(!missing(segmentation)){
    if(length(segmentation)==4){
      ## plot the guide line for the segmentation
      grid.segments(x0=c(0, segmentation[1], W, W+Mlen, S/2-segmentation[1]),
                    x1=c(segmentation[1], W, W+Mlen, S/2-segmentation[1], S/2),
                    y0=c(segmentation[2], segmentation[3], segmentation[4], segmentation[3], segmentation[2]),
                    y1=c(segmentation[2], segmentation[3], segmentation[4], segmentation[3], segmentation[2]),
                    default.units = "native",
                    gp=gpar(lwd =2, col = "red", lty = 2))
    }
  }
  grid.xaxis(at = c(seq(1, W, length.out = 3),
                    W + seq(1, Mlen),
                    W + Mlen + seq(1, W, length.out = 3)),
             label = c(-(W + 1 - seq(1, W + 1, length.out = 3)),
                       rep("", Mlen),
                       seq(0, W, len = 3)))
  grid.yaxis()
  grid.lines(x=c(W, W, 0), y=c(0, max(Profile), ylim[2]),
             default.units="native", gp=gpar(lty=2))
  grid.lines(x=c(W + Mlen + 1, W + Mlen + 1, S/2),
             y=c(0, max(Profile), ylim[2]),
             default.units="native", gp=gpar(lty=2))
  upViewport()
  vp2 <- viewport(y=.9, height=.2,
                  xscale=c(0, S/2+1),
                  name="motif")
  pushViewport(vp2)
  motifStack::plotMotifLogoA(motif)
  upViewport()
  upViewport()
  grid.text(xlab, y=unit(1, 'lines'))
  grid.text(ylab, x=unit(1, 'line'), rot = 90)
  if(missing(legTitle)){
    legvp <- viewport(x=unit(1, "npc")-convertX(unit(1, "lines"), unitTo="npc"),
                      y=unit(1, "npc")-convertY(unit(1, "lines"), unitTo="npc"),
                      width=convertX(unit(14, "lines"), unitTo="npc"),
                      height=convertY(unit(3, "lines"), unitTo="npc"),
                      just=c("right", "top"), name="legendWraper")
    pushViewport(legvp)
    grid.legend(labels=c("For. strand", "Rev. strand"),
                gp=gpar(lwd=2, lty=1, col=c("darkblue", "darkred")))
    upViewport()
  } else {
    grid.text(legTitle,
              y=unit(1, "npc")-convertY(unit(1, "lines"), unitTo="npc"),
              gp=gpar(cex=1.2, fontface="bold"))
  }
  return(invisible())
}
