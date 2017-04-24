#' @title fragment size distribution
#' @description estimate the fragment size of bams
#' @param bamFiles A vector of characters indicates the file names of bams.
#' @param bamFiles.labels labels of the bam files, used for pdf file naming.
#' @param ylim numeric(2). ylim of the histogram.
#' @param logYlim numeric(2). ylim of log-transformed histogram for the insert.
#' @return Invisible fragment length distribution list.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam idxstatsBam
#' @importFrom graphics axis par plot
#' @import GenomicRanges
#' @export
#' @author Jianhong Ou
#' @examples
#' bamFiles <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' bamFiles.labels <- "GL1"
#' fragSizeDist(bamFiles, bamFiles.labels,
#'              ylim=c(0, 1e4), logYlim=log10(c(5e-3, 2)))

fragSizeDist <- function(bamFiles, bamFiles.labels, ylim=NULL,
                         logYlim=NULL){
  opar <- par(c("fig", "mar"))
  on.exit(par(opar))
  summaryFunction <- function(seqname, seqlength, bamFile, ...) {
    param <-
      ScanBamParam(what=c('isize'),
                   which=GRanges(seqname, IRanges(1, seqlength)),
                   flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                    isUnmappedQuery=FALSE,
                                    isNotPassingQualityControls = FALSE))
    table(abs(scanBam(bamFile, ..., param=param)[[1]][["isize"]]))
  }

  idxstats <- unique(do.call(rbind, lapply(bamFiles, function(.ele)
    idxstatsBam(.ele)[, c("seqnames", "seqlength")])))
  seqnames <- as.character(idxstats[, "seqnames"])
  seqlen <- as.numeric(idxstats[, "seqlength"])
  fragment.len <- lapply(bamFiles, summaryFunction,
                         seqname=seqnames, seqlength=seqlen,
                         SIMPLIFY=FALSE)

  names(fragment.len) <- bamFiles.labels

  minor.ticks.axis <- function(ax,n=9,t.ratio=0.5,mn,mx,...){

    lims <- par("usr")
    lims <- if(ax %in% c(1,3)) lims[1:2] else lims[3:4]

    major.ticks <- pretty(lims,n=5)
    if(missing(mn)) mn <- min(major.ticks)
    if(missing(mx)) mx <- max(major.ticks)

    major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

    labels <- sapply(major.ticks,function(i)
      as.expression(bquote(10^ .(i)))
    )
    axis(ax,at=major.ticks,labels=labels,
         las=ifelse(ax %in% c(2, 4), 2, 1), ...)

    n <- n+2
    minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
    minors <- minors[-c(1,n)]

    minor.ticks = c(outer(minors,major.ticks,`+`))
    minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


    axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
  }

  null <- mapply(function(frag.len, frag.name){
    x <- 1:1010
    frag.len <- frag.len[match(x, names(frag.len))]
    frag.len[is.na(frag.len)] <- 0
    y <- 100 * frag.len / sum(frag.len)
    y <- as.numeric(y)
    names(y) <- x
    par(mar=c(5, 5, 4, 2) +.1)
    plot(x, y*10^3, main=paste(frag.name, "fragment sizes"),
         xlim=c(0, 1010), ylim=ylim,
         xlab="Fragment length (bp)",
         ylab=expression(Normalized ~ read ~ density ~ x ~ 10^-3),
         type="l")
    par(fig=c(.4, .95, .4, .95), new=TRUE)
    plot(x, log10(y), xlim=c(0, 1010), ylim=logYlim,
         xlab="Fragment length (bp)", ylab="Norm. read density",
         type="l", yaxt="n")
    minor.ticks.axis(2)
  }, fragment.len, names(fragment.len))

  return(invisible(fragment.len))
}
