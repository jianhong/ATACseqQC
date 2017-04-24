#' @title shift 5' ends
#' @description shift the GAlignmentsLists by 5' ends. 
#' All reads aligning to the positive strand will be offset by +4bp, 
#' and all reads aligning to the negative strand will be offset -5bp by default.
#' @param gal An object of \link[GenomicAlignments]{GAlignmentsList}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @return An object of \link[GenomicAlignments]{GAlignments} with 5' end 
#' shifted reads.
#' @author Jianhong Ou
#' @export
#' @import S4Vectors
#' @import GenomicRanges
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' which <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
#' gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)
#' objs <- shiftGAlignmentsList(gal)
#' export(objs, "shift.bam")
shiftGAlignmentsList <- function(gal, positive=4L, negative=5L){
    stopifnot(is.integer(positive))
    stopifnot(is.integer(negative))
    stopifnot(is(gal, "GAlignmentsList"))
    stopifnot(length(gal)>0)
    stopifnot(all(elementNROWS(gal)<3))
    ## move the 5'end
    ## first reads is 5'end
    gal1 <- unlist(gal)
    gp <- rep(1, length(gal1))
    gp1 <- rep(seq_along(gal), elementNROWS(gal))
    gp[duplicated(gp1)] <- 2
    mcols(gal1)$MD <- NULL
    gal1[gp==1] <- shiftReads(gal1[gp==1], 
                              positive=positive,
                              negative=negative)
    names(gal1) <- mcols(gal1)$qname
    mcols(gal1)$isize[gp==2] <-
        sign(mcols(gal1)$isize[gp==2]) *
        abs(mcols(gal1)$isize[which(gp==2)-1])
    mcols(gal1)$mpos[gp==2] <- start(gal1)[which(gp==2)-1]
    mcols(gal1)$mpos[gp==1] <- start(gal1)[which(gp==1)+1]
    ## till now, gal1 must have mrnm, mpos, names and flag
    stopifnot(length(mcols(gal1)$mrnm)>0)
    stopifnot(length(mcols(gal1)$mpos)>0)
    stopifnot(length(mcols(gal1)$flag)>0)
    stopifnot(length(names(gal1))>0)
    return(gal1)
}