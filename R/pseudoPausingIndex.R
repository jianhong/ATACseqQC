#' Simulation pausing index
#' @description The polII pausing index is the ratio of the reads density at 
#' the 5' end of the gene to that in the gene body. This function will simulate
#' the pausing index by open chromatin coverage instead of PolII signaling.
#' The pausing index is a raio between aggregate distribution of reads in TSS
#' and that elongating gene bodys. The default PI = average coverage of TSS
#' (+20 to +60bp) / the average coverage of in transcripts (+60bp to 100bp).
#' @param obj an object of \link[GenomicAlignments:GAlignments-class]{GAlignments}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param nascentRegion,pausedRegion,elongationRegion numeric(2) or integer(2).
#' The start and end position of the pre-initiation complex, paused complex and
#' productive elongation.
#' @param pseudocount numeric(1) or integer(1). Pseudocount. Default is 1.
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @export
#' @return A object of \link[GenomicRanges:GRanges-class]{GRanges} with pseudo
#' pausing index.
#' @author Jianhong Ou
#' @references https://doi.org/10.1098%2Frsob.210220
#' @examples  
#' library(GenomicRanges)
#' bamfile <- system.file("extdata", "GL1.bam", 
#'                        package="ATACseqQC", mustWork=TRUE)
#' gal1 <- readBamFile(bamFile=bamfile, tag=character(0), 
#'                     which=GRanges("chr1", IRanges(1, 1e6)), 
#'                     asMates=FALSE)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' ppi <- pseudoPausingIndex(gal1, txs)
pseudoPausingIndex <- function(obj, txs,
                      seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                      nascentRegion=c(-20, 20),
                      pausedRegion=c(21, 60),
                      elongationRegion=c(61, 100),
                      pseudocount=1L){
  stopifnot(is(obj, "GAlignments"))
  if(length(obj)==0){
    obj <- loadBamFile(obj, minimal=TRUE)
  }
  stopifnot(is(txs, "GRanges"))
  stopifnot(length(nascentRegion)==2)
  stopifnot(length(pausedRegion)==2)
  stopifnot(length(elongationRegion)==2)
  stopifnot(is.numeric(nascentRegion))
  stopifnot(is.numeric(pausedRegion))
  stopifnot(is.numeric(elongationRegion))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  if(pseudocount!=0) cvg <- cvg + pseudocount
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  TSS <- promoters(txs, upstream = 0L, downstream = 1L)
  str_neg <- as.character(strand(TSS))=="-"
  nascentRegionGR <- pausedRegionGR <- elongationRegionGR <- TSS
  getRegion <- function(x, region){
    a <- min(region)[1]
    b <- max(region)[1]
    if(a>0) {
      w <- diff(region)+1
    }else{
      a <- a+1
      w <- diff(region)
    }
    x[!str_neg] <-
      shift(x[!str_neg], shift = a)
    width(x[!str_neg]) <- w
    start(x[str_neg]) <- start(x[str_neg]) - 
      max(region)[1]
    end(x[str_neg]) <- end(x[str_neg]) - 
      a
    x
  }
  nascentRegionGR <- getRegion(TSS, nascentRegion)
  pausedRegionGR <- getRegion(TSS, pausedRegion)
  elongationRegionGR <- getRegion(TSS, elongationRegion)
  getCvg <- function(x){
    x <- split(x, seqnames(x))
    x <- x[names(cvg)]
    vws <- Views(cvg, x)
    vms <- viewMeans(vws)
    vws <- mapply(vws, vms, FUN=function(vw, vm){
      ir <- ranges(vw)
      mcols(ir)$score <- vm
      ir
    }, SIMPLIFY = FALSE)
    gr <- unlist(x)
    gr$score <- mcols(unlist(IRangesList(vws)))$score
    gr
  }
  nascentRegionVWS <- getCvg(nascentRegionGR)
  pausedRegionVWS <- getCvg(pausedRegionGR)
  elongationRegionVWS <- getCvg(elongationRegionGR)
  ppi <- getCvg(txs)
  ppi$nascentRegion <- nascentRegionVWS$score
  ppi$pausedRegion <- pausedRegionVWS$score
  ppi$elongationRegion <- elongationRegionVWS$score
  ppi$pseudoPI <- (pausedRegionVWS$score+pseudocount)/
    (elongationRegionVWS$score+pseudocount)
  
  return(ppi)
}

