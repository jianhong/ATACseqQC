#' Transcription Start Site (TSS) Enrichment Score
#' @description TSS score is a raio between aggregate distribution of reads centered on TSSs and that flanking 
#' the corresponding TSSs. TSS score = the depth of TSS (each step within 1000 bp each side) / the depth of end flanks (100bp each end).
#' TSSE score = max(mean(TSS score in each window)).
#' @param obj an object of \link[GenomicAlignments:GAlignments-class]{GAlignments}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param upstream,downstream numeric(1) or integer(1). upstream and downstream of TSS. Default is 1000
#' @param endSize numeric(1) or integer(1). the size of the end flanks. Default is 100
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @export
#' @return A object of \link[GenomicRanges:GRanges-class]{GRanges} with TSS scores
#' @author Jianhong Ou
#' @references https://www.encodeproject.org/data-standards/terms/#enrichment
#' @examples  
#' library(GenomicRanges)
#' bamfile <- system.file("extdata", "GL1.bam", 
#'                        package="ATACseqQC", mustWork=TRUE)
#' gal1 <- readBamFile(bamFile=bamfile, tag=character(0), 
#'                     which=GRanges("chr1", IRanges(1, 1e6)), 
#'                     asMates=FALSE)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' tsse <- TSSEscore(gal1, txs)
TSSEscore <- function(obj, txs,
                      seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                      upstream=1000, downstream=1000, endSize=100){
 width <- 100
 step <- width
 pseudocount <- 0
  stopifnot(is(obj, "GAlignments"))
  if(length(obj)==0){
    obj <- loadBamFile(obj, minimal=TRUE)
  }
  stopifnot(is(txs, "GRanges"))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  cvg <- cvg + pseudocount
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  sel.center <- promoters(txs, upstream = upstream, downstream = downstream)
  sel.center$id <- seq_along(sel.center)
  sel.center.sw <- slidingWindows(sel.center, width = width, step = step)
  names(sel.center.sw) <- sel.center$id
  sel.center.s <- unlist(sel.center.sw, use.names = TRUE)
  sel.center.s$idx <- unlist(lapply(lengths(sel.center.sw), seq.int))
  
  sel.center.s <- split(sel.center.s, seqnames(sel.center.s))
  sel.center.s <- sel.center.s[names(cvg)]
  
  vws.center <- Views(cvg, sel.center.s)
  vms.center <- viewMeans(vws.center)
  
  ## do norm  
  #sel.center <- promoters(txs, upstream = upstream - endSize,
  #                        downstream = downstream - endSize
  #)
  sel.left.flank <- flank(sel.center, width=endSize, both=FALSE)
  sel.right.flank <- flank(sel.center, width=endSize, start=FALSE, both = FALSE)
  names(sel.left.flank) <- sel.left.flank$id
  names(sel.right.flank) <- sel.right.flank$id
  
  sel.left.flank <- split(sel.left.flank, seqnames(sel.left.flank))
  sel.left.flank <- sel.left.flank[names(cvg)]
  sel.right.flank <- split(sel.right.flank, seqnames(sel.right.flank))
  sel.right.flank <- sel.right.flank[names(cvg)]
  
  vws.left <- Views(cvg, sel.left.flank)
  vms.left <- viewMeans(vws.left)
  vws.right <- Views(cvg, sel.right.flank)
  vms.right <- viewMeans(vws.right)
  
  vms.m <- mapply(vms.center, sel.center.s, vms.left, vms.right, 
                  FUN = function(v, i, vl, vr){
    i <- i$idx
    id <- sort(union(names(vl), names(vr)))##make sure the left and right paired
    vl <- vl[id]
    vr <- vr[id]
    vl[is.na(vl)] <- vr[is.na(vl)]
    vr[is.na(vr)] <- vl[is.na(vr)]
    vl[is.na(vl)] <- pseudocount
    vr[is.na(vr)] <- pseudocount
    blk <- vl+vr
    names(blk) <- id
    keep <- blk>0 ## in case pseudocount is less than 1
    blk <- blk[keep]
    blk <- blk/2
    keep <- names(v) %in% names(blk)
    v <- v[keep]
    i <- i[keep]
    v <- v*endSize/blk[names(v)]/width
    tt <- table(i)
    rs <- rowsum(v, i) ## possible error for MAX_FLOAT
    if(length(rs)!=length(tt)) return(NULL)
    tt <- tt[rownames(rs)]
    tt[is.na(tt)] <- max(tt, na.rm = TRUE)
    names(tt) <- rownames(rs)
    rs <- cbind(rs, tt)
  }, SIMPLIFY = FALSE)
  
  tt <- do.call(rbind, lapply(vms.m, function(.ele) .ele[, 2]))
  vms.m <- do.call(rbind, lapply(vms.m, function(.ele) .ele[, 1]))
  vms.m <- colSums(vms.m)/colSums(tt)
  
  TSSE <- max(vms.m, na.rm = TRUE)
  return(list(values=vms.m, TSSEscore=TSSE))
}
