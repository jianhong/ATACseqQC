#' Transcription Start Site (TSS) Enrichment Score
#' @description TSS score is a raio between aggregate distribution of reads centered on TSSs and that flanking 
#' the corresponding TSSs. TSS score = the depth of TSS (1000 bp each side) / the depth of end flanks (100bp each end).
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
  stopifnot(is(obj, "GAlignments"))
  stopifnot(is(txs, "GRanges"))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  sel.center <- promoters(txs, upstream = upstream, downstream = downstream)
  sel.center$id <- seq_along(sel.center)
  
  sel.left.flank <- flank(sel.center, width=endSize, both=FALSE)
  sel.right.flank <- flank(sel.center, width=endSize, start=FALSE, both = FALSE)
  
  sel.center <- split(sel.center, seqnames(sel.center))
  sel.center <- sel.center[names(cvg)]
  sel.left.flank <- split(sel.left.flank, seqnames(sel.left.flank))
  sel.left.flank <- sel.left.flank[names(cvg)]
  sel.right.flank <- split(sel.right.flank, seqnames(sel.right.flank))
  sel.right.flank <- sel.right.flank[names(cvg)]
  
  sel.center.id <- unlist(sel.center)$id
  sel.left.id <- unlist(sel.left.flank)$id
  sel.right.id <- unlist(sel.right.flank)$id
  stopifnot(identical(sel.left.id, sel.center.id))
  stopifnot(identical(sel.right.id, sel.center.id))
  
  vws.center <- Views(cvg, sel.center)
  vms.center <- viewMeans(vws.center)
  vms.center <- unlist(vms.center)
  vms.center <- vms.center[order(sel.center.id)]
  
  vws.left.flank <- Views(cvg, sel.left.flank)
  vws.left.flank <- viewSums(vws.left.flank)
  vws.left.flank <- unlist(vws.left.flank)
  vws.left.flank <- vws.left.flank[order(sel.left.id)]
  
  vws.right.flank <- Views(cvg, sel.right.flank)
  vws.right.flank <- viewSums(vws.right.flank)
  vws.right.flank <- unlist(vws.right.flank)
  vws.right.flank <- vws.right.flank[order(sel.right.id)]
  
  vws.flank <- (vws.left.flank + vws.right.flank)/2/endSize
  
  TSS.score <- vms.center/vws.flank
  
  sel.center <- unlist(sel.center)
  sel.center <- sel.center[order(sel.center.id)]
  sel.center$TSS.mean <- vms.center
  sel.center$flank.mean <- vws.flank
  sel.center$TSS.enrichment.score <- TSS.score
  return(sel.center)
}