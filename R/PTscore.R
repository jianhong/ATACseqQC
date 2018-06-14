#' Promoter/Transcript body (PT) score
#' @description PT score is calculated for coverage of promoter divided by the coverage of transcripts body.
#' PT score will show if the signal is enriched in promoters.
#' @param obj an object of \link[GenomicAlignments:GAlignments-class]{GAlignments}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param upstream numeric(1) or integer(1). Start position of promoter. Default is 2000
#' @param downstream numeric(1) or integer(1). End position of promoter. Default is 500
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @export
#' @return A object of \link[GenomicRanges:GRanges-class]{GRanges} with PT scores
#' @author Jianhong Ou
#' @examples  
#' library(GenomicRanges)
#' bamfile <- system.file("extdata", "GL1.bam", 
#'                        package="ATACseqQC", mustWork=TRUE)
#' gal1 <- readBamFile(bamFile=bamfile, tag=character(0), 
#'                     which=GRanges("chr1", IRanges(1, 1e6)), 
#'                     asMates=FALSE)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' pt <- PTscore(gal1, txs)
PTscore <- function(obj, txs,
                     seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                     upstream=2000, downstream=500){
  stopifnot(is(obj, "GAlignments"))
  stopifnot(is(txs, "GRanges"))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  obj <- promoters(obj, upstream = 0, downstream = 1)
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  pro <- promoters(txs, upstream = upstream, downstream = downstream)
  body <- shift(pro, shift = upstream + downstream)
  body[strand(pro)=="-"] <- shift(body[strand(pro)=="-"], 
                                  shift = -1 * (upstream + downstream))
  pro$source <- "promoter"
  body$source <- "transcript"
  pro$oid <- seq_along(pro)
  body$oid <- seq_along(body)
  
  sel.gr <- c(pro, body)
  sel.gr <- sel.gr[order(sel.gr$oid)]
  
  sel.gr <- split(sel.gr, seqnames(sel.gr))
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  sel.gr <- sel.gr[seqlev]
  cvg <- cvg[seqlev]
  vws <- Views(cvg, sel.gr)
  vms <- viewMeans(vws)
  sel.gr <- unlist(sel.gr)
  sel.gr$score <- unlist(vms)
  
  pro <- sel.gr[sel.gr$source %in% "promoter"]
  body <- sel.gr[sel.gr$source %in% "transcript"]
  stopifnot(identical(pro$oid, body$oid))
  sel <- txs
  sel$promoterPart <- ranges(pro)
  sel$transcriptPart <- ranges(body)
  sel$promoter <- pro$score
  sel$transcriptBody <- body$score
  smallNumber <- max(c(1e-6, min(pro$score, na.rm = TRUE), 
                       min(body$score, na.rm = TRUE)), na.rm = TRUE)
  sel$log2meanCoverage <- log2(pro$score + smallNumber) + log2(body$score + smallNumber)
  sel$PT_score <- log2(pro$score + smallNumber) - log2(body$score + smallNumber)
  sel <- sel[order(sel$PT_score, decreasing = TRUE)]
  return(sel)
}