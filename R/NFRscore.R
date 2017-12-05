#' Nucleosome Free Regions (NFR) score
#' @description NFR score is calculated for nuclosome signal divided by nuclosome free signal in 
#' every TSS 400bp windows. 
#' For each TSS window, it is seperated into 3 parts, first 150 bp for nuclesome 1 (n1) and 
#' last 150 bp for nucleosome 2 (n2), and the middle 100 bp is for nucleosome free (nf). 
#' The coverages of 5' ends of the fragments for each part are calculated.
#' The NFR score is calculated by NFR-score = log2(nf) - log2((n1+n2)/2)
#' @param obj an object of \link[GenomicAlignments]{GAlignments}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param nucleosomeSize numeric(1) or integer(1). Default is 150
#' @param nucleosomeFreeSize numeric(1) or integer(1). Default is 100
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @export
#' @return A object of \link[GenomicRanges]{GRanges} with NFR scores
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
#' nfr <- NFRscore(gal1, txs)
NFRscore <- function(obj, txs,
                     seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                     nucleosomeSize=150, nucleosomeFreeSize=100){
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
  sel <- promoters(txs, upstream = nucleosomeSize + floor(nucleosomeFreeSize/2),
                      downstream = nucleosomeSize + ceiling(nucleosomeFreeSize/2))
  
  n1.gr <- promoters(sel, upstream = 0, downstream = nucleosomeSize)
  n2.gr <- shift(n1.gr, shift = nucleosomeSize + nucleosomeFreeSize)
  n2.gr[strand(n2.gr)=="-"] <- shift(n1.gr[strand(n2.gr)=="-"], 
                                     shift = -1 * (nucleosomeSize + nucleosomeFreeSize))
  nf.gr <- shift(n1.gr, shift = nucleosomeSize)
  width(nf.gr) <- nucleosomeFreeSize
  nf.gr[strand(nf.gr)=="-"] <- shift(n1.gr[strand(nf.gr)=="-"], shift = -1 * nucleosomeSize)
  start(nf.gr[strand(nf.gr)=="-"]) <- end(nf.gr[strand(nf.gr)=="-"]) - nucleosomeFreeSize + 1
  
  n1.gr$source <- "n1"
  n2.gr$source <- "n2"
  nf.gr$source <- "nf"
  n1.gr$oid <- seq_along(n1.gr)
  n2.gr$oid <- seq_along(n2.gr)
  nf.gr$oid <- seq_along(nf.gr)
  
  sel.gr <- c(n1.gr, nf.gr, n2.gr)
  sel.gr <- sel.gr[order(sel.gr$oid)]
  
  sel.gr <- split(sel.gr, seqnames(sel.gr))
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  sel.gr <- sel.gr[seqlev]
  cvg <- cvg[seqlev]
  vws <- Views(cvg, sel.gr)
  vms <- viewMeans(vws)
  sel.gr <- unlist(sel.gr)
  sel.gr$score <- unlist(vms)
  n1 <- sel.gr[sel.gr$source %in% "n1"]
  n2 <- sel.gr[sel.gr$source %in% "n2"]
  nf <- sel.gr[sel.gr$source %in% "nf"]
  stopifnot(identical(n1$oid, n2$oid))
  stopifnot(identical(n1$oid, nf$oid))
  sel <- sel[n1$oid]
  sel$n1 <- n1$score
  sel$nf <- nf$score
  sel$n2 <- n2$score
  smallNumber <- min(c(1e-6, nf$score, n1$score, n2$score))
  sel$log2meanCoverage <- log2((3 * (n1$score + n2$score) + 2 * nf$score)/8 + smallNumber)
  sel$NFR_score <- log2(nf$score + smallNumber) + 1 - log2(n1$score + n2$score + smallNumber)
  sel <- sel[order(sel$NFR_score, decreasing = TRUE)]
  return(sel)
}