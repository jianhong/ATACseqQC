#' Nucleosome Free Regions (NFR) score
#' @description NFR score is a raio between cut signal adjacent to TSS and that flanking 
#' the corresponding TSS. Each TSS window of 400 bp is first seperated into 3 sub-regions: 
#' the most upstream 150 bp (n1), the most downstream of 150 bp (n2), and the middle 100 bp (nf). 
#' Then the number of fragments with 5' ends overlapping each region are calculated for each TSS.
#' The NFR score for each TSS is calculated as NFR-score = log2(nf) - log2((n1+n2)/2). 
#' A plot can be generate with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis,
#' very much like a MA plot for gene expression data. 
#' @param obj an object of \link[GenomicAlignments:GAlignments-class]{GAlignments}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param nucleosomeSize numeric(1) or integer(1). Default is 150
#' @param nucleosomeFreeSize numeric(1) or integer(1). Default is 100
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @export
#' @return A object of \link[GenomicRanges:GRanges-class]{GRanges} with NFR scores
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
  smallNumber <- max(c(1e-6, min(nf$score, na.rm=TRUE), 
                       min(n1$score, na.rm=TRUE), min(n2$score, na.rm=TRUE)), na.rm = TRUE)
  sel$log2meanCoverage <- log2((3 * (n1$score + n2$score) + 2 * nf$score)/8 + smallNumber)
  sel$NFR_score <- log2(nf$score + smallNumber) + 1 - log2(n1$score + n2$score + smallNumber)
  sel <- sel[order(sel$NFR_score, decreasing = TRUE)]
  return(sel)
}