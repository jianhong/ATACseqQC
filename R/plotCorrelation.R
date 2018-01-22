#' plot Correlations of multiple samples
#' @description plot PCA or heatmap for multiple bamfiles. The correlation is 
#' calculated by the counts in promoter regions.
#' @param objs an object of \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList}
#' @param txs GRanges of transcripts
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param upstream numeric(1) or integer(1). Start position of promoter. Default is 2000
#' @param downstream numeric(1) or integer(1). End position of promoter. Default is 500
#' @param type Figure type, heatmap or PCA plot.
#' @param ... parameters could be passed to downstream functions such as plot for pca or heatmap for heatmap.
#' @importClassesFrom GenomicAlignments GAlignments GAlignmentsList
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicAlignments GAlignmentsList
#' @importFrom GenomicRanges promoters coverage shift
#' @importFrom IRanges viewMeans Views
#' @importFrom stats prcomp heatmap
#' @importFrom graphics legend
#' @export
#' @return A invisible object of \link[GenomicRanges:GRanges-class]{GRanges} with counts
#' @author Jianhong Ou
#' @examples
#' library(GenomicRanges)
#' library(GenomicAlignments)
#' path <- system.file("extdata", package="ATACseqQC", mustWork=TRUE)
#' bamfiles <- dir(path, "*.bam$", full.name=TRUE)
#' gals <- lapply(bamfiles, function(bamfile){
#'                readBamFile(bamFile=bamfile, tag=character(0), 
#'                            which=GRanges("chr1", IRanges(1, 1e6)), 
#'                            asMates=FALSE)
#'         })
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' plotCorrelation(GAlignmentsList(gals), txs, seqlev="chr1")
plotCorrelation <- function(objs, txs,
                            seqlev=intersect(seqlevels(objs[[1]]), seqlevels(txs)),
                            upstream=2000, downstream=500,
                            type=c("heatmap", "PCA"), ...){
  stopifnot(is(objs, "GAlignmentsList"))
  stopifnot(length(objs)>1)
  stopifnot(is(txs, "GRanges"))
  if(length(names(objs))<1){
    names(objs) <- paste0("input", seq_along(objs))
  }
  type <- match.arg(type)
  objs <- lapply(objs, function(obj) {
    obj <- as(obj, "GRanges")
    mcols(obj) <- NULL
    obj <- promoters(obj, upstream = 0, downstream = 1)
    })
  cvgs <- lapply(objs, coverage)
  cvgs <- lapply(cvgs, function(cvg) cvg[names(cvg) %in% seqlev])
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  sel.gr <- promoters(txs, upstream = upstream, downstream = downstream)
  sel.gr <- split(sel.gr, seqnames(sel.gr))
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  cvgs <- lapply(cvgs, function(cvg) cvg[seqlev])
  sel.gr <- sel.gr[seqlev]
  vms <- lapply(cvgs, function(cvg){
    unlist(viewMeans(Views(cvg, sel.gr)))
  })
  sel.gr <- unlist(sel.gr)
  for(i in seq_along(objs)){
    mcols(sel.gr)[, names(objs)[i]] <- vms[[i]]
  }
  smallNumber <- max(c(1e-6, min(unlist(vms))))
  vms.log2 <- lapply(vms, function(.ele) log2(.ele + smallNumber))
  vms.log2 <- do.call(cbind, vms.log2)
  switch(type,
         PCA={
           vms.pca <- prcomp(vms.log2)
           dots <- list(...)
           if("col" %in% names(dots)){
             col <- dots$col
             plot(vms.pca$rotation, ...)
           }else{
             col <- seq.int(ncol(vms.log2))
             plot(vms.pca$rotation, col=col, ...)
           }
           legend("topright", legend = colnames(vms.log2), box.lwd = NA, col = col, bg=NA, pch=1)
         },
         heatmap={
           vms.cor <- cor(vms.log2, method = "spearman")
           if(all(vms.cor==1)){
             message("all samples are identical in given regions.")
           }else{
             heatmap(vms.cor, ...)
           }
         })
  return(invisible(sel.gr))
}