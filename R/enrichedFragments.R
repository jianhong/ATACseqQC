#' @title enrichment for nucleosome-free fragments and nucleosome signals
#' @description Get the enrichment signals for nucleosome-free fagments and
#'              nucleosomes.
#' @param bamfiles A vector of characters indicates the file names of bams.
#' @param index The names of the index file of the 'BAM' file being processed;
#'        This is given without the '.bai' extension.
#' @param gal A GAlignmentsList object or a list of GAlignmentPairs.
#'        If bamfiles is missing, gal is required.
#' @param TSS an object of \link[GenomicRanges:GRanges-class]{GRanges} indicates
#' the transcript start sites. All the width of TSS should equal to 1.
#' Otherwise, TSS will be reset to the center of input TSS.
#' @param librarySize A vector of numeric indicates the library size. Output of
#'        \link[ChIPpeakAnno]{estLibSize}
#' @param upstream,downstream numeric(1) or integer(1).
#'        Upstream and downstream size from each TSS.
#' @param n.tile numeric(1) or integer(1). The number of tiles to generate
#'        for each element of TSS.
#' @param normal.method character(1). Normalization methods,
#'        could be "none" or "quantile".
#'        See \link[limma]{normalizeBetweenArrays}.
#' @param adjustFragmentLength numeric(1) or integer(1).
#'        The size of fragment to be adjusted to.
#'        Default is set to half of the nucleosome size (80)
#' @param TSS.filter numeric(1). The filter for signal strength of each TSS.
#'        Default 0.5 indicates the average signal strength for the TSS
#'        from upstream to downstream bins should be greater than 0.5.
#' @param seqlev A vector of character indicates the sequence names to be
#'        considered.
#' @return A list of matrixes. In each matrix, each row record the signals for
#' corresponding feature.
#' @author Jianhong Ou
#' @export
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @importFrom ChIPpeakAnno featureAlignedExtendSignal reCenterPeaks
#' @importFrom limma normalizeBetweenArrays
#' @importFrom Rsamtools countBam ScanBamParam
#' @examples
#'
#' bamfiles <- system.file("extdata", "splited",
#'                        c("NucleosomeFree.bam",
#'                          "mononucleosome.bam",
#'                          "dinucleosome.bam",
#'                          "trinucleosome.bam"), package="ATACseqQC")
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' TSS <- promoters(txs, upstream=0, downstream=1)
#' library(ChIPpeakAnno)
#' librarySize <- estLibSize(bamfiles)
#' sigs <- enrichedFragments(bamfiles, TSS=TSS, librarySize=librarySize,
#'                           seqlev="chr1", TSS.filter=0)
#' sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#' featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=2020),
#'                       zeroAt=.5, n.tile=101, upper.extreme=2)
#' featureAlignedDistribution(sigs, reCenterPeaks(TSS, width=2020),
#'                            zeroAt=.5, n.tile=101, type="l")
#'
enrichedFragments <- function(bamfiles, index=bamfiles, gal, TSS, librarySize,
                              upstream = 1010L,
                              downstream = 1010L,
                              n.tile = 101L,
                              normal.method="quantile",
                              adjustFragmentLength=80L,
                              TSS.filter=0.5,
                              seqlev = paste0("chr", c(1:22, "X", "Y"))){
  if(missing(bamfiles)){
    if(missing(gal)){
      stop("gal is required if missing bamfiles")
    }else{
      if(!is(gal, "GAlignmentsList")){
        galInput <- sapply(gal, function(.ele){
          inherits(.ele, c("GAlignments", "GAlignmentPairs"))
        })
        if(any(!galInput)){
          stop("gal must be a GAlignmentsList object or a list of GAlignmentPairs.")
        }
      }
    }
    galInput <- TRUE
    stopifnot(length(gal)==4)
  }else{
    stopifnot(length(bamfiles)==4)
    stopifnot(length(bamfiles)==length(index))
    galInput <- FALSE
  }
  stopifnot(is(TSS, "GRanges"))
  #fragmentLength <- estFragmentLength(bamfiles)
  fragmentLength <- 200 #here we suppose all the fragment of nucleosome is 200.
  TSS <- TSS[seqnames(TSS) %in% seqlev]
  seqlev <- intersect(seqlev, seqlevels(TSS))
  if(length(seqlev)==0){
    stop("No intersection between seqlev and seqlevels of TSS")
  }
  seqlevels(TSS) <- seqlev
  seqinfo(TSS) <- Seqinfo(seqlev, seqlengths = seqlengths(TSS))
  TSS <- unique(TSS)
  if(TSS.filter>0){
    ## filter the TSS before count signals
    TSS.expand <- reCenterPeaks(TSS, width=upstream+downstream+1)
    if(galInput){
      cnt <- lapply(gal, function(.ele){
        countOverlaps(query=TSS.expand, 
                      subject=granges(.ele), 
                      ignore.strand=TRUE)
        })
      cnt <- do.call(cbind, cnt)
    }else{
      cnt <- lapply(bamfiles, countBam, param=ScanBamParam(which=TSS.expand))
      cnt <- do.call(cbind, lapply(cnt, function(.ele) .ele$records))
    }
    cnt <- rowMeans(cnt)>0
    TSS <- TSS[cnt]
  }
  ## Dinucleosome reads will be split into two reads,
  ## and trinucleosome reads will be split into three reads.
  if(galInput){
    sig.pe <-
      featureAlignedExtendSignal(gal = gal,
                                 feature.gr = TSS,
                                 upstream = upstream,
                                 downstream = downstream,
                                 n.tile = n.tile,
                                 fragmentLength = fragmentLength,
                                 librarySize = librarySize,
                                 adjustFragmentLength=adjustFragmentLength,
                                 pe="PE")
    sig.se <-
      featureAlignedExtendSignal(gal = gal,
                                 feature.gr = TSS,
                                 upstream = upstream,
                                 downstream = downstream,
                                 n.tile=n.tile,
                                 fragmentLength = fragmentLength,
                                 librarySize = librarySize,
                                 adjustFragmentLength=adjustFragmentLength,
                                 pe="SE")
  }else{
    sig.pe <-
      featureAlignedExtendSignal(bamfiles,
                                 index = index,
                                 feature.gr = TSS,
                                 upstream = upstream,
                                 downstream = downstream,
                                 n.tile = n.tile,
                                 fragmentLength = fragmentLength,
                                 librarySize = librarySize,
                                 adjustFragmentLength=adjustFragmentLength,
                                 pe="PE")
    sig.se <-
      featureAlignedExtendSignal(bamfiles,
                                 index = index,
                                 feature.gr = TSS,
                                 upstream = upstream,
                                 downstream = downstream,
                                 n.tile=n.tile,
                                 fragmentLength = fragmentLength,
                                 librarySize = librarySize,
                                 adjustFragmentLength=adjustFragmentLength,
                                 pe="SE")
  }
  sig.pe <- lapply(sig.pe, function(.ele){
    .ele[is.na(.ele)] <- 0
    .ele}) ## in case of NA
  sig.se <- lapply(sig.se, function(.ele){
    .ele[is.na(.ele)] <- 0
    .ele}) ## in case of NA
  sig <- sig.pe
  sig[[3]] <- sig.se[[3]]
  sig[[4]] <- sig.pe[[4]] + sig.se[[4]]
  sig[[2]] <- Reduce(`+`, sig[-1])
  sig <- sig[1:2]
  if(normal.method=="none"){
    return(sig)
  }
  ## normalization by mean of top bins values
  bins <- ceiling(fragmentLength/floor(upstream + downstream)/n.tile)
  bins <- 2*bins+1
  sig.rowMeans <- lapply(sig, function(.ele)
    apply(.ele, 1, function(.e) mean(sort(.e, decreasing = TRUE)[1:bins])))
  sig.rowMeans <- do.call(cbind, sig.rowMeans) + 1
  sig.adj.rowMeans <- normalizeBetweenArrays(sig.rowMeans,
                                             method = normal.method)
  sig.factors <- sig.adj.rowMeans / sig.rowMeans
  sig.factors <- as.list(as.data.frame(sig.factors))
  sig.1 <- mapply(function(a, b){
    a * b
  }, sig, sig.factors, SIMPLIFY = FALSE)
  sig.1[[2]] <- sig.1[[2]] - sig.1[[1]]
  sig.1[[2]][sig.1[[2]]<0] <- 0
  sig.rowMeans <- rowMeans(sig.1[[2]])
  sig.1 <- lapply(sig.1, function(.ele) .ele[sig.rowMeans>TSS.filter, ])
  sig.1
}
