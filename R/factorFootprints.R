#' plot ATAC-seq footprints infer factor occupancy genome wide
#' @description Aggregate ATAC-seq footprint for a given motif generated
#'              over binding sites within the genome.
#' @param bamfiles A vector of characters indicates the file names of bams.
#' @param index The names of the index file of the 'BAM' file being processed;
#'        This is given without the '.bai' extension.
#' @param pfm A Position frequency Matrix represented as a numeric matrix
#'        with row names A, C, G and T.
#' @param genome An object of \link[BSgenome]{BSgenome}.
#' @param min.score The minimum score for counting a match.
#'                  Can be given as a character string containing a
#'                  percentage (e.g. "95%") of the highest possible
#'                  score or as a single number.
#'                  See \link[Biostrings]{matchPWM}.
#' @param bindingSites A object of \link[GenomicRanges]{GRanges} indicates
#' candidate binding sites (eg. the output of fimo).
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param upstream,downstream numeric(1) or integer(1).
#'        Upstream and downstream of the binding region for
#'        aggregate ATAC-seq footprint.
#' @importFrom ChIPpeakAnno featureAlignedSignal reCenterPeaks
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Biostrings matchPWM
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @return an invisible list of matrixes with the signals for plot.
#' @export
#' @references Chen, K., Xi, Y., Pan, X., Li, Z., Kaestner, K., Tyler, J.,
#' Dent, S., He, X. and Li, W., 2013.
#' DANPOS: dynamic analysis of nucleosome position and occupancy by sequencing.
#' Genome research, 23(2), pp.341-351.
#' @author Jianhong Ou
#' @examples
#'
#'shiftedBamfile <- system.file("extdata", "GL1.bam",
#'                              package="ATACseqQC")
#'library(MotifDb)
#'CTCF <- query(MotifDb, c("CTCF"))
#'CTCF <- as.list(CTCF)
#'library(BSgenome.Hsapiens.UCSC.hg19)
#'factorFootprints(shiftedBamfile, pfm=CTCF[[1]],
#'                 genome=Hsapiens,
#'                 min.score="95%", seqlev="chr1",
#'                 upstream=100, downstream=100)
#'
factorFootprints <- function(bamfiles, index=bamfiles, pfm, genome, 
                             min.score="95%", bindingSites,
                             seqlev=paste0("chr", c(1:22, "X", "Y")),
                             upstream=100, downstream=100){
  #stopifnot(length(bamfiles)==4)
  stopifnot(is(genome, "BSgenome"))
  stopifnot(all(round(colSums(pfm), digits=4)==1))
  if(missing(bindingSites)){
      pwm <- motifStack::pfm2pwm(pfm)
      mt <- matchPWM(pwm, genome, min.score=min.score, with.score=TRUE,
                     exclude=names(genome)[!names(genome) %in% seqlev])
  }else{
      stopifnot(is(bindingSites, "GRanges"))
      stopifnot(length(bindingSites)>1)
      mt <- bindingSites
  }
  wid <- ncol(pfm)
  mt <- mt[seqnames(mt) %in% seqlev]
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  ## read in bam file
  bamIn <- mapply(readGAlignments, bamfiles, index, SIMPLIFY = FALSE)
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  if(class(bamIn)!="GRangesList") bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn)
  ## keep 5'end as cutting sites
  bamIn <- promoters(bamIn, upstream=0, downstream=1)
  ## split into positive strand and negative strand
  bamIn <- split(bamIn, strand(bamIn))
  ## get coverage
  cvglist <- sapply(bamIn, coverage)
  cvglist <- cvglist[c("+", "-")]
  cvglist <- lapply(cvglist, function(.ele)
    .ele[names(.ele) %in% seqlev])
  ## coverage of mt, must be filtered, otherwise too much
  cvgSum <- cvglist[["+"]] + cvglist[["-"]]
  mt.s <- split(ranges(mt), seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  mt.s <- mapply(function(.cvg, .mt){
    v <- Views(.cvg, .mt)
    .mt[viewApply(v, sum)>0]
  }, cvgSum[seqlev], mt.s[seqlev], SIMPLIFY=FALSE)
  mt <- as(IRangesList(mt.s), "GRanges")
  sig <- featureAlignedSignal(cvglists=cvglist,
                              feature.gr=reCenterPeaks(mt, width=1),
                              upstream=upstream+floor(wid/2),
                              downstream=downstream+ceiling(wid/2),
                              n.tile=upstream+downstream+wid)
  sig <- do.call(cbind, sig)
  sig <- colMeans(sig)
  plotFootprints(sig,
                 Mlen=wid, motif=pwm2pfm(pfm))
  return(invisible(sig))
}

pwm2pfm <- function(pfm, name="motif"){
  if(!all(colSums(pfm)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}
