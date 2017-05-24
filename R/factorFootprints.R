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
#' @importFrom stats cor.test
#' @importFrom ChIPpeakAnno featureAlignedSignal reCenterPeaks
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Biostrings matchPWM
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @return an invisible list of matrixes with the signals for plot.
#' It includes:
#'  - signal    mean values of coverage for positive strand and negative strand
#'              in feature regions
#'  - spearman.correlation    spearman correlations of cleavage counts in the 
#'              highest 10-nucleotide-window and binding prediction score.
#'  - bindingSites    predicted binding sites.
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
  stopifnot(upstream>10 && downstream>10)
  if(missing(bindingSites)){
      pwm <- motifStack::pfm2pwm(pfm)
      mt <- matchPWM(pwm, genome, min.score = "85%", with.score=TRUE,
                            exclude=names(genome)[!names(genome) %in% seqlev])
      mt.user <- matchPWM(pwm, genome, min.score=min.score, with.score=TRUE,
                     exclude=names(genome)[!names(genome) %in% seqlev])
      userdefined <- length(mt.user)
      if(length(mt)<=length(mt.user)){
          mt <- mt.user
          mt$userdefined <- TRUE
      }else{
          ol <- findOverlaps(mt, mt.user, type = "equal")
          mt$userdefined <- FALSE
          mt$userdefined[unique(queryHits(ol))] <- TRUE
      }
      if(length(mt)>10000 && userdefined<10000){## subsample 
          set.seed(seed = 1)
          mt.keep <- seq_along(mt)[!mt$userdefined]
          n <- 10000-userdefined
          if(length(mt.keep)>n){
              mt.keep <- sample(mt.keep, n, replace = FALSE)
              mt.keep <- sort(mt.keep)
              mt.keep <- seq_along(mt) %in% mt.keep
              mt <- mt[mt$userdefined | mt.keep]
          }
      }
  }else{
      stopifnot(is(bindingSites, "GRanges"))
      stopifnot(length(bindingSites)>1)
      stopifnot(length(bindingSites$score))
      mt <- bindingSites
      mt$userdefined <- TRUE
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
  mt.s <- split(mt, seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  mt.v <- Views(cvgSum, mt.s)
  mt.s <- mt.s[viewSums(mt.v)>0]
  mt <- unlist(mt.s)
  sigs <- featureAlignedSignal(cvglists=cvglist,
                              feature.gr=reCenterPeaks(mt, width=1),
                              upstream=upstream+floor(wid/2),
                              downstream=downstream+ceiling(wid/2),
                              n.tile=upstream+downstream+wid)
  cor <- lapply(sigs, function(sig){
      sig.colMeans <- colMeans(sig)
      ## calculate corelation of footprinting and binding score
      windows <- slidingWindows(IRanges(1, ncol(sig)), width = 10, step = 1)[[1]]
      # remove the windows with overlaps of motif binding region
      windows <- windows[end(windows)<=upstream | start(windows)>=upstream+wid]
      sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
      windows.sel <- windows[which.max(sig.windowMeans)][1]
      highest.sig.windows <- 
          rowMeans(sig[, start(windows.sel):end(windows.sel)])
      predictedBindingSiteScore <- mt$score
      cor <- cor.test(x = predictedBindingSiteScore, 
                      y = highest.sig.windows, 
                      method = "spearman")
  })
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  plotFootprints(colMeans(do.call(cbind, sigs)),
                 Mlen=wid, motif=pwm2pfm(pfm))
  return(invisible(list(signal=sigs, 
                        spearman.correlation=cor, 
                        bindingSites=mt)))
}

pwm2pfm <- function(pfm, name="motif"){
  if(!all(colSums(pfm)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}
