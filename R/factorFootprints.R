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
#' @importFrom Biostrings matchPWM maxScore
#' @importFrom Rsamtools ScanBamParam
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
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
#' @author Jianhong Ou, Julie Zhu
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
      maxS <- maxScore(pwm)
      if(!is.numeric(min.score)){
        if(!is.character(min.score)){
          stop("'min.score' must be a single number or string")
        }else{
          nc <- nchar(min.score)
          if (substr(min.score, nc, nc) == "%"){
            min.score <- substr(min.score, 1L, nc - 1L)
          }else{
            stop("'min.score' can be given as a character string containing a percentage",
                 "(e.g. '85%') of the highest possible score")
          }
          min.score <- maxS * as.double(min.score)/100
        }
      }else{
        min.score <- min.score[1]
      }
      predefined.score <- maxS * as.double(0.85)
      mt <- matchPWM(pwm, genome, min.score = min(predefined.score, min.score), 
                     with.score=TRUE,
                     exclude=names(genome)[!names(genome) %in% seqlev])
      if (min.score <= predefined.score){
        mt$userdefined <- TRUE
      } else {
        mt$userdefined <- FALSE
        mt$userdefined[mt$score >= min.score] <- TRUE
      }
      
      if(length(mt)>10000 && sum(mt$userdefined)<10000){## subsample 
          set.seed(seed = 1)
          mt.keep <- seq_along(mt)[!mt$userdefined]
          n <- 10000-sum(mt$userdefined)
          if(length(mt.keep)>n){
              mt.keep <- mt.keep[order(mt[mt.keep]$score, decreasing = TRUE)]
              mt.keep <- mt.keep[seq.int(n)]
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
  if(sum(mt$userdefined)<2){
    stop("less than 2 binding sites by input min.score")
  }
  wid <- ncol(pfm)
  #mt <- mt[seqnames(mt) %in% seqlev]
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which)
  bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                  bamfiles, index, SIMPLIFY = FALSE)
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
  ## too much if use upstream and downstream, just use 3*wid maybe better.
  mt.s.ext <- promoters(mt.s, upstream=wid, downstream=wid+wid)
  stopifnot(all(lengths(mt.s.ext)==lengths(mt.s)))
  mt.v <- Views(cvgSum, mt.s.ext)
  mt.s <- mt.s[viewSums(mt.v)>0] 
  mt <- unlist(mt.s)
  mt.ids <- promoters(reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2)+1)
  mt.ids <- paste0(as.character(seqnames(mt.ids)), ":", start(mt.ids), "-", end(mt.ids))
  sigs <- featureAlignedSignal(cvglists=cvglist,
                              feature.gr=reCenterPeaks(mt, width=1),
                              upstream=upstream+floor(wid/2),
                              downstream=downstream+ceiling(wid/2),
                              n.tile=upstream+downstream+wid)
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  cor <- lapply(sigs, function(sig){
      sig.colMeans <- colMeans(sig)
      ## calculate correlation of footprinting and binding score
      windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, step = 1)[[1]]
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
  plotFootprints(colMeans(do.call(cbind, sigs), na.rm = TRUE),
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
