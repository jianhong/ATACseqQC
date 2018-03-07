#' plot ATAC-seq footprints infer factor occupancy genome wide
#' @description Aggregate ATAC-seq footprint for a given motif generated
#'              over binding sites within the genome.
#' @param bamfiles A vector of characters indicates the file names of bams.  
#' All the bamfiles will be pulled together.
#' @param index The names of the index file of the 'BAM' file being processed;
#'        This is given without the '.bai' extension.
#' @param pfm A Position frequency Matrix represented as a numeric matrix
#'        with row names A, C, G and T.
#' @param genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @param min.score The minimum score for counting a match.
#'                  Can be given as a character string containing a
#'                  percentage (e.g. "95%") of the highest possible
#'                  score or as a single number.
#'                  See \link[Biostrings]{matchPWM}.
#' @param bindingSites A object of \link[GenomicRanges:GRanges-class]{GRanges} indicates
#' candidate binding sites (eg. the output of fimo).
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param upstream,downstream numeric(1) or integer(1).
#'        Upstream and downstream of the binding region for
#'        aggregate ATAC-seq footprint.
#' @param maxSiteNum numeric(1). Maximal number of predicted binding sites.
#'        if predicted binding sites is more than this number, top maxSiteNum binding
#'        sites will be used.
#' @param anchor "cut site" or "fragment center". If "fragment center" is used, 
#'        the input bamfiles must be paired-end.
#' @importFrom stats cor.test
#' @importFrom ChIPpeakAnno featureAlignedSignal reCenterPeaks
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Biostrings matchPWM maxScore
#' @importFrom Rsamtools ScanBamParam testPairedEndBam
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
#'bamfile <- system.file("extdata", "GL1.bam",
#'                        package="ATACseqQC")
#'library(MotifDb)
#'CTCF <- query(MotifDb, c("CTCF"))
#'CTCF <- as.list(CTCF)
#'library(BSgenome.Hsapiens.UCSC.hg19)
#'factorFootprints(bamfile, pfm=CTCF[[1]],
#'                 genome=Hsapiens,
#'                 min.score="95%", seqlev="chr1",
#'                 upstream=100, downstream=100)
#'
factorFootprints <- function(bamfiles, index=bamfiles, pfm, genome, 
                             min.score="95%", bindingSites,
                             seqlev=paste0("chr", c(1:22, "X", "Y")),
                             upstream=100, downstream=100,
                             maxSiteNum=1e6, anchor="cut site"){
  #stopifnot(length(bamfiles)==4)
  stopifnot(is(genome, "BSgenome"))
  stopifnot(all(round(colSums(pfm), digits=4)==1))
  stopifnot(upstream>10 && downstream>10)
  stopifnot(is.numeric(maxSiteNum))
  maxSiteNum <- ceiling(maxSiteNum[1])
  stopifnot(maxSiteNum>1)
  anchor <- match.arg(anchor, choices = c("cut site", "fragment center"))
  if(anchor=="fragment center"){
    null <- mapply(function(.bam, .index){
      suppressMessages(pe <- testPairedEndBam(.bam, .index))
      if(!pe){
        stop("If anchor is fragment center, the bamfiles must be paired end reads.")
      }
    }, bamfiles, index)
  }
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
      suppressWarnings({
        mt <- matchPWM(pwm, genome, min.score = min(predefined.score, min.score), 
                     with.score=TRUE,
                     exclude=names(genome)[!names(genome) %in% seqlev])
      })
      if (min.score <= predefined.score){
        mt$userdefined <- TRUE
      } else {
        mt$userdefined <- FALSE
        mt$userdefined[mt$score >= min.score] <- TRUE
      }
      if(length(mt)>maxSiteNum){## subsample
        mt$oid <- seq_along(mt)
        mt <- mt[order(mt$score, decreasing = TRUE)]
        mt <- mt[seq.int(maxSiteNum)]
        mt <- mt[order(mt$oid)]
        mt$oid <- NULL
      }
      subsampleNum <- min(10000, maxSiteNum)
      if(length(mt)>subsampleNum && sum(mt$userdefined)<subsampleNum){## subsample 
          set.seed(seed = 1)
          mt.keep <- seq_along(mt)[!mt$userdefined]
          n <- subsampleNum-sum(mt$userdefined)
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
      stopifnot(all(!is.na(seqlengths(bindingSites))))
      stopifnot(length(bindingSites)>1)
      stopifnot(length(bindingSites$score)==length(bindingSites))
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
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
  }
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  if(!is(bamIn, "GRangesList")) bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn)
  seqlevelsStyle(bamIn) <- seqlevelsStyle(genome)
  if(anchor=="cut site"){
    ## keep 5'end as cutting sites
    bamIn <- promoters(bamIn, upstream=0, downstream=1)
  }else{
    ## keep fragment center
    bamIn <- reCenterPeaks(bamIn, width=1)
  }
  libSize <- length(bamIn)
  coverageSize <- sum(as.numeric(width(reduce(bamIn, ignore.strand=TRUE))))
  libFactor <- libSize / coverageSize
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
      if(length(predictedBindingSiteScore) == length(highest.sig.windows)){
        suppressWarnings({
          cor <- cor.test(x = predictedBindingSiteScore, 
                          y = highest.sig.windows, 
                          method = "spearman")
        })
      }else{
        cor <- NA
      }
      cor
  })
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  ## segmentation the signals
  ## x2 because stranded.
  Profile <- lapply(sigs, function(.ele) colMeans(.ele, na.rm = TRUE)*2/libFactor)
  ## upstream + wid + downstream
  Profile.split <- lapply(Profile, function(.ele){
    list(upstream=.ele[seq.int(upstream)],
         binding=.ele[upstream+seq.int(wid)],
         downstream=.ele[upstream+wid+seq.int(downstream)])
  })
  optimalSegmentation <- function(.ele){
    .l <- length(.ele)
    short_abun <- cumsum(.ele)/seq.int(.l)
    long_abun <- cumsum(rev(.ele))/seq.int(.l)
    long_abun <- rev(long_abun)
    short_abun <- short_abun[-length(short_abun)]
    long_abun <- long_abun[-1]
    ##long_abun should always greater than short_abun
    long_abun <- long_abun - short_abun
    long_abun[long_abun<0] <- 0
    cov_diff <- numeric(length(short_abun))
    for(i in seq_along(.ele)){
      cov_diff_tmp <- .ele
      cov_diff_tmp <- cov_diff_tmp-short_abun[i]
      cov_diff_tmp[-seq.int(i)] <- cov_diff_tmp[-seq.int(i)] - long_abun[i]
      cov_diff[i] <- mean(cov_diff_tmp^2)
    }
    .ids <- which(cov_diff==min(cov_diff, na.rm = TRUE))
    data.frame(pos=.ids, short_abun=short_abun[.ids], long_abun=long_abun[.ids])
  }
  Profile.seg <- lapply(Profile.split, function(.ele){
    ups <- optimalSegmentation(.ele$upstream)
    downs <- optimalSegmentation(rev(.ele$downstream))
    ## find the nearest pair
    .min <- c(max(rbind(ups, downs)), 0, 0)
    for(i in seq.int(nrow(ups))){
      for(j in seq.int(nrow(downs))){
        tmp <- sum(abs(ups[i, -1] - downs[j, -1]))
        if(tmp < .min[1]){
          .min <- c(tmp, i, j)
        }
      }
    }
    c(colMeans(rbind(ups[.min[2], ], downs[.min[3], ])), binding=mean(.ele$binding, na.rm=TRUE))
  })
  Profile.seg <- colMeans(do.call(rbind, Profile.seg))
  Profile.seg[3] <- Profile.seg[2]+Profile.seg[3]
  names(Profile.seg)[2:3] <- c("distal_abun", "proximal_abun")
  tryCatch({ ## try to avoid the error when ploting.
    ylab <- ifelse(anchor=="cut site", "Cut-site probability", "reads density (arbitrary unit)")
    plotFootprints(c(Profile[["+"]], Profile[["-"]]), 
                   Mlen=wid, motif=pwm2pfm(pfm), ylab=ylab,
                   segmentation=Profile.seg)
  }, error=function(e){
    message(e)
  })
  return(invisible(list(signal=sigs, 
                        spearman.correlation=cor, 
                        bindingSites=mt,
                        Mlen=wid,
                        estLibSize=libSize,
                        Profile.segmentation=Profile.seg)))
}

pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}
