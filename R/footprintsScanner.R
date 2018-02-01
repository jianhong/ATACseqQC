#' scan ATAC-seq footprints infer factor occupancy genome wide
#' @description Aggregate ATAC-seq footprint for a bunch of motifs generated
#'              over binding sites within the genome.
#' @param bamfiles A vector of characters indicates the file names of bams. 
#' All the bamfiles will be pulled together.
#' @param index The names of the index file of the 'BAM' file being processed;
#'        This is given without the '.bai' extension.
#' @param bindingSitesList A object of \link[GenomicRanges:GRangesList-class]{GRangesList} indicates
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
#' @author Jianhong Ou
#' @examples
#'
#'bamfile <- system.file("extdata", "GL1.bam",
#'                        package="ATACseqQC")
#'bsl <- system.file("extdata", "jolma2013.motifs.bindingList.95.rds",
#'                   package="ATACseqQC")
#'bindingSitesList <- readRDS(bsl)
#'footprintsScanner(bamfile, bindingSitesList=bindingSitesList)
#'
footprintsScanner <- function(bamfiles, index=bamfiles, bindingSitesList, 
                                 seqlev="chr1", upstream=100, downstream=100){
  stopifnot(upstream>10 && downstream>10)
  if(missing(bindingSitesList)){
    stop("bindingSitesList is required.")
  }else{
    mts <- lapply(bindingSitesList, function(bindingSites){
      if(!is(bindingSites, "GRanges")){
        stop("bindingSitesList must be a object of GRangesList")
      }
      if(length(bindingSites)<=1){
        stop("The length of elements in bindingSitesList must be greater than 1.")
      }
      bindingSites[seqnames(bindingSites) %in% seqlev]
    })
  }
  mts <- mts[lengths(mts)>100]
  if(length(mts)<1){
    stop("no enough bindingsites in given seqlev.")
  }
  mts.unlist <- unlist(GRangesList(mts), use.names = FALSE)
  mts.unlist$motif <- rep(names(mts), lengths(mts))
  seqlevels(mts.unlist) <- seqlev
  seqinfo(mts.unlist) <- Seqinfo(seqlev, seqlengths = seqlengths(mts.unlist))
  wid <- max(width(mts.unlist))
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mts.unlist), "GRanges")
  param <- ScanBamParam(which=which)
  bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                  bamfiles, index, SIMPLIFY = FALSE)
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  if(!is(bamIn, "GRangesList")) bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn)
  seqlevelsStyle(mts.unlist) <- seqlevelsStyle(bamIn)
  ## keep 5'end as cutting sites
  bamIn <- promoters(bamIn, upstream=0, downstream=1)
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
  mt.s <- split(mts.unlist, seqnames(mts.unlist))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  ## remove the empty regions
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
  ## split the signals
  ## segmentation the signals
  ## x2 because stranded.
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  Profiles <- lapply(sigs, function(.ele){
    aggregate(as.data.frame(.ele), 
              by=list(motif=mt$motif), 
              FUN=function(x) mean(x, na.rm=TRUE)*2/libFactor,
              simplify = FALSE)
  })
  mt.wid <- data.frame(w=width(mt), n=mt$motif)
  mt.wid <- unique(mt.wid)
  mt.width <- mt.wid$w
  names(mt.width) <- as.character(mt.wid$n)
  Profile.splits <- sapply(names(mt.width), function(.name){
    lapply(Profiles, function(Profile){
      Profile <- Profile[Profile[, 1]==.name, -1]
      Profile <- as.numeric(as.matrix(Profile))
      this.wid <- wid - mt.width[.name]
      list(upstream = Profile[ceiling(this.wid/2)+seq.int(upstream)],
           binding = Profile[ceiling(this.wid/2)+upstream+seq.int(mt.width[.name])],
           downstream =Profile[ceiling(this.wid/2)+upstream+mt.width[.name]+seq.int(downstream)])
    })
  },simplify = FALSE)
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
  Profiles.segmentation <- lapply(Profile.splits, function(Profile.split){
    Profile.seg <- lapply(Profile.split, function(.ele){
      ups <- optimalSegmentation(.ele$upstream)
      downs <- optimalSegmentation(rev(.ele$downstream))
      ## find the nearest pair
      tmp <- ups[rep(seq.int(nrow(ups)), each=nrow(downs)), -1] - 
        downs[rep(seq.int(nrow(downs)), nrow(ups)), -1]
      tmp <- rowSums(abs(tmp))
      wm <- which.min(tmp)
      i <- ceiling(wm / nrow(ups))
      j <- wm - (i-1)*nrow(ups)
      .min <- c(tmp[wm], i, j)
      c(colMeans(rbind(ups[.min[2], ], downs[.min[3], ])), 
        binding=mean(.ele$binding, na.rm=TRUE))
    })
    Profile.seg <- colMeans(do.call(rbind, Profile.seg))
    Profile.seg[3] <- Profile.seg[2]+Profile.seg[3]
    names(Profile.seg)[2:3] <- c("distal_abun", "proximal_abun")
    Profile.seg
  })
  Profiles.segmentation <- do.call(rbind, Profiles.segmentation)
  distance <- abs(Profiles.segmentation[, "proximal_abun"] - 
                    Profiles.segmentation[, "binding"])
  Profiles.segmentation[order(distance, decreasing = TRUE), ]
}


#' helper function for preparing the binding list
#' @param pfms A list of Position frequency Matrix represented as a numeric matrix
#'        with row names A, C, G and T.
#' @param genome An object of \link[BSgenome:BSgenome-class]{BSgenome}.
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param expSiteNum numeric(1). Expect number of predicted binding sites.
#'        if predicted binding sites is more than this number, top expSiteNum binding
#'        sites will be used.
#' @importFrom Biostrings matchPWM maxScore
#' @examples
#' library(MotifDb)
#' motifs <- query(MotifDb, c("Hsapiens"))
#' motifs <- as.list(motifs)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' #bindingSitesList <- prepareBindingSitesList(motifs, genome=Hsapiens)
prepareBindingSitesList <- function(pfms, genome, 
                                    seqlev=paste0("chr", c(1:22, "X", "Y")),
                                    expSiteNum=5000){
  mts <- lapply(pfms, function(pfm){
    if(!all(round(colSums(pfm), digits=4)==1)){
      stop("pfms must be list of Position frequency Matrix")
    }
    pwm <- motifStack::pfm2pwm(pfm)
    min.score <- 95
    suppressWarnings({
      mt <- matchPWM(pwm, genome, min.score = paste0(min.score, "%"),
                     with.score=TRUE, exclude=names(genome)[!names(genome) %in% seqlev])
    })
    while(length(mt)<expSiteNum && min.score>80){
      min.score <- min.score - 5
      suppressWarnings({
        mt <- matchPWM(pwm, genome, min.score = paste0(min.score, "%"),
                       with.score=TRUE, exclude=names(genome)[!names(genome) %in% seqlev])
      })
    }
    if(length(mt)>expSiteNum){## subsample
      mt$oid <- seq_along(mt)
      mt <- mt[order(mt$score, decreasing = TRUE)]
      mt <- mt[seq.int(expSiteNum)]
      mt <- mt[order(mt$oid)]
      mt$oid <- NULL
    }
    mt$string <- NULL
    mt
  })
}
