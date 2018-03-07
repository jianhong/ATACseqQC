#' V-plot
#' @description Aggregate ATAC-seq Fragment Midpoint vs. Length for a given motif generated
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
#' @param draw Plot or not. Default TRUE.
#' @param ... parameters could be used by \link[graphics]{smoothScatter}
#' @importFrom ChIPpeakAnno reCenterPeaks
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom Biostrings matchPWM maxScore
#' @importFrom Rsamtools ScanBamParam
#' @importFrom graphics smoothScatter
#' @importFrom KernSmooth bkde2D
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @return an invisible data.frame for plot.
#' @export
#' @references Jorja G. Henikoff, Jason A. Belsky, Kristina Krassovsky, 
#' David M. MacAlpine, and Steven Henikoff.
#' Epigenome characterization at single base-pair resolution.
#' PNAS 2011 108 (45) 18318-18323
#' @author Jianhong Ou
#' @examples
#'
#'bamfile <- system.file("extdata", "GL1.bam",
#'                        package="ATACseqQC")
#'library(MotifDb)
#'CTCF <- query(MotifDb, c("CTCF"))
#'CTCF <- as.list(CTCF)
#'library(BSgenome.Hsapiens.UCSC.hg19)
#'vPlot(bamfile, pfm=CTCF[[1]],
#'      genome=Hsapiens,
#'      min.score="95%", seqlev="chr1",
#'      ylim=c(30, 250))
#'
vPlot <- function(bamfiles, index=bamfiles, pfm, genome, 
                  min.score="95%", bindingSites,
                  seqlev=paste0("chr", c(1:22, "X", "Y")),
                  upstream=200, downstream=200,
                  maxSiteNum=1e6, draw=TRUE, ...){
  stopifnot(is(genome, "BSgenome"))
  stopifnot(upstream>10 && downstream>10)
  stopifnot(is.numeric(maxSiteNum))
  maxSiteNum <- ceiling(maxSiteNum[1])
  stopifnot(maxSiteNum>1)
  if(missing(bindingSites)){
    stopifnot(all(round(colSums(pfm), digits=4)==1))
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
    wid <- ncol(pfm)
  }else{
    stopifnot(is(bindingSites, "GRanges"))
    stopifnot(all(!is.na(seqlengths(bindingSites))))
    stopifnot(length(bindingSites)>1)
    stopifnot(length(bindingSites$score)==length(bindingSites))
    mt <- bindingSites
    mt$userdefined <- TRUE
    wid <- mean(width(mt))
  }
  if(sum(mt$userdefined)<2){
    stop("less than 2 binding sites by input min.score")
  }
  #mt <- mt[seqnames(mt) %in% seqlev]
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which, 
                        flag = scanBamFlag(isProperPair=TRUE, 
                                           isSecondaryAlignment=FALSE, 
                                           isNotPassingQualityControls=FALSE, 
                                           isDuplicate=FALSE))
  
  bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                  bamfiles, index, SIMPLIFY = FALSE)
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  if(!is(bamIn, "GRangesList")) bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn, use.names = FALSE)
  seqlevelsStyle(bamIn) <- seqlevelsStyle(genome)
  mt.ext <- promoters(reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2))
  
  bamIn$w <- width(bamIn)
  bamIn <- reCenterPeaks(bamIn, width = 1)
  
  ol <- findOverlaps(bamIn, mt.ext)
  bam.query <- bamIn[queryHits(ol)]
  mt.ext.subject <- mt.ext[subjectHits(ol)]
  rel <- getRelationship(bam.query, mt.ext.subject)
  rel$FragmentLength <- bam.query$w
  rel$distanceToBindingSite <- rel$distanceToStart - upstream - floor(wid/2)
  rel <- rel[order(rel$distanceToStart, rel$FragmentLength), c("distanceToBindingSite", "FragmentLength")]
  if(draw) smoothScatter(rel, ...)
  return(invisible(rel))
}



getRelationship <- function(queryHits, subjectHits){
  if(!inherits(queryHits, "GRanges")) 
    stop("queryHits must be an object of GRanges")
  if(!inherits(subjectHits, "GRanges")) 
    stop("subjectHits must be an object of GRanges")
  strand <- strand(subjectHits)=="-"
  FeatureStart <- as.numeric(ifelse(strand, 
                                    end(subjectHits), 
                                    start(subjectHits)))
  FeatureEnd <- as.numeric(ifelse(strand, 
                                  start(subjectHits), 
                                  end(subjectHits)))
  PeakStart <- as.numeric(ifelse(strand, end(queryHits), 
                                 start(queryHits)))
  PeakEnd <- as.numeric(ifelse(strand, start(queryHits), end(queryHits)))
  ss <- PeakStart - FeatureStart
  ee <- PeakEnd - FeatureEnd
  se <- PeakStart - FeatureEnd
  es <- PeakEnd - FeatureStart
  shortestDistance <- apply(cbind(ss, ee, se, es), 1,
                            function(.ele) min(abs(.ele)))
  shortestDistanceToStart <- apply(cbind(ss, es), 1, 
                                   function(.ele) min(abs(.ele)))
  data.frame(shortestDistance=shortestDistance, 
             ss=ss,
             distanceToStart=shortestDistanceToStart)
}
