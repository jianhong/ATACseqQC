#' Mapping quality control
#' @description Check the mapping rate, PCR duplication rate,
#' and mitochondria reads contamination.
#' @param bamfile character(1). File name of bam.
#' @param index character(1). File name of index file.
#' @param mitochondria character(1). Sequence name of mitochondria.
#' @param outPath character(1). File name of cleaned bam.
#' @param doubleCheckDup logical(1). Double check duplicates or not if there is no tags for that.
#' @return A list of quality summary.
#' @author Jianhong Ou
#' @export
#' @importFrom S4Vectors FilterRules
#' @importFrom Rsamtools BamFile scanBamFlag ScanBamParam scanBam bamFlagAsBitMatrix filterBam idxstatsBam testPairedEndBam scanBamHeader
#' @import GenomicRanges
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' bamQC(bamfile, outPath=NULL)
bamQC <- function(bamfile, index=bamfile, mitochondria="chrM",
                  outPath=sub(".bam", ".clean.bam", 
                              basename(bamfile)),
                  doubleCheckDup=FALSE){
  stopifnot(length(bamfile)==1)
  stopifnot(is.character(bamfile))
  isPE <- testPairedEndBam(file = bamfile, index = index)
  flag <- scanBamFlag(isSecondaryAlignment = FALSE)
  file <- BamFile(file = bamfile, index = index)
  file.targets <- scanBamHeader(files = file, what="targets")$targets
  if(length(file.targets)<1) return(NA)
  file.targets <- GRanges(names(file.targets), IRanges(rep(1, length(file.targets)), file.targets))
  qname <- NULL
  M1 <- NULL
  M_DISTINCT <- NULL
  M2 <- NULL
  isMitochondria <- Rle()
  isDuplicate <- Rle()
  isProperPair <- Rle()
  isUnmappedQuery <- Rle()
  hasUnmappedMate <- Rle()
  isNotPassingQualityControls <- Rle()
  lenQs <- NULL
  mapqs <- NULL
  for(i in seq_along(file.targets)){
    if(isPE){
      param <- ScanBamParam(what=c("qname", "flag", 
                                   "cigar", "pos", "qwidth",
                                   "mapq", "isize"),
                            flag = flag,
                            which = file.targets[i])
    }else{
      param <- ScanBamParam(what=c("qname", "flag", 
                                   "cigar", "pos", "qwidth",
                                   "mapq"),
                            flag = flag,
                            which = file.targets[i])
    }
    res <- scanBam(file, param = param)[[1L]]
    if(length(res[["qname"]])==0){
      next()
    }
    if(all(res[["qname"]]=="*")){
      next()
    }
    lenQ <- length(res[["qname"]])
    mapq <- as.data.frame(table(res[["mapq"]]))
    if(length(mapqs)>0){
      mapqs <- rbind(mapqs, mapq)
      mapqs <- rowsum(mapqs$Freq, group = mapqs$Var1)
      mapqs <- data.frame(Var1=row.names(mapqs), Freq=mapqs[, 1])
    }else{
      mapqs <- mapq
    }
    qname <- c(qname, res[["qname"]])
    isMitochondria <- c(isMitochondria, rep(as.character(seqnames(file.targets))[i] %in% mitochondria, lenQ))
    isDup <- 
      as.logical(bamFlagAsBitMatrix(res[["flag"]], "isDuplicate"))
    isProperPair <- c(isProperPair, as.logical(bamFlagAsBitMatrix(res[["flag"]], "isProperPair")))
    isUnmappedQuery <- c(isUnmappedQuery, as.logical(bamFlagAsBitMatrix(res[["flag"]], "isUnmappedQuery")))
    hasUnmappedMate <- c(hasUnmappedMate, as.logical(bamFlagAsBitMatrix(res[["flag"]], "hasUnmappedMate")))
    isNotPassingQualityControls <- c(isNotPassingQualityControls, 
                                     as.logical(bamFlagAsBitMatrix(res[["flag"]], "isNotPassingQualityControls")))
    res$mapq <- NULL
    res <- as.data.frame(res, stringsAsFactors=FALSE)
    
    ## PCR Bottlenecking Coefficient 1 (PBC1) = M1/M_DISTINCT where
    ## M1: number of genomic locations where exactly one read maps uniquely
    ## M_DISTINCT: number of distinct genomic locations to which some read maps uniquely
    ## PCR Bottlenecking Coefficient 2 (PBC2) = M1/M2 where
    ## number of genomic locations where only one read maps uniquely
    ## number of genomic locations where two reads map uniquely
    ## Non-Redundant Fraction (NRF) = 
    ## Number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads. 
    dupRate <- sum(isDup)/lenQ
    if(isPE){
      isPaired <- as.logical(bamFlagAsBitMatrix(res$flag, "isPaired"))
      res.se <- res[!isPaired, ]
      res.pe <- res[isPaired, ]
      isFirstMateRead <- as.logical(bamFlagAsBitMatrix(res.pe$flag, "isFirstMateRead"))
      isSecondMateRead <- as.logical(bamFlagAsBitMatrix(res.pe$flag, "isSecondMateRead"))
      firstMate <- res.pe[isFirstMateRead, ]
      secondMate <- res.pe[isSecondMateRead, ]
      rm(res.pe)
      mates.qname <- intersect(firstMate$qname, secondMate$qname)
      res.se <- rbind(res.se, 
                      firstMate[!firstMate$qname %in% mates.qname, ],
                      secondMate[!secondMate$qname %in% mates.qname, ])
      firstMate <- firstMate[match(mates.qname, firstMate$qname), ]
      secondMate <- secondMate[match(mates.qname, secondMate$qname), ]
      if(dupRate==0 & doubleCheckDup){
        ## need to double check duplicate rate
        ## PE, remove the fragments with same cigar and positions.
        mates <- cbind(firstMate[, c("pos", "cigar", "isize")], 
                       secondMate[, c("pos", "cigar", "isize")])
        duplicated.qname <- firstMate[duplicated(mates), "qname"]
        rm(mates)
        duplicated.qname <- c(duplicated.qname, res.se[duplicated(res.se), "qname"])
        isDup <- res$qname %in% duplicated.qname
      }
      
      pos <- cbind(firstMate[, c("pos", "isize")], 
                     secondMate[, c("pos", "isize")])
      rm(firstMate)
      rm(secondMate)
      if(nrow(res.se)>0){
        res.se.cp <- res.se[, c("pos", "isize")]
        res.se.cp$pos <- NA
        res.se.cp$isize <- NA
        pos <- rbind(pos, cbind(res.se[, c("pos", "isize")], res.se.cp))
        rm(res.se)
        rm(res.se.cp)
      }
    }else{ ## SE
      if(dupRate==0 & doubleCheckDup){
        ## need to double check duplicate rate
        ## SE, remove the reads with same cigar and positions.
        isDup <- duplicated(res[, c("cigar", "pos", "qwidth")])
      }
      pos <- res[, c("pos", "qwidth")]
    }
    if(dupRate == 0 & doubleCheckDup){
      dupRate <- sum(isDup)/lenQ
    }
    rm(res)
    gc(reset=TRUE)
    M_DISTINCT <- c(M_DISTINCT, nrow(unique(pos)))
    pos.dup <- duplicated(pos) | duplicated(pos, fromLast = TRUE)
    M1 <- c(M1, sum(!pos.dup))
    pos.dup <- pos[pos.dup, , drop=FALSE]
    pos.dup <- as.data.frame(pos.dup)
    pos.dup <- do.call(paste, as.list(pos.dup))
    stats <- table(pos.dup)
    dup.stats <- table(stats)
    M2 <- c(M2, ifelse("2" %in% names(dup.stats), dup.stats[["2"]], 0))
    isDuplicate <- c(isDuplicate, isDup)
    lenQs <- c(lenQs, lenQ)
  }
  
  totalQNAMEs <- length(unique(qname))
  M1 <- sum(M1)
  NRF <- M1/totalQNAMEs
  PBC1 <- M1/sum(M_DISTINCT)
  M2 <- max(c(1, sum(M2))) ## avoid x/0
  PBC2 <- M1/M2
  
  lenQs <- sum(lenQs)
  mitRate <- sum(isMitochondria)/lenQs
  dupRate <- sum(isDuplicate)/lenQs
  properPairRate <- sum(isProperPair)/lenQs
  unmappedRate <- sum(isUnmappedQuery)/lenQs
  hasUnmappedMateRate <- sum(hasUnmappedMate)/lenQs
  badQualityRate <- sum(isNotPassingQualityControls)/lenQs
  
  if(length(outPath)){
    keepQNAME <- qname[(!as.logical(isMitochondria)) & (!as.logical(isDuplicate)) & 
                         as.logical(isProperPair) & (!as.logical(isNotPassingQualityControls))]
    filter <- FilterRules(list(qn = function(x){ 
      x$qname %in% keepQNAME}))
    filterBam(file = bamfile, 
              index = index,
              destination = outPath,
              filter = filter,
              indexDestination = TRUE)
  }
  return(list(totalQNAMEs=totalQNAMEs,
              duplicateRate=dupRate,
              mitochondriaRate=mitRate,
              properPairRate=properPairRate,
              unmappedRate=unmappedRate,
              hasUnmappedMateRate=hasUnmappedMateRate,
              notPassingQualityControlsRate=badQualityRate,
              nonRedundantFraction=NRF,
              PCRbottleneckCoefficient_1=PBC1,
              PCRbottleneckCoefficient_2=PBC2,
              MAPQ=mapqs,
              idxstats=idxstatsBam(file=bamfile, index=index)))
}
