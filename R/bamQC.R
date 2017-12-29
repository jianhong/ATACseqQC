#' Mapping quality control
#' @description Check the mapping rate, PCR duplication rate,
#' and mitochondria reads contamination.
#' @param bamfile character(1). File name of bam.
#' @param index character(1). File name of index file.
#' @param mitochondria character(1). Sequence name of mitochondria.
#' @param outPath character(1). File name of cleaned bam.
#' @return A list of quality summary.
#' @author Jianhong Ou
#' @export
#' @importFrom S4Vectors FilterRules
#' @importFrom Rsamtools BamFile scanBamFlag ScanBamParam scanBam bamFlagAsBitMatrix filterBam idxstatsBam testPairedEndBam
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @import GenomicRanges
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' bamQC(bamfile, outPath=NULL)
bamQC <- function(bamfile, index=bamfile, mitochondria="chrM",
                  outPath=sub(".bam", ".clean.bam", 
                              basename(bamfile))){
  stopifnot(length(bamfile)==1)
  stopifnot(is.character(bamfile))
  isPE <- testPairedEndBam(file = bamfile, index = index)
  flag <- scanBamFlag(isSecondaryAlignment = FALSE)
  if(isPE){
    param <- ScanBamParam(what=c("qname", "flag", "rname",
                                 "cigar", "pos", "qwidth",
                                 "mapq", "isize"),
                          flag = flag)
  }else{
    param <- ScanBamParam(what=c("qname", "flag", "rname",
                                 "cigar", "pos", "qwidth",
                                 "mapq"),
                          flag = flag)
  }
  file <- BamFile(file = bamfile, index = index)
  res <- scanBam(file, param = param)[[1L]]
  if(all(res[["qname"]]=="*")){
    stop("Can not get qname.")
  }
  lenQ <- length(res[["qname"]])
  mapq <- as.data.frame(table(res[["mapq"]]))
  totalQNAMEs <- length(unique(res[["qname"]]))
  isMitochondria <- res[["rname"]] %in% mitochondria
  mitRate <- sum(isMitochondria)/lenQ
  
  isDuplicate <- 
    as.logical(bamFlagAsBitMatrix(res[["flag"]], "isDuplicate"))
  isProperPair <- as.logical(bamFlagAsBitMatrix(res[["flag"]], "isProperPair"))
  isUnmappedQuery <- as.logical(bamFlagAsBitMatrix(res[["flag"]], "isUnmappedQuery"))
  hasUnmappedMate <- as.logical(bamFlagAsBitMatrix(res[["flag"]], "hasUnmappedMate"))
  isNotPassingQualityControls <- as.logical(bamFlagAsBitMatrix(res[["flag"]], "isNotPassingQualityControls"))
  properPairRate <- sum(isProperPair)/lenQ
  unmappedRate <- sum(isUnmappedQuery)/lenQ
  hasUnmappedMateRate <- sum(hasUnmappedMate)/lenQ
  badQualityRate <- sum(isNotPassingQualityControls)/lenQ
  dupRate <- sum(isDuplicate)/lenQ
  res$mapq <- NULL
  res <- as.data.frame(res, stringsAsFactors=FALSE)
  #res <- res[order(res$qname), ]
  
  ## PCR Bottlenecking Coefficient 1 (PBC1) = M1/M_DISTINCT where
  ## M1: number of genomic locations where exactly one read maps uniquely
  ## M_DISTINCT: number of distinct genomic locations to which some read maps uniquely
  ## PCR Bottlenecking Coefficient 2 (PBC2) = M1/M2 where
  ## number of genomic locations where only one read maps uniquely
  ## number of genomic locations where two reads map uniquely
  ## Non-Redundant Fraction (NRF) = 
  ## Number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads. 
  seqlev <- levels(res$rname)
  if(dupRate == 0){
    isDuplicate <- logical(length=nrow(res))
  }
  pos <- NULL
  if(isPE){ ## Could we calculate dupRate in the meanwhile?
    for(i in seqlev){
      res.sub <- res[res$rname == i, ]
      isPaired <- as.logical(bamFlagAsBitMatrix(res.sub$flag, "isPaired"))
      res.se <- res.sub[!isPaired, ]
      res.pe <- res.sub[isPaired, ]
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
      if(dupRate==0){
        ## need to double check duplicate rate
        ## PE, remove the fragments with same cigar and positions.
        mates <- cbind(firstMate[, c("rname", "pos", "cigar", "isize")], 
                       secondMate[, c("rname", "pos", "cigar", "isize")])
        duplicated.qname <- firstMate[duplicated(mates), "qname"]
        rm(mates)
        duplicated.qname <- c(duplicated.qname, res.se[duplicated(res.se), "qname"])
        isDuplicate[res$rname == i] <- res.sub$qname %in% duplicated.qname
      }
      
      mates <- cbind(firstMate[, c("rname", "pos", "isize")], 
                     secondMate[, c("rname", "pos", "isize")])
      rm(firstMate)
      rm(secondMate)
      pos <- rbind(pos, mates)
      rm(mates)
      if(nrow(res.se)>0){
        res.se.cp <- res.se[, c("rname", "pos", "isize")]
        res.se.cp$rname <- NA
        res.se.cp$pos <- NA
        res.se.cp$isize <- NA
        pos <- rbind(pos, cbind(res.se[, c("rname", "pos", "isize")], res.se.cp))
        rm(res.se)
        rm(res.se.cp)
      }
    }
  }else{ ## SE
    if(dupRate==0){
      ## need to double check duplicate rate
      ## SE, remove the reads with same cigar and positions.
      isDuplicate <- duplicated(res[, c("rname", "cigar", "pos", "qwidth")])
    }
    pos <- res[, c("rname", "pos", "qwidth")]
  }
  if(dupRate == 0){
    dupRate <- sum(isDuplicate)/lenQ
  }
  qname <- res$qname
  rm(res)
  gc(reset=TRUE)
  M_DISTINCT <- nrow(unique(pos))
  pos.dup <- duplicated(pos) | duplicated(pos, fromLast = TRUE)
  M1 <- sum(!pos.dup)
  pos.dup <- pos[pos.dup, , drop=FALSE]
  pos.dup <- apply(pos.dup, 1, paste, collapse=" ")
  stats <- table(pos.dup)
  dup.stats <- table(stats)
  M2 <- ifelse("2" %in% names(dup.stats), dup.stats[["2"]], 1) ## avoid x/0

  NRF <- M1/totalQNAMEs
  PBC1 <- M1/M_DISTINCT
  PBC2 <- M1/M2
  
  if(length(outPath)){
    keepQNAME <- qname[(!isMitochondria) & (!isDuplicate) & 
                                  isProperPair & (!isNotPassingQualityControls)]
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
              MAPQ=mapq,
              idxstats=idxstatsBam(file=bamfile, index=index)))
}
