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
  file <- BamFile(file = bamfile, index = index)
  # flag <- scanBamFlag(isSecondaryAlignment = FALSE, 
  #                     isNotPassingQualityControls = FALSE)
  flag <- scanBamFlag(isSecondaryAlignment = FALSE)
  param <- ScanBamParam(what=c("qname", "flag", "rname",
                               "cigar", "pos", "qwidth",
                               "mapq"),
                        flag = flag)
  res <- scanBam(file, param = param)[[1L]]
  qname <- res[["qname"]]
  if(all(qname=="*")){
    stop("Can not get qname.")
  }
  flag <- res[["flag"]]
  rname <- res[["rname"]]
  mapq <- res[["mapq"]]
  mapq <- as.data.frame(table(mapq))
  totalQNAMEs <- length(unique(qname))
  isDuplicate <- 
    as.logical(bamFlagAsBitMatrix(flag, "isDuplicate"))
  isMitochondria <- rname %in% mitochondria
  dupRate <- sum(isDuplicate)/length(qname)
  if(dupRate==0){
    ## need to double check duplicate rate
    if(testPairedEndBam(file)){
      ## PE, remove the fragments with same cigar and positions.
      file1 <- BamFile(file=bamfile, index=index, yieldSize=50000)
      open(file1)
      gal.names <- gal.ranges <- character(length = ceiling(length(qname)/2) + 1)
      curr.idx <- 0
      while(length(chunk0 <- readGAlignmentPairs(file1, use.names = TRUE))){
        gal <- granges(chunk0, on.discordant.seqnames="drop")
        gal.names[curr.idx+seq_along(gal)] <- names(gal)
        gal.ranges[curr.idx+seq_along(gal)] <- paste(as.character(seqnames(gal)), 
                                                     start(gal), end(gal),
                                                     strand(gal))
      }
      close(file1)
      gal.names <- gal.names[!is.na(gal.names)]
      gal.ranges <- gal.ranges[!is.na(gal.ranges)]
      duplicated.qname <- gal.names[duplicated(gal.ranges)]
      rm(gal.names)
      rm(gal.ranges)
      curr.idx <- qname %in% duplicated.qname
      ## Here, may introduce some bug, if the cigar is switched
      cigar <- do.call(rbind, 
                       split(res[["cigar"]][curr.idx],
                             qname[curr.idx])[duplicated.qname])
      duplicated.qname <- duplicated.qname[duplicated(cigar)]
      isDuplicate <- qname %in% duplicated.qname
    }else{
      ## SE, remove the reads with same cigar and positions.
      se <- do.call(cbind, 
                    res[c("rname", "cigar", "pos", "qwidth")])
      isDuplicate <- duplicated(se)
    }
  }
  lenQ <- length(qname)
  dupRate <- sum(isDuplicate)/lenQ
  mitRate <- sum(isMitochondria)/lenQ
  
  isProperPair <- as.logical(bamFlagAsBitMatrix(flag, "isProperPair"))
  isUnmappedQuery <- as.logical(bamFlagAsBitMatrix(flag, "isUnmappedQuery"))
  hasUnmappedMate <- as.logical(bamFlagAsBitMatrix(flag, "hasUnmappedMate"))
  isNotPassingQualityControls <- as.logical(bamFlagAsBitMatrix(flag, "isNotPassingQualityControls"))
  properPairRate <- sum(isProperPair)/lenQ
  unmappedRate <- sum(isUnmappedQuery)/lenQ
  hasUnmappedMateRate <- sum(hasUnmappedMate)/lenQ
  badQualityRate <- sum(isNotPassingQualityControls)/lenQ
  
  ## PCR Bottlenecking Coefficient 1 (PBC1) = M1/M_DISTINCT where
  ## M1: number of genomic locations where exactly one read maps uniquely
  ## M_DISTINCT: number of distinct genomic locations to which some read maps uniquely
  ## PCR Bottlenecking Coefficient 2 (PBC2) = M1/M2 where
  ## number of genomic locations where only one read maps uniquely
  ## number of genomic locations where two reads map uniquely
  ## Non-Redundant Fraction (NRF) = 
  ## Number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads. 
  if(testPairedEndBam(file)){ ## Could we calculate dupRate in the meanwhile?
    pe <- do.call(cbind, res[c("qname", "rname", "pos", "qwidth")])
    pe <- split(pe[, -1], pe[, 1])
    pe.len <- lengths(pe)
    pe <- split(pe, pe.len)
    pe <- lapply(pe, function(.ele){
      .ele <- do.call(rbind, .ele)
      Reduce(paste, as.data.frame(.ele))
    })
    pos <- unlist(pe, use.names = FALSE)
  }else{ ## SE
    pos <- Reduce(paste, res[c("rname", "pos", "qwidth")])
  }
  stats <- table(pos)
  dup.stats <- table(stats)
  M1 <- ifelse("1" %in% names(dup.stats), dup.stats[["1"]], 0)
  M2 <- ifelse("2" %in% names(dup.stats), dup.stats[["2"]], 1) ## avoid x/0
  M_DISTINCT <- sum(dup.stats)
  NRF <- M1/totalQNAMEs
  PBC1 <- M1/M_DISTINCT
  PBC2 <- M1/M2
  
  if(length(outPath)){
    keepQNAME <- qname[(!isMitochondria) & (!isDuplicate) & isProperPair & (!isNotPassingQualityControls)]
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
              PCRbottleneckingCoefficient_1=PBC1,
              PCRbottleneckingCoefficient_2=PBC2,
              MAPQ=mapq,
              idxstats=idxstatsBam(file=bamfile, index=index)))
}
