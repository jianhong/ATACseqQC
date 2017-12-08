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
      gal <- readGAlignmentPairs(file, use.names = TRUE, 
                                 param=param)
      gal <- granges(gal, on.discordant.seqnames="drop")
      ## Here, may introduce some bug, if the cigar is switched
      cigar <- do.call(rbind, 
                       split(res[["cigar"]], qname)[names(gal)])
      duplicated.qname <- 
        names(gal)[duplicated(gal) & duplicated(cigar)]
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
              MAPQ=mapq,
              idxstats=idxstatsBam(file=bamfile, index=index)))
}
