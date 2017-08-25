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
  flag <- scanBamFlag(isSecondaryAlignment = FALSE, 
                      isNotPassingQualityControls = FALSE)
  param <- ScanBamParam(what=c("qname", "flag", "rname",
                               "cigar", "pos", "qwidth"),
                        flag = flag)
  res <- scanBam(file, param = param)[[1L]]
  qname <- res[["qname"]]
  flag <- res[["flag"]]
  rname <- res[["rname"]]
  totalQNAMEs <- length(unique(qname))
  isDuplicate <- 
    as.logical(bamFlagAsBitMatrix(flag, "isDuplicate"))
  isMitochondria <- rname=="chrM"
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
  dupRate <- sum(isDuplicate)/length(qname)
  mitRate <- sum(isMitochondria)/length(qname)
  if(length(outPath)){
    keepQNAME <- qname[(!isMitochondria) & (!isDuplicate)]
    filter <- FilterRules(list(qn = function(x){ 
      x$qname %in% keepQNAME}))
    filterBam(file = bamfile, 
              index = index,
              destination = outPath,
              filter = filter,
              indexDestination = TRUE)
  }
  return(c(totalQNAMEs=totalQNAMEs,
           duplicateRate=dupRate,
           mitochondriaRate=mitRate,
           idxstats=idxstatsBam(file=bamfile, index=index)))
}
