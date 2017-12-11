#' reads duplicates frequency
#' @description Calculate the frequency of read duplication based on alignment coordinates
#' determined by rname, strand, pos, cigar, mrnm, and isize. 
#' @param bamFile character. bam file name.
#' @param index character. name of index file.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam testPairedEndBam
#' @export
#' @return a 2 column matrix with reads duplicates frequency
#' @author Haibo Liu
#' @examples 
#' bamFile <- system.file("extdata", "GL1.bam", package = "ATACseqQC")
#' freq <- readsDupFreq(bamFile)

readsDupFreq <- function(bamFile, index = bamFile){
  if(!testPairedEndBam(bamFile)){
    stop("This is not Paired-End BAM file.")
  }
  
  param <-ScanBamParam(what=c("rname", "strand", "cigar", "pos", "mrnm", "mpos", "isize"),
                       flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                        isUnmappedQuery = FALSE,
                                        isNotPassingQualityControls = FALSE))
  
  bam <- scanBam(bamFile, index=index, param=param, asMate=FALSE)
  
  lst <- lapply(names(bam[[1]]), function(elt) {
    do.call(c, unname(lapply(bam, "[[", elt)))
  })
  
  names(lst) <- names(bam[[1]])
  rm("bam")
  df <- do.call(data.frame, lst)
  rm("lst")
  df <-df[!is.na(df$isize) & df$isize >0, ]
  freqByPos <-
    as.data.frame(table(paste(df$rname, df$strand, df$cigar,
                              df$pos, df$mrnm, df$isize,
                              sep="_")))
  rm("df")
  freqByDuplication <- as.data.frame(table(freqByPos$Freq))
  freqByDuplication$Var1 <- as.numeric(as.character(freqByDuplication$Var1))
  freqByDuplication <- as.matrix(freqByDuplication)
  return(freqByDuplication)
}
