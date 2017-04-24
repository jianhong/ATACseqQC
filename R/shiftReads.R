#' @title shift read for 5'end
#' @description shift reads for 5'ends
#' @param x an object of GAlignments
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @return an object of GAlignments
#' @import S4Vectors
#' @importFrom stringr str_extract_all
#' @importFrom GenomicAlignments cigar qwidth
#' @importFrom Biostrings DNAStringSet PhredQuality
#' @author Jianhong Ou
shiftReads <- function(x, positive=4L, negative=5L){
  strds <- as.character(strand(x)) == "-"
  ns <- ifelse(strds, negative, positive)
  cigars <- str_extract_all(cigar(x), "\\d+|\\D")
  x@cigar <- mapply(function(cigar, strd, n){
    if(cigar[1]!="*"){
      cigar.num <- suppressWarnings(as.numeric(cigar))
      cigar.num[is.na(cigar.num)] <- 0
      cigar.num[which(cigar %in% c("D", "P")) - 1] <- 0
      if(strd){
        cigar.cumsum <- cumsum(rev(cigar.num))
        id.rev <- which(cigar.cumsum >= n)[1]
        id <- length(cigar) - id.rev + 1
        cigar.num[id:length(cigar.num)] <- 0
        cigar.num[id] <- cigar.cumsum[id.rev] - n
        cigar[id] <- as.character(cigar.num[id])
        cigar <- cigar[1:(length(cigar) -
                            which(rev(cigar.num)!=0)[1] +
                            2)]
      }else{
        cigar.cumsum <- cumsum(cigar.num)
        id <- which(cigar.cumsum >= n)[1]
        cigar.num[1:id] <- 0
        cigar.num[id] <- cigar.cumsum[id] - n
        cigar[id] <- as.character(cigar.num[id])
        cigar <- cigar[which(cigar.num!=0)[1]:length(cigar)]
      }
      paste(cigar, collapse="")
    }else{
      cigar
    }
  }, cigars, strds, ns)

  width.x <- qwidth(x)
  mcols(x)$isize <- sign(mcols(x)$isize) * (abs(mcols(x)$isize) - ns)
  ##x end will auto matic change if the cigar changed.
  mcols(x)$seq[strds] <- DNAStringSet(substr(mcols(x)$seq[strds], 1,
                                             width.x[strds]))
  mcols(x)$qual[strds] <- PhredQuality(substr(mcols(x)$qual[strds], 1,
                                              width.x[strds]))

  x@start[!strds] <- as.integer(start(x)[!strds] + ns[!strds])
  mcols(x)$seq[!strds] <- DNAStringSet(substr(mcols(x)$seq[!strds],
                                              ns[!strds]+1,
                                              nchar(mcols(x)$seq[!strds])))
  mcols(x)$qual[!strds] <- PhredQuality(substr(mcols(x)$qual[!strds],
                                               ns[!strds]+1,
                                               nchar(mcols(x)$qual[!strds])))
  x
}
