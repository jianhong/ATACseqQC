#' @title shift read for 5'end
#' @description shift reads for 5'ends
#' @param x an object of GAlignments
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @return an object of GAlignments
#' @import S4Vectors
#' @importFrom GenomicAlignments cigar qwidth sequenceLayer cigarNarrow 
#' cigarQNarrow cigarWidthAlongQuerySpace
#' @importFrom Biostrings DNAStringSet PhredQuality
#' @author Jianhong Ou
shiftReads <- function(x, positive=4L, negative=5L){
  strds <- as.character(strand(x)) == "-"
  ns <- ifelse(strds, negative, positive)
  ## fix soft-clipping
  cigars <- cigar(x)
  mcols(x)$seq <- sequenceLayer(mcols(x)$seq, cigars, 
                                from="query", 
                                to="query-after-soft-clipping")
  mcols(x)$qual <- sequenceLayer(mcols(x)$qual, cigars, 
                                 from="query", 
                                 to="query-after-soft-clipping")
  cigars <- as.character(cigarNarrow(cigars))
  cigar_width <- cigarWidthAlongQuerySpace(cigars)
  seq_width <- width(mcols(x)$seq)
  withInsertionsAt5Ends <- which(cigar_width!=seq_width)
  if(length(withInsertionsAt5Ends)>0){
    ## clip from 3'end
    mcols(x)$seq[withInsertionsAt5Ends] <- 
      substr(mcols(x)$seq[withInsertionsAt5Ends], 
             start= (seq_width+1-cigar_width)[withInsertionsAt5Ends],
             stop = seq_width[withInsertionsAt5Ends])
  }
  cigars <- cigarQNarrow(cigars, 
                         start=ifelse(strds, 1, positive+1), 
                         end=ifelse(strds, -negative-1, -1))
  x@cigar <- as.character(cigars)
  x@start <- x@start + attributes(cigars)$rshift
  
  width.x <- qwidth(x)
  if(length(mcols(x)$isize)==length(x)){
    mcols(x)$isize <- sign(mcols(x)$isize) * 
      (abs(mcols(x)$isize) - positive - negative)
  }
  ##x end will auto matic change if the cigar changed.
  mcols(x)$seq[strds] <- DNAStringSet(substr(mcols(x)$seq[strds], 1,
                                             width.x[strds]))
  mcols(x)$qual[strds] <- PhredQuality(substr(mcols(x)$qual[strds], 1,
                                              width.x[strds]))
  mcols(x)$seq[!strds] <- DNAStringSet(substr(mcols(x)$seq[!strds],
                                              ns[!strds]+1,
                                              nchar(mcols(x)$seq[!strds])))
  mcols(x)$qual[!strds] <- PhredQuality(substr(mcols(x)$qual[!strds],
                                               ns[!strds]+1,
                                               nchar(mcols(x)$qual[!strds])))
  x
}