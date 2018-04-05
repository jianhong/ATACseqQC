#' Calculating duplication frequency
#' @description Calculating the frequency of read duplication based on alignment
#' status determined by rname, strand, pos, cigar, mrnm, mpos and isize. 
#' @param bamFile A character vector of length 1L containing the name of a BAM file.
#' Only a BAM file with duplication reads are meaningful for estimating the library
#' complexity. For example, a raw BAM file output by aligners, or a BAM file with 
#' mitochondrial reads removed.
#' @param index A character vector of length 1L containing the name of a BAM index file.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam testPairedEndBam
#' @export
#' @return A two-column matrix of integers. The 1st column is the frequency 
#' j = 1,2,3,.... The 2nd column is the number of genomic regions with the same 
#' fequency (j) of duplication. The frequency column is in ascending order.
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
    
    lst <- lapply(names(bam[[1]]), function(.elt) {
        do.call(c, unname(lapply(bam, "[[", .elt)))
    })
    
    names(lst) <- names(bam[[1]])
    rm("bam")
    df <- do.call(data.frame, lst)
    rm("lst")
    df <-df[!is.na(df$isize) & df$isize >0, ]
    freqByPos <-
        as.data.frame(table(paste(df$rname, df$strand, df$cigar,
                                  df$pos, df$mrnm, df$isize, df$mpos,
                                  sep="_")))
    rm("df")
    freqByDuplication <- as.data.frame(table(freqByPos$Freq))
    freqByDuplication$Var1 <- as.numeric(as.character(freqByDuplication$Var1))
    freqByDuplication <- as.matrix(freqByDuplication)
    
    if (nrow(freqByDuplication) <= 3)
    {
        warning("There is not much information for estimating library complexity.\n",
             "Are you sure that you used a BAM file without remove duplication reads?", 
             call. = FALSE)
    }
    return(freqByDuplication)
}