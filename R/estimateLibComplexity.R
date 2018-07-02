#' Library complexity estimation
#' @description Estimating the library complexity.
#' @param histFile A two-column matrix of integers. The 1st column is the frequency 
#' j = 1,2,3,.... The 2nd column is the number of genomic regions with the same 
#' fequency (j) of duplication. This file should be sorted by the first column 
#' in ascending order. For example, one row of a histogram file:
#' 10    20
#' means there are 10 genomic regions, each of which is covered by 20 identical 
#' fragments at a given sequencing depth of a sequencing library.
#' @param times An positive integer representing the minimum required number of 
#' successful estimation. Default is 100.
#' @param interpolate.sample.sizes A numeric vector with values between (0, 1].
#' @param extrapolate.sample.sizes A numeric vector with values greater than 1.
#' @importFrom preseqR ds.rSAC.bootstrap
#' @export
#' @author Haibo Liu
#' @return invisible estimates, a data frame of 3 columns: relative sequence depth, 
#' number of distinct fragments, number of putative sequenced reads.
#' @seealso \link{readsDupFreq}
#' @examples 
#' library(preseqR)
#' data(FisherButterflyHist)
#' estimateLibComplexity(histFile=FisherButterflyHist, times=100)

estimateLibComplexity <- function(histFile, times=100, 
                                  interpolate.sample.sizes=seq(0.1, 1, by=0.1),
                                  extrapolate.sample.sizes=seq(5, 20, by=5)){
    total <- histFile[,1] %*% histFile[,2]  ## added
    suppressWarnings({result = ds.rSAC.bootstrap(histFile, r=1, times=times)})
    sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
    estimates <- data.frame(relative.size=sequences, values=rep(NA, length(sequences)))
    for ( i in seq_along(sequences))
    {
        estimates$values[i] <- result$f(sequences[i])
    }
    estimates$reads <- estimates$relative.size * total ## added
    plot(x=estimates$reads/10^6, y=estimates$values/10^6,  
         type="o", xlab =expression(Putative~sequenced~fragments~x~10^6), 
         ylab=expression(Distinct~fragments~x~10^6),
         main="Estimation of ATAC-seq\nlibrary complexity")
    return(invisible(estimates))  
}
