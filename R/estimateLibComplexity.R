#' library complexity estimation
#' @description estimate the library complexity
#' @param histFile A two-column matrix. The 1st column is the frequency j = 1,2,3,.... 
#' The 2nd column is the number of genomic regions with the same 
#' fequency (j) of duplication. This file should be sorted by the first column in 
#' ascending order. For example, one row of a histogram file
#' 10    20
#' means there are 10 genomic regions, each of which is covered by 20 identical 
#' fragments at a given sequencing depth of a sequencing library.
#' @param times An positive integer representing the minimum required number of successful
#' estimation. Default is 100.
#' @param interpolate.sample.sizes a numeric vector with values between (0, 1]
#' @param extrapolate.sample.sizes a numeric vector with values greater than 1
#' @importFrom preseqR ds.mincount.bootstrap
#' @export
#' @author Haibo Liu
#' @return invisible estimates
#' @seealso \link{readsDupFreq}
#' @examples 
#' library(preseqR)
#' data(FisherButterflyHist)
#' estimateLibComplexity(histFile=FisherButterflyHist, times=100)

estimateLibComplexity <- function(histFile, times=100, 
                                  interpolate.sample.sizes=seq(0.1, 1, by=0.1),
                                  extrapolate.sample.sizes=seq(5, 20, by=5)){
    result = ds.mincount.bootstrap(histFile, r=1, times=times)
    sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
    estimates <- data.frame(relative.size=sequences, values=rep(NA, length(sequences)))
    for ( i in seq_along(sequences))
    {
        estimates$values[i] <- result$FUN.bootstrap(sequences[i])
    }
    
    plot(x=estimates$relative.size, y=estimates$values/10^6,  
         type="o", xlab ="Relative sequencing depth", 
         ylab=expression(Distinct~fragments~x~10^6),
         main="Estimation of ATAC-seq\nlibrary complexity")
    return(invisible(estimates))  
}
