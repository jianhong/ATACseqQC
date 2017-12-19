#' Plotting Saturation curves
#' @description Plotting the saturation curves.
#' @param subsamplingPeakFiles A character vector containing peak files from peak
#' calling tools, such as MACS2. Currently only MACS2 output is supported.
#' @param subsamplingSizes A named vector of integers, which are the sizes of 
#' subsamples for peak calling. The names of subsamplingPeakFiles should be identical
#' to the basenames of subsamplingPeakFiles.
#' @param sep A character vector of length 1L, which is the column separator used 
#' in peak files.
#' @param header A boolean (TRUE or FASLE) vector of length 1L, showing whether 
#' there are column headers in the peak files.
#' @param fdr A decimal between 0 and 1, a cutoff of statistical significance of 
#' peak detection.
#' @param fdrCol An integer, column index for fdr.
#' @param startCol An integer, column index for start positions of peak regions.
#' @param endCol An integer, column index for end positions of peak regions.
#' @param skipLines An integer, the number of lines (comments or instruction) to 
#' skip when peak files are read into R.
#' @param peakCaller A character vector of length 1L contatining the name of the 
#' peak caller used to generate the peak files, such as "MACS2". Currently only 
#' MACS2 output (XXX.narrowPeak or XXX.broadPeak) is support.
#' @param outPrefix A character vector of length 1L, the file prefix for outputting
#' saturation plots.
#' @param span An integer, the span parameter for loess smoothing to fit a smoothed 
#' saturation curve.
#' @param degree An integer, the degree of local polynomial used for loess.
#' @importFrom utils read.delim
#' @importFrom stats loess.smooth
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics lines
#' @return  A data frame of three columns:
#' subsamplingSizes, the number of subsampled fragments;
#' numPeaks, the number of peaks with fdr less than a given threshold when a given
#' number of fragmetns are subsampled;
#' breadth, the total breadth of peaks with fdr less than a given threshold for 
#' give subsampling when a given number of fragmetns are subsampled.
#' @export
#' @author Haibo Liu
#' @examples 
#' if(interactive()){
#' }


saturationPlot <- function(subsamplingPeakFiles, subsamplingSizes, sep="\t", 
                           header=FALSE, fdr=0.05, fdrCol=9, startCol=2, 
                           endCol=3, skipLines=1, peakCaller="MACS2", outPrefix, span=2, degree=2)
{
    ## make sure subsamplingPeakFiles positionally match subsamplingSizes
    if(is.null(names(subsamplingSizes)) || ! identical(subsamplingPeakFiles,names(subsamplingSizes)))
    {
        stop("Subsampling peak file names should match the name of subsamplingSize!")
    }
    
    peaks<- lapply(subsamplingPeakFiles, function(.ele){
        read.delim(file=.ele, sep=sep, header=header, as.is=TRUE, skip=skipLines)
    })
    
    names(peaks) <- subsamplingPeakFiles
    
    ## saturation plotting based on the number of peak regions
    numPeaks <- sapply(peaks, function(.ele){
        if (peakCaller =="MACS2")
        {
            npeaks <- sum(10^(-.ele[fdrCol]) <= fdr)  
        }
        npeaks
    }, simplify = TRUE)
    names(numPeaks) <- names(peaks)
    
    
    ## saturation plotting based on the total breadth of peak region 
    breadth <- sapply(peaks, function(.ele){
        sum(.ele[endCol] - .ele[startCol] + 1)
    }, simplify = TRUE)
    names(breadth) <- names(peaks)
    
    peakStat <- data.frame(subsamplingSizes=subsamplingSizes, numPeaks =numPeaks, breadth=breadth)
    
    ## changing scales of x and y for plotting
    subsamplingSizes <- subsamplingSizes/10^6
    numPeaks <- numPeaks/10^3
    breadth <- breadth/10^6
    
    ## plotting saturation curves based on the number of peaks and the overall peak breadth
    pdf(paste0(outPrefix, " peak number-based saturation plot.pdf"), width=5, height=5)
    plot(x= subsamplingSizes, y= numPeaks, pch=16, xlab=expression(Effective~fragments~x~10^6),
         ylab=expression(Peaks~x~10^3))
    lines(loess.smooth(x=subsamplingSizes, y=numPeaks, span = span, degree = degree,
                       family =  "gaussian", evaluation = 50))
    dev.off()
    
    pdf(paste0(outPrefix, " peak breadth-based saturation plot.pdf"), width=5, height=5)
    plot(x= subsamplingSizes, y= breadth, pch=16, xlab=expression(Effective~fragments~x~10^6),
         ylab=expression(Total~peak~breadth~(Mb)))
    lines(loess.smooth(x=subsamplingSizes, y=breadth, span = span, degree = degree,
                       family =  "gaussian", evaluation = 50))
    dev.off()
    
    return(invisible(peakStat))
}


