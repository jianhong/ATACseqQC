#' Saturation plot
#' @description plot the saturations
#' @param subsamplingPeakFiles a list of files from peak calling tools
#' @param subsamplingSizes a vector of integers, which are the sizes of subsamples for 
#' peak calling. subsamplingPeakFiles should match subsamplingSizes
#' in a correct order.
#' @param sep a string, column separator
#' @param header a boolean (TRUE or FASLE), showing whether there are headers in the 
#' peak files
#' @param fdr decimal between 0 and 1, a cutoff of statistical significance of peak detection
#' @param fdrCol integer, column index for fdr
#' @param startCol integer, column index for start positions of peak regions
#' @param endCol integer, column index for end positions of peak regions
#' @param skipLines integer, the number of lines (comment or instruction) to skip when peak 
#' files are read into R
#' @param ... parameters could be passed to plot.
#' @param peakCaller "MACS2", "ZINBA", "SPP", or "F-seq". Current only support MACS2 output
#' @importFrom utils read.delim
#' @export
#' @author Haibo Liu
#' @examples 
#' if(interactive()){
#' }


saturationPlot <- function(subsamplingPeakFiles, 
                           subsamplingSizes, 
                           peakCaller = "MACS2",
                           sep="\t", header=c(TRUE, FALSE), 
                           fdr=0.05, fdrCol, startCol=2, 
                           endCol=3, skipLines, ...){
  stopifnot(peakCaller %in% c("MACS2"))
  ## subsamplingPeakFiles positionally match subsamplingSizes
  if(!all(subsamplingPeakFiles==names(subsamplingSizes))){
    stop("Subsampling peak file names should match the name of subsamplingSize!")
  }
  
  peaks<- lapply(subsamplingPeakFiles, function(.ele){
    read.delim(file=.ele, sep=sep, header=header, as.is=TRUE, skip=skipLines)
  })
  
  names(peaks) <- subsamplingPeakFiles
  ## saturation plot based on the number of peak regions
  numPeaks <- sapply(peaks, function(.ele){
    if (peakCaller =="MACS2"){
      npeaks <- sum(10^(-.ele[fdrCol]) <= fdr)  
    }
    npeaks
  }, simplify = TRUE)
  
  names(numPeaks) <- names(peaks)
  
  plot(x= subsamplingSizes/10^6, y=numPeaks, xlab="Subsampled Fragments(M)", 
       ylab ="Significant peaks", type ="o", ...)
  peakNumber <- data.frame(x=subsamplingSizes, y=numPeaks)
  
  ## total breadth of peak region 
  breadth <- sapply(peaks, function(.ele){
    cumsum(.ele[endCol] - .ele[startCol])
  }, simplify = TRUE)
  names(breadth) <- names(peaks)
  
  plot(x= subsamplingSizes/10^6, y= breadth/10^6, xlab="Subsampled Fragments(M)", 
       ylab ="Breadth of peak regions (Mb)", type="o", ...)
  peakBreath <- data.frame(x=subsamplingSizes, y=breadth)
  return(invisible(list(peakNumber=peakNumber, peakBreath=peakBreath)))
}