#' export list of GAlignments into bam files
#' @description wraper for \link[rtracklayer]{export} to export list of 
#' GAlignment into bam files.
#' @param objs A list of \link[GenomicAlignments]{GAlignments}.
#' @param outPath character(1). Output file path.
#' @return status of export.
#' @export
#' @importFrom rtracklayer export
#' @author Jianhong Ou
#' @examples 
#' library(GenomicAlignments)
#' gal1 <- GAlignments(seqnames=Rle("chr1"), pos=1L, cigar="10M",
#'                     strand=Rle(strand(c("+"))), names="a", score=1)
#' galist <- GAlignmentsList(a=gal1)
#' writeListOfGAlignments(galist)

writeListOfGAlignments <- function(objs, outPath="."){
    stopifnot(inherits(objs, c("list", "GAlignmentsList")))
    null <- sapply(objs, function(.ele){
        if(!is(.ele, "GAlignments")){
            stop("All elements in objs must be GAlignments.")
        }
    })
    if(is.null(outPath)){
        stop("invalid outPath.")
    }
    stopifnot(length(outPath)==1)
    if(!file.exists(outPath)){
        dir.create(outPath, showWarnings = FALSE, recursive = TRUE)
    }
    mapply(function(data, n){
        if(length(data)>0) try({
            export(data, file.path(outPath, paste0(n, ".bam")))
        })
    }, objs, names(objs))
}