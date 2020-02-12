#' export list of GAlignments into bam files
#' @description wraper for \link[rtracklayer]{export} to export list of 
#' GAlignment into bam files.
#' @param objs A list of \link[GenomicAlignments:GAlignments-class]{GAlignments}.
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
        if(length(data)>0){
          ## remove the NA values
          mc <- mcols(data)
          checkpoint <- sapply(mc, function(.ele){
            any(is.na(.ele))
          })
          if(any(checkpoint)){
            for(i in which(checkpoint&nchar(colnames(mc))==2)){
              mc[, i] <- NULL
            }
            mcols(data) <- mc
          }
          try({
            export(data, file.path(outPath, paste0(n, ".bam")))
          })
        }else{
          meta <- metadata(data)
          if("file" %in% names(meta)){
            file.copy(from = meta$file, to = file.path(outPath, paste0(n, ".bam")))
            file.copy(from = paste0(meta$file, ".bai"), 
                      to = file.path(outPath, paste0(n, ".bam.bai")))
          }
        }
    }, objs, names(objs))
}