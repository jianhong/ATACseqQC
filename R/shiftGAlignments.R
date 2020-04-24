#' @title shift 5' ends for single end reads
#' @description shift the GAlignmentsLists by 5' ends. 
#' All reads aligning to the positive strand will be offset by +4bp, 
#' and all reads aligning to the negative strand will be offset -5bp by default.
#' @param gal An object of 
#' \link[GenomicAlignments:GAlignments-class]{GAlignments}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @param outbam file path to save shift reads. 
#' If missing, no file will be write.
#' @return An object of 
#' \link[GenomicAlignments:GAlignments-class]{GAlignments} with 5' end 
#' shifted reads.
#' @author Jianhong Ou
#' @export
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom Rsamtools mergeBam
#' @importFrom rtracklayer export
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' which <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
#' gal <- readBamFile(bamfile, tag=tags, 
#'                    what=c("qname", "flag", "mapq", "seq", "qual"),
#'                    which=which, asMates=FALSE, bigFile=TRUE)
#' objs <- shiftGAlignments(gal)
#' export(objs, "shift.bam")
shiftGAlignments <- function(gal, positive=4L, negative=5L, outbam){
  stopifnot(is.integer(positive))
  stopifnot(is.integer(negative))
  stopifnot(is(gal, "GAlignments"))
  if(length(gal)==0){
    ## big file mode
    meta <- metadata(gal)
    if(!all(c("file", "param") %in% names(meta))){
      stop("length of gal could not be 0.")
    }
    ow <- getOption("warn")
    on.exit(options(warn = ow))
    options(warn=-1)
    chunk <- 100000
    index <- ifelse(length(meta$index)>0, meta$index, meta$file)
    bamfile <- BamFile(meta$file, index=index, 
                       yieldSize=chunk, asMates = meta$asMates)
    outfile <- NULL
    open(bamfile)
    while (length(chunk0 <- readGAlignments(bamfile, param=meta$param))) {
      gal1 <- shiftGAlignments(chunk0, positive = positive, negative = negative)
      outfile <- c(tempfile(fileext = ".bam"), outfile)
      export(gal1, outfile[1], format="BAM")
      rm(gal1)
    }
    close(bamfile)
    if(length(outfile)>1){
      mergedfile <- mergeBam(outfile, 
                             destination=tempfile(fileext = ".bam"), 
                             indexDestination=TRUE,
                             header=meta$file)
      unlink(outfile)
      unlink(paste0(outfile, ".bai"))
    }else{
      if(length(outfile)==1){
        mergedfile <- outfile
      }else{
        stop("Can not get any proper mapped reads from  your inputs.")
      }
    }
    if(!missing(outbam)){
      file.copy(from=mergedfile, to=outbam)
      file.copy(from=paste0(mergedfile, ".bai"), 
                to=paste0(outbam, ".bai"))
      gal1 <- GAlignments()
      meta$file <- outbam
      meta$asMates <- FALSE
      metadata(gal1) <- meta
    }else{
      gal1 <- readGAlignments(mergedfile, param = meta$param)
      mcols(gal1)$MD <- NULL
      names(gal1) <- mcols(gal1)$qname
      gal1 <- gal1[order(names(gal1))]
      unlink(mergedfile)
      unlink(paste0(mergedfile, ".bai"))
    }
    return(gal1)
  }
  stopifnot(length(gal)>0)
  mcols(gal)$mrnm <- NULL
  mcols(gal)$mpos <- NULL
  gal <- shiftReads(gal, 
                    positive=positive,
                    negative=negative)
  names(gal) <- mcols(gal)$qname
  if(length(gal)<1){
    return(gal)
  }
  ## till now, gal must have names and flag
  stopifnot(length(mcols(gal)$flag)>0)
  stopifnot(length(names(gal))>0)
  mcols(gal)$mrnm <- NULL
  mcols(gal)$mpos <- NULL
  if(!(missing(outbam))){
    tryCatch(export(gal, outbam),
             error=function(e){
               message(e)
             })
  }
  return(gal)
}