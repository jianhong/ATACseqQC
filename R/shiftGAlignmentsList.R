#' @title shift 5' ends
#' @description shift the GAlignmentsLists by 5' ends. 
#' All reads aligning to the positive strand will be offset by +4bp, 
#' and all reads aligning to the negative strand will be offset -5bp by default.
#' @param gal An object of \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @param outbam file path to save shift reads. If missing, no file will be write.
#' @return An object of \link[GenomicAlignments:GAlignments-class]{GAlignments} with 5' end 
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
#' gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)
#' objs <- shiftGAlignmentsList(gal)
#' export(objs, "shift.bam")
shiftGAlignmentsList <- function(gal, positive=4L, negative=5L, outbam){
    stopifnot(is.integer(positive))
    stopifnot(is.integer(negative))
    stopifnot(is(gal, "GAlignmentsList"))
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
      bamfile <- BamFile(meta$file, index=index, yieldSize=chunk, asMates = meta$asMates)
      outfile <- NULL
      mpos <- NULL
      open(bamfile)
      while (length(chunk0 <- readGAlignmentsList(bamfile, param=meta$param))) {
        gal1 <- shiftGAlignmentsList(chunk0, positive = positive, negative = negative)
        outfile <- c(tempfile(fileext = ".bam"), outfile)
        this.mpos <- mcols(gal1)$mpos
        names(this.mpos) <- paste(mcols(gal1)$qname, start(gal1))
        mpos <- c(mpos, this.mpos)
        export(gal1, outfile[1], format="BAM")
        rm(gal1)
      }
      close(bamfile)
      if(length(outfile)>1){
        mergedfile <- mergeBam(outfile, 
                               destination=tempfile(fileext = ".bam"), 
                               indexDestination=TRUE)
        unlink(outfile)
        unlink(paste0(outfile, ".bai"))
      }else{
        mergedfile <- outfile
      }
      gal1 <- readGAlignments(mergedfile, param = meta$param)
      if(missing(outbam)){
        unlink(mergedfile)
        unlink(paste0(mergedfile, ".bai"))
      }else{
        tryCatch({
          file.rename(from=mergedfile, to=outbam)
          file.rename(from=paste0(mergedfile, ".bai"), 
                      to=paste0(outbam, ".bai"))
        }, error=function(e){
          message(e)
        })
      }
      mcols(gal1)$MD <- NULL
      names(gal1) <- mcols(gal1)$qname
      gal1 <- gal1[order(names(gal1))]
      mcols(gal1)$mpos <- mpos[paste(mcols(gal1)$qname, start(gal1))]
      return(gal1)
    }
    stopifnot(length(gal)>0)
    stopifnot(all(elementNROWS(gal)<3))
    ## move the 5'end
    ## first reads is 5'end
    gal1 <- unlist(gal)
    gp <- rep(seq_along(gal), elementNROWS(gal))
    k <- width(gal1) <= max(c(positive, negative))
    if(any(k)){
      gal <- gal[-gp[k]]
      gal1 <- unlist(gal)
    }
    gp <- rep(1, length(gal1))
    gp1 <- rep(seq_along(gal), elementNROWS(gal))
    gp[duplicated(gp1)] <- 2
    mcols(gal1)$MD <- NULL
    gal1 <- shiftReads(gal1, 
                       positive=positive,
                       negative=negative)
    names(gal1) <- mcols(gal1)$qname
    mcols(gal1)$mpos[gp==2] <- start(gal1)[which(gp==2)-1]
    mcols(gal1)$mpos[gp==1] <- start(gal1)[which(gp==1)+1]
    ## till now, gal1 must have mrnm, mpos, names and flag
    stopifnot(length(mcols(gal1)$mrnm)>0)
    stopifnot(length(mcols(gal1)$mpos)>0)
    stopifnot(length(mcols(gal1)$flag)>0)
    stopifnot(length(names(gal1))>0)
    if(!(missing(outbam))){
      tryCatch(export(gal1, outbam),
               error=function(e){
                 message(e)
               })
    }
    return(gal1)
}