#' @title shift 5' ends
#' @description shift the GAlignmentsLists by 5' ends. 
#' All reads aligning to the positive strand will be offset by +4bp, 
#' and all reads aligning to the negative strand will be offset -5bp by default.
#' @param gal An object of \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @param outbam file path to save shift reads. If missing, no file will be write.
#' @return An object of \link[GenomicAlignments:GAlignments-class]{GAlignments} with 5' end 
#' shifted reads. The PCR duplicated will be removed unless there is metadata
#' keepDuplicates set to TRUE.
#' @author Jianhong Ou
#' @export
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom Rsamtools mergeBam
#' @importFrom rtracklayer export
#' @importFrom utils read.csv write.csv
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
      on.exit({options(warn = ow);close(bamfile)})
      df4Duplicates <- NULL
      while (length(chunk0 <- readGAlignmentsList(bamfile, param=meta$param))) {
        metadata(chunk0)$df4Duplicates <- df4Duplicates
        gal1 <- shiftGAlignmentsList(chunk0, positive = positive, 
                                     negative = negative)
        this.mpos <- mcols(gal1)$mpos
        if(length(this.mpos)!=length(gal1)){
          stop("Can not get mpos info from the reads.")
        }
        names(this.mpos) <- paste(mcols(gal1)$qname, start(gal1))
        mpos <- c(mpos, this.mpos)
        if(length(metadata(gal1)$df4Duplicates)){
          df4Duplicates <- read.csv(metadata(gal1)$df4Duplicates)
        }
        outfile <- c(tempfile(fileext = ".bam"), outfile)
        if(length(meta$header)>0) metadata(gal1)$header <- meta$header
        exportBamFile(gal1, outfile[1])
        rm(gal1)
      }
      close(bamfile)
      on.exit()
      if(length(outfile)>1){
        BAI <- paste0(outfile[1], ".bai")
        mergedfile <- mergeBam(outfile, 
                               destination=tempfile(fileext = ".bam"), 
                               indexDestination=file.exists(BAI),
                               header=meta$file)
        unlink(outfile)
        if(file.exists(BAI)) unlink(paste0(outfile, ".bai"))
      }else{
        if(length(outfile)==1){
          mergedfile <- outfile
        }else{
          stop("Can not get any proper mapped reads from  your inputs.")
        }
      }
      if(!missing(outbam)){
        file.copy(from=mergedfile, to=outbam)
        gal1 <- GAlignments()
        meta$file <- outbam
        if(file.exists(paste0(mergedfile[1], ".bai"))){
          file.copy(from=paste0(mergedfile, ".bai"), 
                    to=paste0(outbam, ".bai"))
          meta$index <- outbam
        }else{
          meta$index <- NULL
        }
        meta$asMates <- FALSE
        meta$mpos <- mpos
        metadata(gal1) <- meta
      }else{
        gal1 <- readGAlignments(mergedfile, param = meta$param)
        mcols(gal1)$MD <- NULL
        names(gal1) <- mcols(gal1)$qname
        gal1 <- gal1[order(names(gal1))]
        mcols(gal1)$mpos <- mpos[paste(mcols(gal1)$qname, start(gal1))]
        unlink(mergedfile)
        if(file.exists(paste0(mergedfile, ".bai"))){
          unlink(paste0(mergedfile, ".bai"))
        }
      }
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
    id2 <- which(gp==2)
    id1 <- id2-1
    checkDup <- FALSE
    if(length(metadata(gal)$keepDuplicates)){
      if(length(mcols(gal1)$isize)>0 && !metadata(gal)$keepDuplicates){
        checkDup <- TRUE
      }
    }else{
      checkDup <- TRUE
    }
    if(checkDup){
      if(length(metadata(gal)$df4Duplicates)){
        df0 <- metadata(gal)$df4Duplicates
      }else{
        df0 <- NULL
      }
      df <- DataFrame(seqnames=seqnames(gal1), 
                      cigar=cigar(gal1), 
                      start=start(gal1), 
                      isize=mcols(gal1)$isize)
      df <- rbind(df0, df)
      dupids <- duplicated(df)
      df2 <- dupids[id2+nrow(df0)]
      df1 <- dupids[id1+nrow(df0)]
      if(any(df1 & df2)){
        tobeRemoved <- sort(c(id2[df1 & df2], id1[df1 & df2]))
        gal1 <- gal1[-tobeRemoved]
        gp <- gp[-tobeRemoved]
        id2 <- which(gp==2)
        id1 <- id2-1
        df <- df[-(tobeRemoved+nrow(df0)), , drop=FALSE]
        if(nrow(df)){
          df4Duplicates <- tempfile()
          write.csv(df, df4Duplicates, row.names = FALSE)
          rm(df)
          metadata(gal1)$df4Duplicates <- df4Duplicates
        }
      }
    }
    gal1 <- shiftReads(gal1, 
                       positive=positive,
                       negative=negative)
    names(gal1) <- mcols(gal1)$qname
    mcols(gal1)$mpos[gp==2] <- start(gal1)[id1]
    mcols(gal1)$mpos[id1] <- start(gal1)[id2]
    if(length(mcols(gal1)$MC)>0){
      mcols(gal1)$MC[gp==2] <- cigar(gal1[id1])
      mcols(gal1)$MC[id1] <- cigar(gal1[id2])
    }
    gal1 <- gal1[!is.na(mcols(gal1)$mpos)] ## remove the single ends
    if(length(gal1)<1){
      message("Only single end reads")
      return(gal1)
    }
    ## till now, gal1 must have mrnm, mpos, names and flag
    stopifnot(length(mcols(gal1)$mrnm)>0)
    stopifnot(length(mcols(gal1)$mpos)>0)
    stopifnot(length(mcols(gal1)$flag)>0)
    stopifnot(length(names(gal1))>0)
    if(!(missing(outbam))){
      tryCatch(exportBamFile(gal1, outbam),
               error=function(e){
                 message(e)
               })
    }
    return(gal1)
}