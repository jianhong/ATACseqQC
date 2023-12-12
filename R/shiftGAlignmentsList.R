#' @title shift 5' ends
#' @description shift the GAlignmentsLists by 5' ends. 
#' All reads aligning to the positive strand will be offset by +4bp, 
#' and all reads aligning to the negative strand will be offset -5bp by default.
#' @param gal An object of \link[GenomicAlignments:GAlignmentsList-class]{GAlignmentsList}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @param outbam file path to save shift reads. If missing, no file will be write.
#' @param slidingWindowSize The width of each tile when the input is big file.
#' By default 50e6, the memory cost will be about 6GB for each thread for a 5GB bam file.
#' Increase the value will increase the memory cost but may speed up the process.
#' @param BPPARAM The parallel parameters used by BiocParrallel.
#' @return An object of \link[GenomicAlignments:GAlignments-class]{GAlignments} with 5' end 
#' shifted reads. The PCR duplicated will be removed unless there is metadata
#' keepDuplicates set to TRUE.
#' @author Jianhong Ou
#' @export
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom Rsamtools mergeBam bamWhich `bamWhich<-` filterBam bamTag
#' `bamTag<-`
#' @importFrom rtracklayer export
#' @importFrom utils read.csv write.csv combn txtProgressBar setTxtProgressBar
#' @importFrom BiocParallel bptry bplapply
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' which <- as(seqinfo(Hsapiens)["chr1"], "GRanges")
#' gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)
#' objs <- shiftGAlignmentsList(gal)
#' export(objs, "shift.bam")
#' \dontrun{
#'   bamfile <- 'a.big.file.bam'
#'   tags <- c("NM", "MD")
#'   which <- GRanges(c('chr1:1-249250621:*', 'chr2:1-243199373:*'))
#'   gal <- readBamFile(bamfile, tag=tags, which=which, 
#'   asMates=TRUE, bigFile=TRUE)
#'   library(BiocParallel)
#'   BPPARAM <- MulticoreParam(workers = 2, progress=TRUE)
#'   objs <- shiftGAlignmentsList(gal, BPPARAM = BPPARAM,
#'    outbam="shift.bam")
#' }
shiftGAlignmentsList <- function(gal, positive=4L, negative=5L, outbam,
                                 BPPARAM = NULL,
                                 slidingWindowSize = 50e6){
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
      index <- ifelse(length(meta$index)>0, meta$index, meta$file)
      bamW <- bamWhich(meta$param)
      bamW <- as(bamW, "GRanges")
      bamW <- slidingWindows(bamW,
                             width=slidingWindowSize,
                             step=slidingWindowSize)
      bamW <- unlist(bamW)
      rmBam <- function(bam, bai=paste0(bam, ".bai")){
        unlink(bam)
        if(file.exists(bai[1])) unlink(bai)
      }
      mposTag <- vapply(combn(LETTERS, 2, simplify = FALSE), function(.ele){
        paste(.ele, collapse = '')
      }, FUN.VALUE = character(1L))
      if(any(!mposTag %in% bamTag(meta$param))){
        mposTag <- mposTag[!mposTag %in% bamTag(meta$param)][1]
      }else{
        stop('Can not find a proper tag to save mpos! Please remove some tags.')
      }
      processBamByTile <- function(meta, i, bamW, rmBam, mposTag){
        sbp <- meta$param
        bamWhich(sbp) <- bamW[i]
        destination <- tempfile(fileext = ".bam")
        null <- filterBam(file=meta$file, index = index,
                          destination=destination,
                          indexDestination=TRUE,
                          param=sbp)
        chunk <- 1000000
        bamfile <- BamFile(destination,
                           index=destination,
                           yieldSize=chunk,
                           asMates = meta$asMates)
        open(bamfile)
        on.exit({
          close(bamfile)
          rmBam(destination)
          })
        outfile <- NULL
        SE <- GAlignments()
        df4Duplicates <- NULL
        while (length(chunk0 <-
                      readGAlignmentsList(bamfile, param=meta$param))) {
          isSE <- lengths(chunk0)==1
          if(sum(isSE)){
            SE <- c(SE, unlist(chunk0[isSE]))
            chunk0 <- chunk0[!isSE]
            ## pair SE
            SE <- split(SE, mcols(SE)$qname)
            isSE <- lengths(SE)==1
            chunk0 <- c(chunk0, SE[!isSE])
            SE <- SE[isSE]
            SE <- unlist(SE)
            if(length(chunk0)==0){
              next
            }
          }
          startpos <- vapply(start(chunk0), min,
                             FUN.VALUE = integer(1L))
          chunk0 <- chunk0[order(startpos)]
          metadata(chunk0)$df4Duplicates <- df4Duplicates
          gal1 <- shiftGAlignmentsList(chunk0, positive = positive, 
                                       negative = negative)
          rm(chunk0)
          gc(verbose=FALSE)
          if(length(metadata(gal1)$df4Duplicates)){
            df4Duplicates <- read.csv(metadata(gal1)$df4Duplicates)
          }
          if(length(mcols(gal1)$mpos)!=length(gal1)){
            stop("Can not get mpos info from the reads.")
          }
          mcols(gal1)[, mposTag] <- mcols(gal1)$mpos
          outfile <- c(tempfile(fileext = ".bam"), outfile)
          if(length(meta$header)>0) metadata(gal1)$header <- meta$header
          exportBamFile(gal1, outfile[1])
          rm(gal1)
          gc(verbose=FALSE)
        }
        close(bamfile)
        rmBam(destination)
        on.exit()
        if(length(SE)>0){
          singleEnds <- 
            tempfile(fileext = paste0('.',
                                      as.character(seqnames(SE[1])),
                                      '.bam'))
          SE <- renameMcol(SE, "mrnm", tagNewName=mposTag)
          exportBamFile(SE, singleEnds)
          rm(SE)
          gc(verbose=FALSE)
        }else{
          singleEnds <- NULL
        }
        return(list(outfile=outfile,
                    singleEnds = singleEnds))
      }
      if(is.null(BPPARAM)){
        pb <- txtProgressBar(max = length(bamW), style = 3)
        res <- lapply(
          seq_along(bamW),
          function(i){
            .res <- processBamByTile(meta=meta,
                                     i=i,
                                     bamW=bamW,
                                     rmBam=rmBam,
                                     mposTag=mposTag)
            setTxtProgressBar(pb, i)
            .res
            })
        close(pb)
      }else{
        retryBatch <- 3
        retry <- TRUE
        res <- list()
        while(retryBatch>0){
          retryBatch <- retryBatch - 1
          if(any(retry)){
            res <- bptry(bplapply(seq_along(bamW),
                                  processBamByTile,
                                  meta=meta, bamW=bamW,
                                  rmBam=rmBam, mposTag=mposTag,
                                  BPREDO = res,
                                  BPPARAM = BPPARAM))
            retry <- vapply(res, FUN = function(.ele){
              is(.ele, 'bperror')
            }, FUN.VALUE = logical(1L))
          }else{
            retryBatch <- -1
          }
        }
      }
      options(warn = ow)
      on.exit()
      getRes <- function(res, key){
        unlist(lapply(res, function(.ele) .ele[[key]]))
      }
      outfile <- getRes(res, 'outfile')
      SE <- getRes(res, 'singleEnds')
      SE <- SE[lengths(SE)>0]
      SE_chr <- vapply(
        strsplit(SE, split = '.', fixed = TRUE),
        FUN = function(.ele){
          .ele[length(.ele)-1]
        },
        FUN.VALUE = character(1L)
      )
      SE <- split(SE, SE_chr)
      param <- meta$param
      bamTag(param) <- c(bamTag(param), mposTag)
      SE_outfile <- lapply(SE, function(se){
        chr_SE <- lapply(se, readGAlignments, param=param)
        chr_SE <- unlist(GAlignmentsList(chr_SE))
        chr_SE <- renameMcol(chr_SE, mposTag, tagNewName="mrnm")
        rmBam(se)
        this_outfile <- NULL
        if(length(chr_SE)){
          chr_SE <- split(chr_SE, mcols(chr_SE)$qname)
          isSE <- lengths(chr_SE)==1
          chunk0 <- chr_SE[!isSE]
          rm(chr_SE)
          if(length(chunk0)){
            startpos <- vapply(start(chunk0), min,
                               FUN.VALUE = integer(1L))
            chunk0 <- chunk0[order(startpos)]
            gal1 <- shiftGAlignmentsList(chunk0, positive = positive,
                                         negative = negative)
            rm(chunk0)
            gc(verbose=FALSE)
            if(length(gal1)){
              if(length(mcols(gal1)$mpos)!=length(gal1)){
                stop("Can not get mpos info from the reads.")
              }
              mcols(gal1)[, mposTag] <- mcols(gal1)$mpos
              this_outfile <- tempfile(fileext = ".bam")
              if(length(meta$header)>0) metadata(gal1)$header <- meta$header
              exportBamFile(gal1, this_outfile)
            }
            rm(gal1)
            gc(verbose=FALSE)
          }
        }
        return(this_outfile)
      })
      outfile <- c(outfile, unlist(SE_outfile))
      
      if(length(outfile)>1){
        BAI <- paste0(outfile[1], ".bai")
        mergedfile <- mergeBam(outfile, 
                               destination=tempfile(fileext = ".bam"), 
                               indexDestination=file.exists(BAI),
                               header=meta$file)
        rmBam(outfile)
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
        ## add mpos tag to the bamTag to let the readBamFile known
        meta$mposTag <- mposTag
        bamTag(meta$param) <- c(bamTag(meta$param), mposTag)
        metadata(gal1) <- meta
      }else{
        param <- meta$param
        bamTag(param) <- c(bamTag(param), mposTag)
        gal1 <- readGAlignments(mergedfile, param = param)
        mcols(gal1)$MD <- NULL
        names(gal1) <- mcols(gal1)$qname
        gal1 <- gal1[order(names(gal1))]
        gal1 <- renameMcol(gal1, mposTag, tagNewName="mpos")
        rmBam(mergedfile)
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