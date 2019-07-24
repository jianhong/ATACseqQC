#' @title split bams into nucleosome free, mononucleosome,
#' dinucleosome and trinucleosome
#' @description use random forest to split the reads into nucleosome free,
#'              mononucleosome, dinucleosome and trinucleosome. 
#'              The features used in random forest including 
#'              fragment length, GC content, and 
#'              UCSC phastCons conservation scores.
#' @param obj an object of \link[GenomicAlignments:GAlignments-class]{GAlignments}
#' @param txs GRanges of transcripts
#' @param genome an object of BSgenome
#' @param conservation an object of \link[GenomicScores:GScores-class]{GScores}.
#' @param outPath folder to save the splitted alignments. If outPath is setting, 
#' the return of the function will not contain seq and qual fields.
#' @param breaks a numeric vector for fragment size of nucleosome free,
#' mononucleosome, dinucleosome and trinucleosome. The breaks pre-defined
#' here is following the description of Greenleaf's paper (see reference).
#' @param labels a character vector for labels of the levels 
#' of the resulting category.
#' @param labelsOfNucleosomeFree,labelsOfMononucleosome character(1). The label
#' for nucleosome free and mononucleosome.
#' @param trainningSetPercentage numeric(1) between 0 and 1. Percentage of 
#' trainning set from top coverage.
#' @param cutoff numeric(1) between 0 and 1. cutoff value for prediction.
#' @param halfSizeOfNucleosome numeric(1) or integer(1). Thre read length will
#' be adjusted to half of the nucleosome size to enhance the signal-to-noise
#' ratio.
#' @return a list of GAlignments
#' @author Jianhong Ou
#' @import BiocGenerics
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @importFrom Biostrings letterFrequency
#' @importFrom BSgenome getSeq
#' @importFrom ChIPpeakAnno annotatePeakInBatch
#' @importFrom GenomicAlignments coverage 
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @importFrom GenomicScores gscores
#' @export
#' @references Buenrostro, J.D., Giresi, P.G., Zaba, L.C., Chang, H.Y. and 
#' Greenleaf, W.J., 2013. Transposition of native chromatin for fast and 
#' sensitive epigenomic profiling of open chromatin, DNA-binding proteins and 
#' nucleosome position. Nature methods, 10(12), pp.1213-1218.
#' 
#' Chen, K., Xi, Y., Pan, X., Li, Z., Kaestner, K., Tyler, J., Dent, S., 
#' He, X. and Li, W., 2013. DANPOS: dynamic analysis of nucleosome position 
#' and occupancy by sequencing. Genome research, 23(2), pp.341-351.
#' @examples 
#' library(GenomicRanges)
#' bamfile <- system.file("extdata", "GL1.bam", 
#'                        package="ATACseqQC", mustWork=TRUE)
#' tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
#' gal1 <- readBamFile(bamFile=bamfile, tag=tags, 
#'                     which=GRanges("chr1", IRanges(1, 1e6)), 
#'                     asMates=FALSE)
#' names(gal1) <- mcols(gal1)$qname
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(phastCons100way.UCSC.hg19)
#' splitGAlignmentsByCut(gal1, txs=txs, genome=Hsapiens, 
#'                       conservation=phastCons100way.UCSC.hg19)
#' 
splitGAlignmentsByCut <- function(obj, txs, genome, conservation,
                                  outPath,
                          breaks=c(0, 100, 180, 247, 315, 473, 558, 615, Inf),
                          labels = c("NucleosomeFree", "inter1",
                                     "mononucleosome", "inter2",
                                     "dinucleosome", "inter3",
                                     "trinucleosome", "others"),
                          labelsOfNucleosomeFree="NucleosomeFree",
                          labelsOfMononucleosome="mononucleosome",
                          trainningSetPercentage=.15,
                          cutoff = .8, halfSizeOfNucleosome=80L){
  stopifnot(length(labels)+1==length(breaks))
  stopifnot(length(labelsOfMononucleosome)==1)
  stopifnot(length(labelsOfNucleosomeFree)==1)
  stopifnot(labelsOfMononucleosome %in% labels)
  stopifnot(labelsOfNucleosomeFree %in% labels)
  conservationFlag <- FALSE
  if(!missing(conservation)){
    if(length(conservation)){
      stopifnot(is(conservation, "GScores"))
      conservationFlag <- TRUE ## conservation is supplied.
    }
  }
  stopifnot(is(obj, "GAlignments"))
  
  if(length(obj)==0 && (!missing(outPath)) && (!conservationFlag)){
    ## big file mode
    meta <- metadata(obj)
    if(!all(c("file", "param") %in% names(meta))){
      stop("length of obj could not be 0.")
    }
    
    ow <- getOption("warn")
    on.exit(options(warn = ow))
    options(warn=-1)
    chunk <- 100000
    index <- ifelse(length(meta$index)>0, meta$index, meta$file)
    bamfile <- BamFile(meta$file, index=index, yieldSize=chunk, asMates = meta$asMates)
    outfile <- list()
    mergedfile <- list()
    mpos <- NULL
    open(bamfile)
    if(length(meta$mpos)>0){
      mpos <- meta$mpos
    }
    while (length(chunk0 <- readGAlignments(bamfile, param=meta$param))) {
      if(length(mpos)){
        mcols(chunk0)$MD <- NULL
        names(chunk0) <- mcols(chunk0)$qname
        chunk0 <- chunk0[order(names(chunk0))]
        mcols(chunk0)$mpos <- mpos[paste(mcols(chunk0)$qname, start(chunk0))]
      }
      gal1 <- splitGAlignmentsByCut(chunk0, txs=txs, genome = genome,
                                    outPath = outPath, breaks = breaks,
                                    labels = labels, labelsOfNucleosomeFree = labelsOfNucleosomeFree,
                                    labelsOfMononucleosome = labelsOfMononucleosome,
                                    trainningSetPercentage = trainningSetPercentage,
                                    cutoff = cutoff, halfSizeOfNucleosome = halfSizeOfNucleosome)
      ## gal1 is a list of GAlignments
      ## save the gal1 to tmpfiles
      for(i in seq_along(gal1)){
        outfile[[names(gal1)[i]]] <- c(tempfile(fileext = ".bam"), outfile[[names(gal1)[i]]])
        export(gal1[[i]], outfile[[names(gal1)[i]]][1], format="BAM")
      }
      rm(gal1)
    }
    close(bamfile)
    for(i in seq_along(outfile)){
      if(length(outfile[[i]])>1){
        mergedfile[[names(outfile)[i]]] <- mergeBam(outfile[[i]], 
                                                    destination=file.path(outPath, 
                                                                          paste0(names(outfile)[i], ".bam")), 
                                                    indexDestination=TRUE)
        unlink(outfile[[i]])
        unlink(paste0(outfile[[i]], ".bai"))
      }else{
        file.copy(outfile[[i]], 
                  file.path(outPath, paste0(names(outfile)[i], ".bam")))
        file.copy(paste0(outfile[[i]], ".bai"), 
                  file.path(outPath, paste0(names(outfile)[i], ".bam.bai")))
        mergedfile[[names(outfile)[i]]] <- 
          file.path(outPath, paste0(names(outfile)[i], ".bam"))
        unlink(outfile[[i]])
        unlink(paste0(outfile[[i]], ".bai"))
      }
    }
    meta$param <- ScanBamParam(flag=meta$param@flag, what=c("qname", "flag", "mapq", "isize"))
    objs <- lapply(mergedfile, function(.ele){
      ## read the bam files
      chunk0 <- readGAlignments(.ele, param=meta$param)
      names(chunk0) <- mcols(chunk0)$qname
      chunk0 <- chunk0[order(names(chunk0))]
      mcols(chunk0)$mpos <- mpos[paste(mcols(chunk0)$qname, start(chunk0))]
      chunk0
    })
    return(objs)
  }
  
  if(length(obj)==0){
    obj <- loadBamFile(obj)
  }
  stopifnot(length(names(obj))==length(obj))
  
  seqlev <- unique(seqlevels(obj))
  seqlev <- seqlev[seqlev %in% unique(seqnames(obj))]
  if(length(seqlev)==0){
    stop("no reads in the obj.")
  }
  objs <- split(obj,
                cut(abs(mcols(obj)$isize),
                    breaks = breaks,
                    labels = labels))
  if(!conservationFlag){
    if(!missing(outPath)){
      writeListOfGAlignments(objs, outPath = outPath)
    }
    return(objs)
  }
  stopifnot(is(txs, "GRanges"))
  stopifnot(is(genome, "BSgenome"))
  nf.old <- objs[[labelsOfNucleosomeFree]]
  nc.old <- objs[[labelsOfMononucleosome]]
  rm(objs)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  nf.cvg <- coverage(nf.old)
  nc.cvg <- coverage(nc.old)
  nf.cvg <- nf.cvg[seqlev]
  nc.cvg <- nc.cvg[seqlev]
  nf.cvg.quantile <- sapply(nf.cvg, function(.ele) {
    .ele <- table(.ele)
    .ele <- .ele[-1]
    .x <- cumsum(.ele)/sum(.ele)
    as.numeric(names(.x[which(.x>=1-trainningSetPercentage)[1]]))
  })
  nc.cvg.quantile <- sapply(nc.cvg, function(.ele) {
      .ele <- table(.ele)
      .ele <- .ele[-1]
      .x <- cumsum(.ele)/sum(.ele)
      as.numeric(names(.x[which(.x>=1-trainningSetPercentage)[1]]))
  })
  cvg.quantile <- ifelse(nf.cvg.quantile < nc.cvg.quantile,
                         nf.cvg.quantile, nc.cvg.quantile)
  #names(cvg.quantile) <- seqlev
  nf.cvg.view <-
    as(IRangesList(mapply(function(x, y) as(Views(x, x>y), "IRanges"),
                          nf.cvg, cvg.quantile)), "GRanges")
  nc.cvg.view <-
    as(IRangesList(mapply(function(x, y) as(Views(x, x>y), "IRanges"),
                          nc.cvg, cvg.quantile)), "GRanges")
  #annotation the GRanges to get the strand info
  TSS <- promoters(txs,
                   upstream = 0,
                   downstream = 1)
  rm(txs)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  TSS <- TSS[seqnames(TSS) %in% seqlev]
  seqlev <- seqlev[seqlev %in% unique(seqnames(TSS))]
  seqlevels(TSS) <- seqlev
  seqinfo(TSS) <- Seqinfo(seqlev, seqlengths = seqlengths(TSS))
  TSS <- unique(TSS)
  nf.anno <- annotatePeakInBatch(nf.cvg.view, AnnotationData=TSS,
                                 output="nearestLocation")
  nc.anno <- annotatePeakInBatch(nc.cvg.view, AnnotationData=TSS,
                                 output="nearestLocation")
  nf.anno$feature_strand[is.na(nf.anno$feature_strand)] <- "*"
  nc.anno$feature_strand[is.na(nc.anno$feature_strand)] <- "*"
  strand(nf.anno) <- nf.anno$feature_strand
  strand(nc.anno) <- nc.anno$feature_strand
  mcols(nf.anno) <- DataFrame(score=1)
  mcols(nc.anno) <- DataFrame(score=1)
  nx <- GRoperator(nf.anno, nc.anno, operator = "-")
  nx <- nx[width(nx)>=40]
  nf <- nx[nx$score>0]
  nc <- nx[nx$score<0]
  if(length(nc)<1){
    stop("not enough mononucleosome reads for training! Just try without conservation score.")
  }
  if(length(nf)<1){
    stop("not enough nucleosome free reads for training! Just try without conservation score.")
  }
  block <- 100000
  nf <- reCenterPeaks(nf, width=halfSizeOfNucleosome)
  nc <- reCenterPeaks(nc, width=halfSizeOfNucleosome)
  coverage.nf <-
    Views(nf.cvg[seqlev], split(ranges(nf),
                                as.character(seqnames(nf)))[seqlev])
  coverage.nf <- IRangesList(lapply(coverage.nf, function(.ele){
    ir <- as(.ele, "IRanges")
    mcols(ir)$coverage <- viewMeans(.ele)
    ir
  }))
  nf <- as(coverage.nf, "GRanges")
  nf$score <- unlist(lapply(coverage.nf, function(.ele)
    as.numeric(mcols(.ele)$coverage)))
  nf <- nf[order(nf$score, decreasing=TRUE)]
  nf <- nf[seq_len(min(block, length(nf)))]
  coverage.nc <-
    Views(nc.cvg[seqlev], split(ranges(nc),
                                as.character(seqnames(nc)))[seqlev])
  coverage.nc <- IRangesList(lapply(coverage.nc, function(.ele){
    ir <- as(.ele, "IRanges")
    mcols(ir)$coverage <- viewMeans(.ele)
    ir
  }))
  nc <- as(coverage.nc, "GRanges")
  nc$score <- unlist(lapply(coverage.nc, function(.ele)
    as.numeric(mcols(.ele)$coverage)))
  nc <- nc[order(nc$score, decreasing=TRUE)]
  nc <- nc[seq_len(min(block, length(nc)))]
  obj <- obj[seqnames(obj) %in% seqlev]
  newdata <- split(as(obj, "GRanges"), names(obj))
  newdata <- range(newdata, ignore.strand=TRUE)
  newdata <- unlist(newdata)
  seqlevels(newdata) <- seqlev
  seqinfo(newdata) <- seqinfo(TSS)
  nd <- GenomicRanges::trim(newdata)
  rm(newdata)
  nf.seq <- getSeq(genome, nf)
  nc.seq <- getSeq(genome, nc)
  # gc ratios
  nf.gc <- letterFrequency(nf.seq, letters="CG", as.prob = TRUE)
  nc.gc <- letterFrequency(nc.seq, letters="CG", as.prob = TRUE)
  
  nd.gc <- lapply(split(nd, rep(as.character(seq.int(ceiling(length(nd)/block))),
                                each=block)[seq_along(nd)]), 
                   function(.ele){
                     letterFrequency(getSeq(genome, .ele), letters="CG", as.prob = TRUE)
                   })
  nd.gc <- nd.gc[order(as.numeric(names(nd.gc)))]
  nd.gc <- do.call(rbind, nd.gc)
  rm(nf.seq)
  rm(nc.seq)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  # conservation
  getScoresFromCons <- function(cons, gr){# this step too slow
      gr.cons <- gscores(x=cons, ranges=gr)
      stopifnot(identical(ranges(gr), ranges(gr.cons)))
      gr.cons <- gr.cons$scores
      gr.cons[is.na(gr.cons)] <- 0
      gr.cons
  }
  nf.conservation <- getScoresFromCons(conservation, nf)
  nc.conservation <- getScoresFromCons(conservation, nc)
  nd.conservation <- getScoresFromCons(conservation, nd)
  # median fragment length
  nf.old <- as(nf.old, "GRanges")
  nc.old <- as(nc.old, "GRanges")
  nf.ol <- findOverlaps(nf, nf.old)
  nc.ol <- findOverlaps(nc, nc.old)
  nf.ol <- split(width(nf.old)[subjectHits(nf.ol)], queryHits(nf.ol))
  nc.ol <- split(width(nc.old)[subjectHits(nc.ol)], queryHits(nc.ol))
  rm(nf.old)
  rm(nc.old)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  stopifnot(length(nf.ol)==length(nf))
  stopifnot(length(nc.ol)==length(nc))
  nf.frag.len <- sapply(nf.ol, median)
  nc.frag.len <- sapply(nc.ol, median)
  nf.frag.len <- nf.frag.len[order(as.numeric(names(nf.frag.len)))]
  nc.frag.len <- nc.frag.len[order(as.numeric(names(nc.frag.len)))]
  nd.frag.len <- width(nd)
  x.nf <- cbind(nf.frag.len, nf.conservation, nf.gc)
  x.nc <- cbind(nc.frag.len, nc.conservation, nc.gc)
  x <- rbind(x.nf, x.nc)
  testdata <- cbind(nd.frag.len, nd.conservation, nd.gc)
  colnames(testdata) <- colnames(x) <- NULL
  y <- as.factor(rep(c("f", "n"), c(length(nf), length(nc))))
  fit <- randomForest(x, y,
                      ntree=2*ceiling(sqrt(length(y))),
                      keep.inbag=TRUE)
  pred <- predict(fit, newdata=testdata, type='prob')
  pred.free <- pred[, "f"]
  pred.bind <- pred[, "n"]
  rm(nd.gc, nf.gc, nc.gc)
  rm(pred)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  nucleosome <- names(obj) %in% names(nd[pred.bind >= cutoff])
  nucleosomefree <- names(obj) %in%
    names(nd[pred.free>= cutoff & pred.bind < cutoff])
  NucleosomeFree <- obj[nucleosomefree]
  Nucleosome <- obj[nucleosome]
  left <- obj[!(nucleosome | nucleosomefree)]
  rm(nd)
  rm(obj)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  objs <- split(left,
                cut(abs(mcols(left)$isize),
                    breaks = breaks,
                    labels = labels))
  rm(left)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  Nucleosome <- split(Nucleosome,
                      cut(abs(mcols(Nucleosome)$isize),
                          breaks = breaks,
                          labels = labels))
  Nucleosome[[labelsOfNucleosomeFree]] <- NucleosomeFree
  rm(NucleosomeFree)
  gc(verbose = FALSE, reset = TRUE, full = TRUE)
  for(i in intersect(names(Nucleosome), names(objs))){
    objs[[i]] <- c(Nucleosome[[i]], objs[[i]])
  }
  if(!missing(outPath)){
    writeListOfGAlignments(objs, outPath = outPath)
  }
  objs
}


GRbin <- function(A, B, col){
    A1 <- B1 <- C <- disjoin(c(A, B), ignore.strand=FALSE)
    ol1 <- findOverlaps(C, A, minoverlap=1L,
                        type="any", ignore.strand=FALSE)
    ol2 <- findOverlaps(C, B, minoverlap=1L,
                        type="any", ignore.strand=FALSE)
    A1$value <- B1$value <- 0
    A1$value[queryHits(ol1)] <- mcols(A)[subjectHits(ol1), col]
    B1$value[queryHits(ol2)] <- mcols(B)[subjectHits(ol2), col]
    GRangesList(A=A1, B=B1)
}
GRoperator <- function(A, B, col="score", 
                       operator=c("+", "-", "*", "/", "^", "%%")){
    if(!is(A, "GRanges") || !is(B, "GRanges")){
        stop("A and B must be objects of GRanges")
    }
    operator <- match.arg(operator)
    if(!(col %in% colnames(mcols(A)))){
        stop("col is not in metadata of A")
    }
    if(!(col %in% colnames(mcols(B)))){
        stop("col is not in metadata of B")
    }
    C <- GRbin(A, B, col)
    A <- C$A
    B <- C$B
    out <- A
    operator <- .Primitive(operator)
    out$value <- operator(A$value, B$value)
    colnames(mcols(out)) <- col
    out
}
