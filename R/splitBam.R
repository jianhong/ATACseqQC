#' @title prepare bam files for downstream analysis
#' @description shift the bam files by 5'ends and split the bam files.
#' @param bamfile character(1). File name of bam.
#' @param tags A vector of characters indicates the tags in bam file.
#' @param outPath Output file path.
#' @param txs \link[GenomicRanges]{GRanges} of transcripts.
#' @param genome An object of \link[BSgenome]{BSgenome}
#' @param conservation An object of \link[GenomicScores]{GScores}.
#' @param positive integer(1). the size to be shift for positive strand
#' @param negative integer(1). the size to be shift for negative strand
#' @param breaks A numeric vector for fragment size of nucleosome freee,
#' mononucleosome, dinucleosome and trinucleosome
#' @param labels A vector of characters indicates the labels for the levels
#' of the resulting category.
#' The length of labels = length of breaks - 1
#' @param seqlev A vector of characters indicates the sequence levels.
#' @param cutoff numeric(1). Cutoff value for prediction by
#' \link[randomForest]{randomForest}.
#' @return an invisible list of \link[GenomicAlignments]{GAlignments}
#' @author Jianhong Ou
#' @export
#' @import GenomeInfoDb
#' @importFrom Rsamtools scanBamFlag
#' @seealso \link{shiftGAlignmentsList}, \link{splitGAlignmentsByCut}, and 
#' \link{writeListOfGAlignments}
#' @examples
#' bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
#' tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(phastCons100way.UCSC.hg19)
#' objs <- splitBam(bamfile, tags,
#'                  txs=txs, genome=Hsapiens,
#'                  conservation=phastCons100way.UCSC.hg19,
#'                  seqlev="chr1")
#'

splitBam <- function(bamfile, tags, outPath=NULL,
                     txs, genome, conservation,
                     positive=4L, negative=5L,
                     breaks=c(0, 100, 180, 247, 315, 473, 558, 615, Inf),
                     labels = c("NucleosomeFree", "inter1",
                                "mononucleosome", "inter2",
                                "dinucleosome", "inter3",
                                "trinucleosome", "others"),
                     seqlev=paste0("chr", c(1:22, "X", "Y")),
                     cutoff = .8){
  stopifnot(length(labels)+1==length(breaks))
  stopifnot(is(txs, "GRanges"))
  stopifnot(is(conservation, "GScores"))
  stopifnot(is(genome, "BSgenome"))
  stopifnot(length(seqlev)>0)
  stopifnot(is.integer(positive))
  stopifnot(is.integer(negative))
  ## prepare for output
  if(!is.null(outPath)){
    stopifnot(length(outPath)==1)
    if(!file.exists(outPath)){
      dir.create(outPath, showWarnings = FALSE, recursive = TRUE)
    }
    stopifnot(file.exists(outPath))
    shiftPath <- file.path(outPath, "shifted.bam")
    if(file.exists(shiftPath)){
      stop("File ", shiftPath, " exits!")
    }
    for(i in labels){
      if(file.exists(file.path(outPath, paste0(i, ".bam")))){
        stop("File ", file.path(outPath, paste0(i, ".bam")), " exits!")
      }
    }
  }

  ## shift
  which <- as(seqinfo(genome)[seqlev], "GRanges")
  gal <- readBamFile(bamfile, tag=tags, which=which, 
                     what=c("qname", "flag", "mapq", "isize", 
                            "seq", "qual", "mrnm"),
                     flag=scanBamFlag(),
                     asMates=TRUE)
  gal1 <- shiftGAlignmentsList(gal,
                   positive=positive, negative=negative)
  ## split
  objs <- splitGAlignmentsByCut(gal1, breaks=breaks, labels=labels,
                        txs=txs, genome=genome, conservation=conservation,
                        cutoff=cutoff)
  ## output
  if(!is.null(outPath)){
      null <- writeListOfGAlignments(objs, outPath)
  }
  objs$all <- gal1
  return(invisible(objs))
}
