test_that("bamQC works not correct", {
  gal <- GAlignments(seqnames=Rle(c("chr1", "chrM"), c(990, 10)), 
                     pos=c(rep(c(1L, 241L), 491L) + rep(seq.int(491), each=2),
                           rep(c(1L, 241L), 9L)),
                     cigar=rep("10M", 1000),
                     strand=Rle(strand(rep(c("+", "-"), 500))), 
                     isize=rep(250, 1000) * c(1, -1),
                     flag=rep(c(99, 147), 500))
  names(gal) <- rep(seq.int(500), each=2)
  seqlengths(gal) <- c("chr1"=2000, "chrM"=1000)
  tmpfile <- "test.bam"
  export(gal, tmpfile)
  x <- bamQC(tmpfile)
  exp <- list(totalQNAMEs=500,
              mitochondriaRate=0.01,
              nonRedundantFraction=491/500,
              PCRbottleneckCoefficient_1=491/493,
              PCRbottleneckCoefficient_2=491)
  for(i in names(exp)){
    expect_equal(x[[i]], exp[[i]])
  }
})