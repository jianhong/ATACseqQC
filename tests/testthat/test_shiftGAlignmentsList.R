test_that("shiftGAlignmentsList works not correct", {
  bamfile <- system.file("extdata", "GL1.bam", package="ATACseqQC")
  tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
  which <- GRanges('chr1:1-249250621:*')
  gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)
  objs1 <- shiftGAlignmentsList(gal)
  gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE,
                     bigFile = TRUE)
  objs2 <- shiftGAlignmentsList(gal)
  x <- paste(mcols(objs1)$qname, start(objs1), end(objs1))
  y <- paste(mcols(objs2)$qname, start(objs2), end(objs2))
  objs2 <- objs2[match(x, y)]
  expect_equal(cigar(objs1), cigar(objs2))
  expect_equal(mcols(objs1)$seq, mcols(objs1)$seq)
  expect_equal(mcols(objs1)$qual, mcols(objs1)$qual)
  expect_equal(mcols(objs1)$mpos, mcols(objs1)$mpos)
})
