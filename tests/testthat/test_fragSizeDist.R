test_that("fragSizeDist works not correct", {
  gal <- GAlignments(seqnames=Rle("chr1", 2020), pos=rep(1L, 2020),
                     cigar=rep("10M", 2020),
                     strand=Rle(strand(rep(c("+", "-"), 1010))), 
                     isize=rep(seq.int(1010), each=2) * c(1, -1),
                     qname=rep(seq.int(1010), each=2),
                     flag=rep(c(99, 147), 1010))
  seqlengths(gal) <- c("chr1"=2000)
  tmpfile <- "test.bam"
  export(gal, tmpfile)
  size <- fragSizeDist(tmpfile, "test")
  expect_true(all(size$test==2))
  expect_equal(length(size$test), 1010)
})