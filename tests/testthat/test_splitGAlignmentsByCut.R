test_that("splitGAlignmentsByCut works not correct", {
  gal <- GAlignments(seqnames=Rle("chr1", 6), pos=1:6,
                     cigar=c("50M", "50M", "50M", "50M", "50M", "50M"),
                     strand=Rle(strand(c("+", "-", "+", "-")), c(1, 2, 2, 1)),
                     qname=tail(letters, 6),
                     isize=c(180, -180, -265, 265, 185, -185))
  expect_error(splitGAlignmentsByCut(gal))
  names(gal) <- mcols(gal)$qname
  objs <- splitGAlignmentsByCut(gal)
  expect_equal(unname(lengths(objs)), c(0, 2, 2, 2, 0, 0, 0, 0))
})
