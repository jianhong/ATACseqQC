test_that("shiftReads works not correct", {
  gal <- GAlignments(seqnames=Rle("chr1", 6), pos=1:6,
    cigar=c("47M", "50M", "50M", "50M", "50M", "50M"),
    strand=Rle(strand(c("+", "-", "+", "-")), c(1, 2, 2, 1)),
    qname=tail(letters, 6),
    qual=PhredQuality(c(
      "CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ",
      "JJJJJJJHHJJJJIJJIIIJIIIIIIIJJJIJJIJIHFHHHHFFFFFCCC",
      "FFFFFFGHHHHJJJJJJJIJIIJJJIHG@JJJJJJJHHHFFHFDEDFCB@",
      "@CCFFFFEHGGFHGIGGIJIJJJIJJJJJJJJIJJJJJJJJJIGJIJJIJ",
      "CC@FFFFFHHHHHJJJJIIJGIIJJIIGIJIJJJIIJGHHJIJJIJJIJE",
      "DHHIIIIJIHJJJJHIIIGJJGIJJJIJIIIIIIGIGHHHHHFFFFFCCC"
    )),
    isize=c(180, -180, -265, 265, 185, -185))
  mcols(gal)$seq=DNAStringSet(c(
    "ATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGAAATGATCTG",
    "GACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGT",
    "TGAATCTCTCAGGACTCAGCAGCAGCCCACCCCCTCAACGACCTCTCCAT",
    "GACCAAGTCCGCCTCTCAGCCTCACGTTCTTAAGCAGCATCTAAGCAAAA",
    "GCAGCAGCCCACCCCCTCAACGACCTCTCCATTGGCCAGAATCTGGGACA",
    "AAAGCCACTGCCCCTCAGACCGCAGTTGCGAGTAGCACGAGGATCTGGGG"
  ))
  shifted <- ATACseqQC:::shiftReads(gal, positive=4L, negative=5L)
  shifted2 <- ATACseqQC:::shiftReads(gal, positive=5L, negative=4L)
  ns <- as.character(strand(gal))=="+"
  expect_equal(start(shifted)[ns], start(gal)[ns]+4)
  expect_equal(end(shifted)[!ns], end(gal)[!ns]-5)
  expect_equal(width(shifted)[ns], width(gal)[ns]-4)
  expect_equal(width(shifted)[!ns], width(gal)[!ns]-5)
  expect_equal(abs(mcols(shifted)$isize), abs(mcols(gal)$isize)-9)
  expect_equal(width(mcols(shifted)$seq), qwidth(shifted))
  expect_equal(width(mcols(shifted)$qual), qwidth(shifted))
  expect_equal(start(shifted2)[ns], start(gal)[ns]+5)
  expect_equal(end(shifted2)[!ns], end(gal)[!ns]-4)
  expect_equal(width(shifted2)[ns], width(gal)[ns]-5)
  expect_equal(width(shifted2)[!ns], width(gal)[!ns]-4)
  expect_equal(abs(mcols(shifted2)$isize), abs(mcols(gal)$isize)-9)
})