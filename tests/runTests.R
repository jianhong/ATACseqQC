require("ATACseqQC") || stop("unable to load Package:ATACseqQC")
require("GenomicAlignments") || stop("unable to load Package:GenomicAlignments")
require("testthat") || stop("unable to load testthat")
test_check("ATACseqQC")
