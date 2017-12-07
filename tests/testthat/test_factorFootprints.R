test_that("factorFootprints works not correct", {
  if(interactive()){
    set.seed(0)
    fa <- sample(c("A", "C", "G", "T"), 500, replace = TRUE)
    rg <- 245:255
    fa[rg] <- "A"
    fa <- c(">chr1", paste(fa, collapse = ""))
    dir.create("unknownGenome")
    writeLines(fa, con = file.path("unknownGenome", "chr1.fa"))
    seedsFile <- c("Package: BSgenome.unknown",
                   "Title: Full genome for test use",
                   "Description: Full genome sequences for test use.",
                   "Version: 0.0.1",
                   "organism: unknown",
                   "common_name: unknown",
                   "provider: unknown",
                   "provider_version: unknown",
                   "release_date: Apr. 2017",
                   "release_name: unknown",
                   "source_url: http://test.org",
                   "organism_biocview: unknown",
                   "BSgenomeObjname: unknown",
                   "seqnames: \"chr1\"",
                   "SrcDataFiles: chr1.fa from code",
                   "PkgExamples: genome$chr1",
                   "seqs_srcdir: unknownGenome")
    writeLines(seedsFile, con = "genome.seed")
    
    forgeBSgenomeDataPkg("genome.seed")
    
    install.packages("BSgenome.unknown", repos = NULL, type = "source")
    on.exit({
      remove.packages("BSgenome.unknown")
      unlink(c("test.bam", "test.bam.bai", "BSgenome.unknown", "unknownGenome", "genome.seed"),
             recursive = TRUE)
    })
    library("BSgenome.unknown")
    genome <- BSgenome.unknown
    
    gal <- GAlignments(seqnames=Rle("chr1", 50), 
                       pos=as.integer(c(seq(150, 240, by=10), seq(200, 245, by=5), seq(205, 245, by=10),
                                        seq(255, 295, by=10), seq(255, 300, by=5), seq(260, 350, by=10))),
                       cigar=rep("10M", 50),
                       strand=Rle(strand("+"), 50))
    seqlengths(gal) <- c("chr1"=1000)
    tmpfile <- "test.bam"
    export(gal, tmpfile)
    pfm <- matrix(rep(c(1, 0, 0, 0), 10), nrow=4)
    rownames(pfm) <- c("A", "C", "G", "T")
    ffp <- factorFootprints(tmpfile, index=tmpfile, pfm, genome, 
                            min.score="95%", seqlev="chr1", 
                            upstream=100, downstream=100)
    expect_equal(ffp$estLibSize, 50)
    expect_true(all(colMeans(ffp$signal[["-"]])==0))
    expect_equal(sum(colMeans(ffp$signal[["+"]])), 50)
  }
})
