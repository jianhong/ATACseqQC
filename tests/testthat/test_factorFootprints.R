test_that("factorFootprints works not correct", {
    # data(annotatedPeak)
    # x<-addGeneIDs(annotatedPeak[1:6],"org.Hs.eg.db",c("symbol","omim"))
    # expect_is(x$symbol, "character")
    # expect_is(x$omim, "character")
    # expect_equal(x$feature[1:4], c("ENSG00000202254", "ENSG00000186086",
    #                           "ENSG00000065135", "ENSG00000197106"))
    # expect_equal(unname(x$symbol[1:4]), c(NA, "NBPF6", "GNAI3", "SLC6A17"))
    # x <- addGeneIDs(c("ENSG00000065135", "ENSG00000116396", 
    #                   "ENSG00000197106", "ENSG00000186086", 
    #                   "ENSG00000202254"), org.Hs.eg.db, c("symbol"))
    # expect_equal(unname(x$symbol[match(c("ENSG00000065135", "ENSG00000116396", 
    #                               "ENSG00000197106", "ENSG00000186086", 
    #                               "ENSG00000202254"), x$ensembl_gene_id)]), 
    #              c('GNAI3', 'KCNC4', 'SLC6A17', 'NBPF6', NA))
    # expect_error(addGeneIDs("ENSG00000065135", orgAnn="generate.error"))
})