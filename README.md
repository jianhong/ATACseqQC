# ATACseqQC

[![platforms](http://bioconductor.org/shields/availability/3.5/ATACseqQC.svg)](http://bioconductor.org/packages/devel/bioc/html/ATACseqQC.html)
[![build](http://bioconductor.org/shields/build/devel/bioc/ATACseqQC.svg)](http://bioconductor.org/packages/devel/bioc/html/ATACseqQC.html)
[![test coverage](https://codecov.io/github/Bioconductor-mirror/ATACseqQC/branch/master)](http://bioconductor.org/packages/devel/bioc/html/ATACseqQC.html)

ATAC sequencing Quality Control


ATAC-seq, an assay for Transposase-Accessible Chromatin using 
sequencing, is a rapid and sensitive method for chromatin accessibility 
analysis. It was developed as an alternative method to MNase-seq, FAIRE-seq 
and DNAse-seq. Comparing to the other methods, ATAC-seq requires less amount 
of the biological samples and time to process. In the process of analyzing 
several ATAC-seq dataset produced in our labs, we learned some of the unique 
aspects of the quality assessment for ATAC-seq data.To help users to quickly 
assess whether their ATAC-seq experiment is successful, we developed 
ATACseqQC package partially following the guideline published in Nature 
Method 2013 (Greenleaf et al.), including diagnostic plot of fragment size 
distribution, proportion of mitochondria reads, nucleosome positioning 
pattern, and CTCF or other Transcript Factor footprints.

## Installation

To install this package, start R and enter:

```r
library(BiocInstaller)
biocLite("ATACseqQC")
```

## Documentation

To view documentation of ATACseqQC, start R and enter:
```r
browseVignettes("ATACseqQC")
```

