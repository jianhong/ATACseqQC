# CHANGES IN VERSION 1.2.6

* improve bamQC

# CHANGES IN VERSION 1.2.5

* output NRF, PBC1, and PBC2 from bamQC.
* add unit test for bamQC.
* fix a bug in factorFootprints when motif is not exactly pfm.

# CHANGES IN VERSION 1.2.4

* add properPairRate, unmappedRate, hasUnmappedMateRate, 
  notPassingQualityControlsRate in output of bamQC.

# CHANGES IN VERSION 1.2.3

* Add mapq summary in output of bamQC.
* Add unit test for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.2.2

* Fix a bug in factorFootprints when bindingSites is supplied.
* adjust the Cut-site probability by the depth=librarySize/libraryCoverageSize for output of factorFootprints
* Modified the vignettes.
* Expand the functionality of the bamQC function.
* Import motifStack.
* add unit test.

# CHANGES IN VERSION 1.2.1

* fix a bug in factorFootprints

# CHANGES IN VERSION 1.1.17

* add gal argument for enrichedFragments to improve the efficency.

# CHANGES IN VERSION 1.1.16

* Fix a bug in fragSizeDist

# CHANGES IN VERSION 1.1.15

* Fix a bug in bamQC

# CHANGES IN VERSION 1.1.14

* Fix the error: when 'type' is "any", at least one of 'maxgap' and 
'minoverlap' must be set to its default value

# CHANGES IN VERSION 1.1.13

* improve the efficency of factorFootprints.

# CHANGES IN VERSION 1.1.12

* fix the soft clipping and improve the efficency.

# CHANGES IN VERSION 1.1.11

* fix a bug in fragSizeDist.

# CHANGES IN VERSION 1.1.10

* Check is PE before run fragSizeDist.

# CHANGES IN VERSION 1.1.9

* Remove duplicates for outputs of bamQC.

# CHANGES IN VERSION 1.1.8

* Add function bamQC.

# CHANGES IN VERSION 1.1.7

* Fix a bug in shiftGAlignmentsList.

# CHANGES IN VERSION 1.1.6

* Change default behavior of splitBam

# CHANGES IN VERSION 1.1.5

* update the vignette.

# CHANGES IN VERSION 1.1.4

* change author list.

# CHANGES IN VERSION 1.1.3

* change output of factorFootprints.

# CHANGES IN VERSION 1.1.2

* add cumulativePercentage in vignette

# CHANGES IN VERSION 1.1.1

* add new author.

# CHANGES IN VERSION 0.99.8

* reduce the memory cost.

# CHANGES IN VERSION 0.99.7

* update phastCons to GScores.

# CHANGES IN VERSION 0.99.6

* update documentation
* change package Name from ATACqc to ATACseqQC
* remove unused dependence
* change class checking from `class` to `is`

# CHANGES IN VERSION 0.99.5

* add bindingSites argument for factorFootprints

# CHANGES IN VERSION 0.99.4

* change shiftBam to shiftGAlignmentsList
* drop seqlevs arguments from functions

# CHANGES IN VERSION 0.99.3

* try to avoid error in windows that splitBam ask too much memory

# CHANGES IN VERSION 0.99.2

* change function name from splitBamByCut to splitGAlignmentsByCut
* add function shiftBam and writeGAlignmentsList
* rewrite splitBam function by just call shiftBam, splitGAlignmentsByCut
and writeGAlignmentsList.

# CHANGES IN VERSION 0.99.1

* fix the bug "object inverted is not exported by 'namespace:Biostrings'"

# CHANGES IN VERSION 0.99.0

* Submit to Bioconductor.

# CHANGES IN VERSION 0.1.0

* Create the package.
