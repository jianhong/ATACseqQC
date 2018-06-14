# CHANGES IN VERSION 1.5.1

* avoid NA in PTscores

# CHANGES IN VERSION 1.3.26

* Add index parameter for fragSizeDist and splitBam.

# CHANGES IN VERSION 1.3.25

* Add warning message for readsDupFreq when bam files without duplicates.

# CHANGES IN VERSION 1.3.24

* Fix a bug in NFRscore.

# CHANGES IN VERSION 1.3.23

* Move IGVSnapshot to extdata because it is not support windows.
* add seqlenghts check for footprints.

# CHANGES IN VERSION 1.3.22

* Try to reduce the memory cost for bamQC.

# CHANGES IN VERSION 1.3.21

* Add doubleCheckDup parameter for bamQC.

# CHANGES IN VERSION 1.3.20

* fix the missing links in documentation.

# CHANGES IN VERSION 1.3.19

* fix the missing links in documentation.

# CHANGES IN VERSION 1.3.18

* Remove the dependence of SRAdb.

# CHANGES IN VERSION 1.3.17

* update vPlot.

# CHANGES IN VERSION 1.3.16

* copy getRelationship from ChIPpeakAnno.

# CHANGES IN VERSION 1.3.15

* add new function distanceDyad.

# CHANGES IN VERSION 1.3.14

* add new function vPlot.

# CHANGES IN VERSION 1.3.13

* fix a bug in footprintsScanner.
* update the vignette.

# CHANGES IN VERSION 1.3.12

* improve the efficiency of bamQC.
* add new function footprintsScanner.

# CHANGES IN VERSION 1.3.11

* update the documentation for function estimateLibComplexity, readsDupFreq, saturationPlot 
* fix a bug in saturationPlot.R ( using sum instead of cumsum for calculate the overall peak breadth)
* improve the efficiency of bamQC.
* add new function IGVSnapshot.

# CHANGES IN VERSION 1.3.10

* add new function plotCorrelation

# CHANGES IN VERSION 1.3.9

* add new functions readFreq, estimateLibComplexity and saturationPlot
* output NRF, PBC1, and PBC2 from bamQC

# CHANGES IN VERSION 1.3.8

* add properPairRate, unmappedRate, hasUnmappedMateRate, 
  notPassingQualityControlsRate in output of bamQC.

# CHANGES IN VERSION 1.3.7

* add mapq summary in output of bamQC.
* add unit test for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.3.6

* add unit test for factorFootprints and fragSizeDist.
  
# CHANGES IN VERSION 1.3.5

* adjust the Cut-site probability by the 
  depth=librarySize/libraryCoverageSize for output of factorFootprints

# CHANGES IN VERSION 1.3.4

* Add function PTscore and NFRscore
* add Profile.segmentation in output of factorFootprints
* add unit test.

# CHANGES IN VERSION 1.3.3

* Fix a bug in factorFootprints when bindingSites is supplied.
* Modified the vignettes.
* Expand the functionality of the bamQC function.
* Import motifStack.

# CHANGES IN VERSION 1.3.2

* Fix a bug in factorFootprints when bindingsite is less than 2.

# CHANGES IN VERSION 1.3.1

* Fix a bug in factorFootprints

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
