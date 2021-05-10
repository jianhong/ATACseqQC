# CHANGES IN VERSION 1.15.11

* Break the limitation of sequence length must have ends less than or equal to 536870912.

# CHANGES IN VERSION 1.15.10

* fix the issue that idxstatsBam return values with "*"

# CHANGES IN VERSION 1.15.9

* Add rmarkdown as suggest package.

# CHANGES IN VERSION 1.15.8

* update documentation for the case when no BSgenome object is available.

# CHANGES IN VERSION 1.15.7

* fix the NA values for TSSEscore when infinite value is in the data.

# CHANGES IN VERSION 1.15.6

* fix the missing link of documentation for rtracklyaer:import.
    
# CHANGES IN VERSION 1.15.5

* remove duplicates when shift reads.

# CHANGES IN VERSION 1.15.4

* Fix the issue when empty object input into exportBamFile.

# CHANGES IN VERSION 1.15.3

* Reuse header when exportBamFile in splitGAlignmentsByCut function.

# CHANGES IN VERSION 1.15.2

* Fix the tag MC in exportBamFile function.

# CHANGES IN VERSION 1.15.1

* write exportBamFile function to replace rtracklayer::export.bam.

# CHANGES IN VERSION 1.13.9

* fix the issue that plotCorrelation heatmap is scaled by row.

# CHANGES IN VERSION 1.13.8

* throw an error if not enought nucleosome free read nor mononucleosome reads for training.

# CHANGES IN VERSION 1.13.7

* fix a bug introduced by matchPWM by paste ^ and $ into exclude sequence name.

# CHANGES IN VERSION 1.13.6

* update documentation of plotFootprints.

# CHANGES IN VERSION 1.13.4

* fix a formular for TSSE score.

# CHANGES IN VERSION 1.13.3

* fix a bug the after shift, the index is not changed.

# CHANGES IN VERSION 1.13.2

* change the normalization method by library size for factorFootprints for user defined group samples.

# CHANGES IN VERSION 1.13.1

* Add documentation to decribe the format of bindingSites for factorFootprints.

# CHANGES IN VERSION 1.11.9

* Fix a issue for splitGAlignmentsList, splitGAlignments by supply the header info to mergeBam.

# CHANGES IN VERSION 1.11.8

* Fix a issue 'mergeBam' 'destination' exists, 'overwrite' is FALSE for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.11.7

* fix a typo in doc for readBamFile.

# CHANGES IN VERSION 1.11.6

* fix the sample for shiftGAlignments.

# CHANGES IN VERSION 1.11.5

* add function shiftGAlignments for single end reads.

# CHANGES IN VERSION 1.11.4

* Fix the issue of if there is no reads in bam file for shiftGAlignmentsList.

# CHANGES IN VERSION 1.11.2

* Fix the issue of "[E::sam_parse1] unrecognized type N".

# CHANGES IN VERSION 1.11.1

* Add flag parameter for splitBam.

# CHANGES IN VERSION 1.9.9

* export prepareBindingSitesList function.
* Add rownames for footprintsScanner counts data.

# CHANGES IN VERSION 1.9.8

* Add error message for vPlot when no paired reads in bam file.

# CHANGES IN VERSION 1.9.7

* Fix the bug that gscore changed the output for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.9.6

* Try to decrease the memory cost for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.9.5

* Try to decrease the memory cost for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.9.4

* Add the error handle if not enough mononucleosome reads for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.9.3

* Try to decrease the memory cost for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.9.2

* Fix the bug if the bam file containsupplementary alignments.

# CHANGES IN VERSION 1.9.1

* Fix the bug if the bam file contain mix of single ends and paired ends.

# CHANGES IN VERSION 1.7.9

* add maximalBindingWidth parameter to footprintsScanner

# CHANGES IN VERSION 1.7.8

* change the function of footprintsScanner.

# CHANGES IN VERSION 1.7.7

* add parameter outPath for splitGAlignmentsByCut.

# CHANGES IN VERSION 1.7.6

* fix a bug in PTscore.

# CHANGES IN VERSION 1.7.5

* use `file.copy` instead of `file.rename`.

# CHANGES IN VERSION 1.7.4

* add parameter outbam for shiftGAlignmentsList.

# CHANGES IN VERSION 1.7.3

* Update documentation for Transcription start site (TSS) enrichment values

# CHANGES IN VERSION 1.7.2

* add the new biocViews tag "ImmunoOncology"

# CHANGES IN VERSION 1.7.1

* add more documentation for bigFile parameter in readBamFile.R

# CHANGES IN VERSION 1.5.7

* fix the memory issue of big bam file.

# CHANGES IN VERSION 1.5.6

* fix the bug when the reads length is smaller than 5 for shiftGAlignmentsList

# CHANGES IN VERSION 1.5.4

* export plotFootprints.

# CHANGES IN VERSION 1.5.3

* add Feng Yan as an author for function of estimateLibComplexity.

# CHANGES IN VERSION 1.5.2

* replace 'ds.mincount.bootstrap' with 'ds.rSAC.bootstrap'
* add Transcription Start Site (TSS) Enrichment Score: TSSEscore

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
