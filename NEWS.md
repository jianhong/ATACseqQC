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
