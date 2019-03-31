# somefxns

R package. Collection of custom functions that I commonly apply.

The following functions are used for analyzing MIC assays for measuring antibiotic levels in fecal water
1. importVictor
2. checkCtrl
3. plateAdjust
4. growthCurve
5. survivalCurve
6. calcSurvival
7. calcabx

The following functions are used for data visualization
1. perform_pca
2. pca_score
3. pca_loading
4. pca_biplot
5. generate_heatmap

The following functions are used for analyzing Sanger sequences (first iteration of isolationWorkflow)
1. abif_fasta
2. parse_blast

Other useful functions
1. genbymultiple
..* generate number sequence by multiples. Useful for analysis of dilution series


Install Instructions:
```
# if devtools package not yet installed, install devtools first
install_packages('devtools')

# loading devtools library
library(devtools)

# download and install isolationWorkflow
install_github("schyen/somefxns")

# load somefxns package
library('somefxns')
```
