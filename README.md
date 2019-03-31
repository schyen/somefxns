# somefxns

R package. Collection of custom functions that I commonly apply.

The following functions are used for analyzing MIC assays for measuring antibiotic levels in fecal water
* importVictor
* checkCtrl
* plateAdjust
* growthCurve
* survivalCurve
* calcSurvival
* calcabx

The following functions are used for data visualization
* perform_pca
* pca_score
* pca_loading
* pca_biplot
* generate_heatmap

The following functions are used for analyzing Sanger sequences (first iteration of isolationWorkflow)
* abif_fasta
* parse_blast

Other useful functions
* genbymultiple (generate number sequence by multiples. Useful for analysis of dilution series)


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
