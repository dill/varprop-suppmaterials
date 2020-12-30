#  Supplementary materials for the paper "Variance propagation for density surface models"

Two examples from the paper:

A. Re-analysis of the Island Scrub-jay survey from Sillet et al. (2012)
B. Group-size factor-smooth interaction for harbour porpoise from Hammond et al (2013).

## Data

Island Scrub-jay data is from the paper:

*Sillett, T. Scott, Richard B. Chandler, J. Andrew Royle, Marc Kery, and Scott A. Morrison. ‘Hierarchical Distance-Sampling Models to Estimate Population Size and Habitat-Specific Abundance of an Island Endemic’. Ecological Applications 22, no. 7 (2012): 1997-2006*

Data is available at: https://figshare.com/articles/Supplement_1_R_code_data_and_grid_covariates_ used_in_the_analyses_/3517754

Harbour porpoise data are from the SCANS-II project supported by the EU LIFE Nature programme  under project LIFE04NAT/GB/000245 and by the governments of range states:  Belgium, Denmark, France, Germany, Ireland, Netherlands, Norway, Poland, Portugal, Spain, Sweden and UK. Thanks to Louise Burt and Phil Hammond for assistance with the data.

## Software

We recommend using the latest versions of `dsm`, `mrds` and `Distance` available from github, using the following commands:

```r
#install.packages("devtools")
devtools::install.github("DistanceDevelopment/mrds")
devtools::install.github("DistanceDevelopment/Distance")
devtools::install.github("DistanceDevelopment/dsm")
```

