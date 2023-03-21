# ----------------------------------------------------------------
# .R script for renv package to pick up package dependencies
# READ ONLY, NOT NECESSARY TO RUN TO REPRODUCE RESULTS
# ----------------------------------------------------------------

# below libraries are required for Step 1: simulation 1
library(tsDyn)
library(betapart)
library(vegan)
library(entropy)
library(ppcor)
library(doSNOW)
library(parallel)

# below libraries are required for Step 2: simulation 2
library(betapart)
library(vegan)
library(crqa)

# below libraries are required for Step 3: reanalysis
library(betapart)
library(vegan)
library(crqa)
library(esmpack)
library(performance)
library(misty) # for calculating ICC in descriptive statistics
library(nlme) # MLM estimation
library(brms) # MLM bayesian method
library(rstan) # required for brms
