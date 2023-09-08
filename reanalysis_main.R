# These the following package isn't on CRAN. Install manually if needed.
# remotes::install_github("wviechtb/esmpack")

library(betapart)
library(vegan)
library(crqa)
library(esmpack)
library(performance)
library(misty) # for calculating ICC in descriptive statistics
library(nlme) # MLM estimation
library(boot) # for bootstrapping RMSEs

# ===== load functions from other .R ================

# relative SD function
# from: https://ppw.kuleuven.be/okp/software/relative_variability/
source("RSD/checkInput.R")
source("RSD/maximumVAR.R")
source("RSD/relativeSD.R")
source("RSD/checkOutput.R")

# load function to calculate moment-level ER variability indices
source("func_rean_calculateERV.R")
# load function to calculate descriptive statistics
source("func_rean_desStat.R")
# load function to calculate multilevel modeling
source("func_rean_MLM.R")
# load function to bootstrap RMSEs
source("func_bootRMSE.R")


# =======================

# Calculate ER variability indices
    # warnings from calculating relative SD:
    # all-zero ratings return NaN due to division of zero
dfERV1<-loadESMcalculateERV(1)
dfERV2<-loadESMcalculateERV(2)
dfERV3<-loadESMcalculateERV(3)

# Prepare descriptive statistics for supplemental material 5
desstat1 <- summarydesstat(dfERV1)
desstat2 <- summarydesstat(dfERV2)
desstat3 <- summarydesstat(dfERV3)
desstat_all <- rbind((cbind(dataset = 1, desstat1)),
                     (cbind(dataset = 2, desstat2)),
                     (cbind(dataset = 3, desstat3)))

# Calculate multilevel modeling results
    # completeIndices = FALSE: include all available observations (results in manuscript)
    # completeIndices = TRUE : only include observations with no missing/NA variability indices
resMLM1 <- MLMresults(dfERV1,1,completeIndices = FALSE)
resMLM2 <- MLMresults(dfERV2,2,completeIndices = FALSE)
    # modelmoment.betweenRSD threw a warning on singular precision matrix.
    # Fixed effects estimated were almost the same as running lme not using 'optim' optimizer
    # and did not affect the the conclusion of the paper.
resMLM3 <- MLMresults(dfERV3,3,completeIndices = FALSE)
resMLM_all<-rbind(resMLM1,resMLM2,resMLM3)

# output .csv if needed:
# write.csv(desstat_all,"desstat_summary.csv")
# write.csv(resMLM_all,"reanalysis_MLMsummary.csv")

bootrep <- 1000
resbootRMSE <- bootRMSEresults(bootrep = bootrep, saveRdata = TRUE)
# singular precision matrix in some of the bootstrap resamplings but was not critical
# to obtaining RMSE estimate for that occasion of resampling.
