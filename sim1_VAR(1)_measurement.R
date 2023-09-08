library(doSNOW)
library(parallel)
library(ppcor) # partial correlation

options(scipen=999)

#======================
# simulation setup
#======================

# set seed to reproduce exact results
set.seed(1999)

# simulation data generation + ER variability calculation:
# print to text file (TRUE), or store at R environment variable (FALSE)
printtxtresult <- FALSE

# study design
simrep<- 1000

siminput <- expand.grid(
  # reps
  rep = 1:simrep,
  # N
  n = c(30,70,100),
  # autocorrelation
  autoregressive = c(-0.09,0.12,0.33),
  #  mean ER endorsement: now it is used to move the dataset up to avoid negative values
  ER_mean = c(3),
  # within-strategy SD
  ER_withinSD = c(0.10,0.19,0.28),
  # Number of ER strategies
  ERn = c(2,3,5,6),
  # correlation
  correlation = c(-0.11,0.18,0.47),
  missingness = c(0),
  rounding = c(Inf)
)


# ==================================
# define functions
# ==================================
# Load all metric functions
source("func_indices.R")

# correction factor for SD: SD gets inflated when autocorrelation is high,
# So it is adjusted according to Beran (1994). Beran, J. (1994). Statistics for long-memory processes.
correctSD <- function(a,n){
  # see section 1.1 of Beran (1994)
  a<- abs(a)
  d <- 1+2*a/(1-a) # eq. 1.14
  cf <- d*(1-(1/(1-a)/n+ a^n/(1-a)/n)) # eq 1.12
  cf^0.25
}

# Predetermine seed values and store siminput -----------------------------
# remove the rep column
siminput[["rep"]] <- NULL
# correct SD by the Beran formula
siminput[["adjSD"]] <- correctSD(siminput[["autoregressive"]],siminput[["n"]])*siminput[["ER_withinSD"]]

siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))

siminput_missingness <- rbind((replace(siminput, "missingness", 0.1)),
                              (replace(siminput, "missingness", 0.2)),
                              (replace(siminput, "missingness", 0.3)),
                              (replace(siminput, "missingness", 0.4)),
                              (replace(siminput, "missingness", 0.5))
)
siminput_rounding <- rbind((replace(siminput, "rounding", 0)),
                           (replace(siminput, "rounding", 1)),
                           (replace(siminput, "rounding", 2)),
                           (replace(siminput, "rounding", 3)),
                           (replace(siminput, "rounding", 4))
)


siminput <- rbind(siminput, siminput_missingness, siminput_rounding)
saveRDS(siminput, file = "sim1_input_reviseresubmit.RData")


# VAR(1) data generation function
simulate_data <- function(n = 50, ERn = 2, autoregressive = 1, correlation = 0, ER_withinSD = 1, ER_mean = 0, ...){
  # Assign autoregressive regression parameter to diagonal;
  # VAR.sim needs at least 2 columns so +1 if ERn == 1.
  Bmat <- diag(autoregressive, ERn+(ERn==1))
  # create correlation matrix
  Cmat <- diag(1, nrow(Bmat))
  Cmat[Cmat == 0] <- correlation
  # convert into var-covariance matrix.
  Cmat <- Cmat * ER_withinSD^2
  # Generate data
  out <- VAR.sim(B = Bmat, n = n, lag = 1, include = "none", varcov = Cmat)
  # Output a matrix/vector with specified ERn; Set desired mean
  out[,1:ERn, drop = FALSE] + ER_mean
}


# Function to calculate partial correlation between indices and parameters
# for extra analyses in supplemental materials

print_result_measurement <- function(lParam, tabinput){
  nc <- length(lParam)
  lIndex <- colnames(tabinput[,c(8:13)])
  # remove the simulation input parameter columns
  output <- tabinput
  # colnames(output) <- c(lParam,lIndex)

  res_pcor <- data.frame()
  res_pcpvalue <- data.frame()
  res_cor <- data.frame()
  for (i in 1:length(lIndex)){
    dftemp <- cbind(output[,lIndex[i]],output[,1:length(lParam)])
    # tryCatch({
    pcortemp <- round(pcor(dftemp)$estimate[1,],4)
    pcorpvalue <- round(pcor(dftemp)$p.value[1,],4)
    # }, error = function(e){
    #   # assign NA to pcor results as pcor throws error upon NA input
    #   pcortemp <- NA
    #   pcorpvalue <- NA
    # })
    tryCatch({
      cortemp <- round(cor(dftemp)[1,],4)
    }, error = function(e){
      # assign NA to cor results if throwing error
      cortemp <- NA
    })
    res_pcor <- rbind.data.frame(res_pcor,data.frame(as.list(pcortemp)))
    res_pcpvalue  <- rbind.data.frame(res_pcpvalue,data.frame(as.list(pcorpvalue)))
    res_cor <- rbind.data.frame(res_cor,data.frame(as.list(cortemp)))
    res_pcor[i,1] <- lIndex[i]
    res_pcpvalue[i,1] <- lIndex[i]
    res_cor[i,1] <- lIndex[i]
  }
  list(partial_correlation=res_pcor,
       pcor.pvalue=res_pcpvalue,
       correlation = res_cor)
}



#======================
# Simulation: Data generation
#======================


# prepare parallel processing
nclust <- parallel::detectCores()-1
cl <- makeCluster(nclust)
registerDoSNOW(cl)
paths <- .libPaths()
# run simulation
time_start <- Sys.time()
tab <- foreach(rownum = 1:nrow(siminput), .packages = c("tsDyn", "betapart", "vegan", "missMethods"), .combine = rbind) %dopar% {
  # Set seed
  suppressMessages(attach(siminput[rownum, ]))
  set.seed(seed)
  .libPaths(paths)
  # Simulate data
  df <- simulate_data(n = n, ERn = ERn,
                      autoregressive = autoregressive,
                      correlation = correlation,
                      ER_withinSD = ER_withinSD,
                      ER_mean = ER_mean)
  # scalefactor <- scalemax/6
  # df <- df * scalefactor
  df <- delete_MCAR(df, missingness, 1)
  df[!!rowSums(is.na(df)),] <- NA
  df <- round(df, digits = rounding)
  # df[df<0] <- 0
  # df[df>scalemax] <- scalemax
  tempbray.suc <- calc.bray.suc(df)
  tempbray.amm <- calc.bray.amm(df)
  out <- c(
    Autoregression = autoregressive,
    withinStgy_SD = adjSD,
    Correlation = correlation,
    NER = ERn,
    nobs = n,
    missingness = missingness,
    rounding = rounding,
    sucf = tempbray.suc[3],        # Bray-Curtis dissimilarity: Full index
    sucr = tempbray.suc[1],        # Bray-Curtis dissimilarity: replacement
    sucn = tempbray.suc[2],        # Bray-Curtis dissimilarity: nestedness
    ammf = tempbray.amm[3],        # Bray-Curtis dissimilarity: Full index
    ammr = tempbray.amm[1],        # Bray-Curtis dissimilarity: replacement
    ammn = tempbray.amm[2]        # Bray-Curtis dissimilarity: nestedness
    # metric_person_within_SD(df),
    # metric_person_between_SD(df),
    # metric_person_mssd(df),
    # metric_person_sdSD(df),
    # metric_person_SD(df),
    # metric_person_beta(df, "bray"),        # Bray-Curtis dissimilarity: Full index
    # metric_person_beta(df, "bray",".bal"), # Bray-Curtis dissimilarity: Replacement subscomponent
    # metric_person_beta(df, "bray",".gra"), # Bray-Curtis dissimilarity: Nestedness subcomponent
    # metric_person_beta(df, "ruz"),         # Jaccard dissimilarity: Full index
    # metric_person_beta(df, "ruz",".bal"),  # Jaccard dissimilarity: Replacement subscomponent
    # metric_person_beta(df, "ruz",".gra"),  # Jaccard dissimilarity: Nestedness subcomponent
    # metric_person_vegan(df, "chord"),
    # metric_person_vegan(df, "chisq")
  )

  if (printtxtresult){
  # Option 1: print text files
  write.table(x = t(out), file = sprintf("simresults_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  NULL
  }else{
  # Option 2: print to R environment variable
  out
  }
}
end_time <- Sys.time()
#Close cluster
stopCluster(cl)

time_per_row <- as.numeric(end_time - time_start) / nrow(siminput)
writeLines(paste0("Time per row: ", time_per_row), "time_per_row.txt")

# End of simulation -------------------------------------------------------
print("End of simulation")

#======================
# Simulation: Evaluate indices' performance
#======================

# if results are stored at .txt, read and combine the .txt output
if (is.null(tab)){
  txt_files_ls = list.files(pattern="simresults_*")
  txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = F, sep ="\t")})
  tab <-do.call("rbind", lapply(txt_files_df, as.matrix))
}

# the below list is in the same order as what was inserted each row during loop
list_parameters <- c("Autoregression",
                     "withinStgy_SD",
                     "Correlation",
                     "NER",
                     "nobs"
                     )

res.missingness.00 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0)
                                                   ,])
res.missingness.01 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0.1)
                                                   ,])
res.missingness.02 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0.2)
                                                   ,])
res.missingness.03 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0.3)
                                                   ,])
res.missingness.04 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0.4)
                                                   ,])
res.missingness.05 <- print_result_measurement(list_parameters,
                                               tab[(tab[,"rounding"] == Inf &
                                                      tab[,"missingness"] == 0.5)
                                                   ,])
res.rounding.Inf <- print_result_measurement(list_parameters,
                                              tab[(tab[,"rounding"] == Inf &
                                                     tab[,"missingness"] == 0)
                                                  ,])
res.rounding.0 <- print_result_measurement(list_parameters,
                                              tab[(tab[,"rounding"] == 0 &
                                                     tab[,"missingness"] == 0)
                                                  ,])
res.rounding.1 <- print_result_measurement(list_parameters,
                                           tab[(tab[,"rounding"] == 1 &
                                                  tab[,"missingness"] == 0)
                                               ,])
res.rounding.2 <- print_result_measurement(list_parameters,
                                           tab[(tab[,"rounding"] == 2 &
                                                  tab[,"missingness"] == 0)
                                               ,])
res.rounding.3 <- print_result_measurement(list_parameters,
                                           tab[(tab[,"rounding"] == 3 &
                                                  tab[,"missingness"] == 0)
                                               ,])
res.rounding.4 <- print_result_measurement(list_parameters,
                                           tab[(tab[,"rounding"] == 4 &
                                                  tab[,"missingness"] == 0)
                                               ,])
res.missingness.00
res.missingness.01
res.missingness.02
res.missingness.03
res.missingness.04
res.missingness.05
res.rounding.Inf
res.rounding.4
res.rounding.3
res.rounding.2
res.rounding.1
res.rounding.0

write.csv(res.missingness.00,paste0("sim1_measurement_m0rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.missingness.01,paste0("sim1_measurement_m1rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.missingness.02,paste0("sim1_measurement_m2rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.missingness.03,paste0("sim1_measurement_m3rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.missingness.04,paste0("sim1_measurement_m4rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.missingness.05,paste0("sim1_measurement_m5rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.Inf,paste0("sim1_measurement_m0rI_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.4,paste0("sim1_measurement_m0r4_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.3,paste0("sim1_measurement_m0r3_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.2,paste0("sim1_measurement_m0r2_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.1,paste0("sim1_measurement_m0r1_",simrep," ",Sys.Date(),".csv"))
write.csv(res.rounding.0,paste0("sim1_measurement_m0r0_",simrep," ",Sys.Date(),".csv"))
