options(scipen=999)

#======================
# simulation setup
#======================

# set seed to reproduce exact results
set.seed(1999)

# study design
siminput <- expand.grid(
  # reps
  rep = 1:1000,
  # N
  n = c(30,70,100),
  # autocorrelation
  autoregressive = c(-0.09,0.12,0.33),
  #  mean ER endorsement: now it is used to move the dataset up to avoid -ve values
  ER_mean = c(3),
  # within-strategy SD
  ER_withinSD = c(0.10,0.19,0.28),
  # Number of ER strategies]
  ERn = c(2,3,5,6),
  # max of scale: expects no relationship
  scalemax = c(100),
  # correlation
  correlation = c(-0.11,0.18,0.47)
)

# ==================================
# define functions
# ==================================
# Load all metric functions
source("func_indices.R")

# correction factor for SD
correctSD <- function(a,n){
  # see section 1.1 of Beran (1994)
  a<- abs(a)
  d <- 1+2*a/(1-a) # eq. 1.14
  cf <- d*(1-(1/(1-a)/n+ a^n/(1-a)/n)) # eq 1.12
  cf^0.25
}

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

print_result <- function(timemode="successive", lParam, lIndex){
  nc <- length(lParam)+n_singleoutput_metric
  # remove the simulation input parameter columns
  output <- (tab[,-c(1:nc)])

  # 0 = suc.second, 1 = mom.all
  if (timemode == "successive"){
    output <- output[, seq_len(ncol(output)) %% 3 == 0]
  } else { #"allmoment"
    output <- output[, seq_len(ncol(output)) %% 3 == 1]
  }
  output <- cbind((tab[,  1:nc]),output)
  colnames(output) <- c(lParam,lIndex)

  res_pcor <- data.frame()
  res_pcpvalue <- data.frame()
  for (i in 1:length(lIndex)){
    dftemp <- cbind(output[,lIndex[i]],output[,1:length(lParam)])
    tryCatch({
      pcortemp <- round(pcor(dftemp)$estimate[1,],4)
      pcorpvalue <- round(pcor(dftemp)$p.value[1,],4)
    }, error = function(e){
      # assign NA to pcor results as pcor throws error upon NA input
      pcortemp <- NA
      pcorpvalue <- NA
    })
    res_pcor <- rbind.data.frame(res_pcor,data.frame(as.list(pcortemp)))
    res_pcpvalue  <- rbind.data.frame(res_pcpvalue,data.frame(as.list(pcorpvalue)))
    res_pcor[i,1] <- lIndex[i]
    res_pcpvalue[i,1] <- lIndex[i]
  }
  list(partial_correlation=res_pcor,pcor.pvalue=res_pcpvalue)
}



#======================
# Simulation: Data generation
#======================

# remove the rep column
siminput[["rep"]] <- NULL
# correct SD by the Beran formula
siminput[["adjSD"]] <- correctSD(siminput[["autoregressive"]],siminput[["n"]])*siminput[["ER_withinSD"]]

# Predetermine seed values and store siminput -----------------------------
siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))
saveRDS(siminput, file = "siminput.RData")


# prepare parallel processing
library(doSNOW)
library(parallel)
nclust <- parallel::detectCores()-1
cl <- makeCluster(nclust)
registerDoSNOW(cl)
paths <- .libPaths()
# run simulation
time_start <- Sys.time()
# EL: not sure why adding "paths = paths," bring error in my machine. Temporarily remove
tab <- foreach(rownum = 1:nrow(siminput), .packages = c("tsDyn", "betapart", "vegan", "entropy"), .combine = rbind) %dopar% {
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

  out <- c(
    autoregressive,
    adjSD,
    correlation,
    ERn,
    n,
    metric_person_within_SD(df),
    metric_person_between_SD(df),
    metric_person_SD(df),
    metric_person_beta(df, "bray"),
    metric_person_beta(df, "bray",".bal"),
    metric_person_beta(df, "bray",".gra"),
    metric_person_beta(df, "ruz"),
    metric_person_beta(df, "ruz",".bal"),
    metric_person_beta(df, "ruz",".gra"),
    metric_person_vegan(df, "chord"),
    metric_person_vegan(df, "chisq")
  )

  # Option 1: print text files
  write.table(x = t(out), file = sprintf("results_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  NULL

  # Option 2: print to R environment variable
  # out
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

if (is.null(tab)){
  txt_files_ls = list.files(pattern="result*")
  txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = F, sep ="\t")})
  tab <-do.call("rbind", lapply(txt_files_df, as.matrix))
}

library(ppcor) # partial correlation
# source("func_sim1performance.R")
source("func_sim1performance_new.R")



n_singleoutput_metric <- 2
# the below list is in the same order as what was inserted each row during loop
list_parameters <- c("Autoregression",
                     "withinStgy_SD",
                     "Correlation",
                     "NER",
                     "nobs")
# the below list has to match what specified in sim1_VAR(1).R
list_metrics <- c("withinSD",
                  "betweenSD",
                  # person-level/non-temporal indices go above
                  # remember to amend n_singleoutput_metric = number of such indices
                  "SD",
                  "bray",
                  "bray.bal",
                  "bray.gra",
                  "jaccard",
                  "jaccard.bal",
                  "jaccard.gra",
                  "chord",
                  "chisq"
)

# allmoment OR successive
res.suc <- print_result("successive", list_parameters, list_metrics)
res.amm <- print_result("allmoment", list_parameters, list_metrics)

# print results - successive difference
res.suc
# print results - all-moment comparison
res.amm
