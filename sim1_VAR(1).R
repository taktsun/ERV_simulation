options(scipen=999)

#======================
# simulation setup
#======================

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

# correction factor for SD
correctSD <- function(a,n){
  # see section 1.1 of Beran (1994)
  a<- abs(a)
  d <- 1+2*a/(1-a) # eq. 1.14
  cf <- d*(1-(1/(1-a)/n+ a^n/(1-a)/n)) # eq 1.12
  cf^0.25
}

# remove the rep column
siminput[["rep"]] <- NULL
# correct SD by the Beran formula
siminput[["adjSD"]] <- correctSD(siminput[["autoregressive"]],siminput[["n"]])*siminput[["ER_withinSD"]]

# Predetermine seed values and store siminput -----------------------------
set.seed(1999)
siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))
saveRDS(siminput, file = "siminput.RData")

# ==================================
# define functions
# ==================================
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
# Load all metric functions
source("metric_functions.R")

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
  # CJ: If the sim takes very long or uses a lot of memory, print to text files.
  # CJ: If not, just return the output.
  # CJ: We have to see which works best when we start running a larger chunk!
  # write.table(x = t(out), file = sprintf("results_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  # NULL
  out
}
end_time <- Sys.time()
#Close cluster
stopCluster(cl)

time_per_row <- as.numeric(end_time - time_start) / nrow(siminput)
writeLines(paste0("Time per row: ", time_per_row), "time_per_row.txt")

# End of simulation -------------------------------------------------------
stop("End of simulation")

#======================
# evaluate indices' performance
#======================


library(ppcor) # partial correlation
source("func_sim1performance.R")

n_singleoutput_metric <- 2
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
# allmoment, successive, or composite
print_result("successive")
