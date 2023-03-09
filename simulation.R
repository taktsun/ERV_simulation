options(scipen=999)

#======================
# simulation setup
#======================

# study design
siminput <- expand.grid(
  # reps
  rep = 1:100,
  # N
  n = c(30,70,100),
  # meanshift in multiples of SD
  meanshift = c(0.0),
  # autocorrelation
  autoregressive = c(-0.09,0.12,0.33),
  #  mean ER endorsement: now it is used to move the dataset up to avoid -ve values
  ER_mean = c(3),
  # within-strategy SD
  ER_withinSD = c(0.10,0.19,0.28),
  # Number of ER strategies
  ERn = c(2,3,5,6),
  # max of scale: expects no relationship
  scalemax = c(100),
  # correlation
  correlation = c(-0.11,0.18,0.47),
  # cross-lagged association
  cross = c(0),
  # measurement correction
  measurementcorrection = FALSE
)

# Load all correction functions
source("correction_functions.R")

siminput[["rep"]] <- NULL
siminput[["adjSD"]] <- correctSD(siminput[["autoregressive"]],siminput[["n"]])*siminput[["ER_withinSD"]]
# Predetermine seed values and store siminput -----------------------------

set.seed(1999)
siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))
saveRDS(siminput, file = "siminput.RData")

# ==================================
# define functions
# ==================================
simulate_data <- function(n = 50, ERn = 2, autoregressive = 1, correlation = 0, cross = 0, ER_withinSD = 1, ER_mean = 0, ...){
  # Assign autoregressive regression parameter to diagonal;
  # VAR.sim needs at least 2 columns so +1 if ERn == 1.
  Bmat <- diag(autoregressive, ERn+(ERn==1))
  # Assign cross-strategy regression coefficient to off-diagonal
  Bmat[Bmat == 0] <- cross
  # create correlation matrix
  Cmat <- diag(1, nrow(Bmat))
  Cmat[Cmat == 0] <- correlation
  # convert into var-covariance matrix.
  Cmat <- Cmat * ER_withinSD^2
  # Generate data
  out <- VAR.sim(B = Bmat, n = n, lag = 1, include = "none", varcov = Cmat)
  # Center to zero; this is to facilitate your argument ER_mean, but I'm not
  # convinced that this is meaningful
  # EL: indeed this is not needed. The VAR model has zero mean with sufficient number of repetitions.
  # out <- out - matrix(colMeans(out), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
  # Output a matrix/vector with specified ERn; Set desired mean
  out[,1:ERn, drop = FALSE] + ER_mean
}
# Load all metric functions
source("metric_functions.R")

# prepare parallel processing
library(doSNOW)
library(parallel)
nclust <- parallel::detectCores()
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
  df <- simulate_data(n = n, ERn = ERn, autoregressive = autoregressive,
                      correlation = correlation, cross = cross,
                      ER_withinSD = ER_withinSD, ER_mean = ER_mean)
  # Mean shift
  df <- df + get_sign_vector(n = n, ERn = ERn)*meanshift*ER_withinSD
  # Measurement corrections
  if (measurementcorrection){
    # CJ: THis introduces bias, that's probably not what you want to do
    df <- correct_range_bound(df)
  }

  out <- c(
    beta.multi.abund(df)$beta.BRAY,
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
    metric_person_vegan(df, "chisq"),
    metric_person_KLdiv(df)
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
