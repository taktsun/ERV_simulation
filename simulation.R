options(scipen=999)

#======================
# simulation setup
#======================

# study design
siminput <- expand.grid(
  # reps
  rep = 1:20,
  # N
  n = c(10,100),
  # meanshift in multiples of SD. High meanshift, lower ERV expected
  meanshift = c(0.0,1,2),
  # higher auto correlation, lower ERV expected
  autoregressive = c(-0.5,0.25,0.75),
  #  mean ER endorsement: only accept values between 0 - 1
  ER_mean = c(0.2,0.3,0.4),
  # within-strategy SD: only accept values between 0 - 1
  ER_withinSD = c(0.14,0.20,0.26),
  # Number of ER strategies: expects no relationship
  ERn = c(2,3,5),
  # max of scale: expects no relationship
  scalemax = c(100),
  # cross-lagged association: expects no relationship
  cross = c(0),
  # composite metric? (include successive comparison?)
  composite = TRUE,
  # measurement correction?
  measurementcorrection = TRUE
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
simulate_data <- function(n = 50, ERn = 2, autoregressive = 1, cross = 0, ER_withinSD = 1, ER_mean = 0, ...){
  # Assign autoregressive regression parameter to diagonal;
  # VAR.sim needs at least 2 columns so +1 if ERn == 1.
  Bmat <- diag(autoregressive, ERn+(ERn==1))
  # Assign cross-strategy regression coefficient to off-diagonal
  Bmat[Bmat == 0] <- cross
  # Generate data
  out <- VAR.sim(B = Bmat, n = n, lag = 1, include = "none", varcov = diag(ER_withinSD^2, nrow(Bmat)))
  # Center to zero; this is to facilitate your argument ER_mean, but I'm not
  # convinced that this is meaningful
  out <- out - matrix(colMeans(out), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
  # Output a matrix/vector with specified ERn; Set desired mean
  out[,1:ERn] + ER_mean
}
# Load all metric functions
source("metric_functions.R")

# prepare parallel processing
library(doSNOW)
library(parallel)
nclust <- parallel::detectCores()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

# run simulation
time_start <- Sys.time()
tab <- foreach(rownum = 1:nrow(siminput), .packages = c("tsDyn", "betapart", "vegan", "entropy"), .combine = rbind) %dopar% {
  # Set seed
  suppressMessages(attach(siminput[rownum, ]))
  set.seed(seed)

  # Simulate data
  df <- simulate_data(n = n, ERn = ERn, autoregressive = autoregressive, cross = cross, ER_withinSD = ER_withinSD, ER_mean = ER_mean)
  # Mean shift
  df <- df + get_sign_vector(n = n, ERn = ERn)*meanshift*ER_withinSD
  # Measurement corrections
  if (measurementcorrection){
    df <- correct_range_bound(df)
  }

  out <- c(
    metric_mssd(df),
    metric_mean_euclidean(df),
    metric_person_between_SD(df),
    metric_person_within_SD(df),
    metric_person_SD(df,successive = composite),
    metric_person_vegan(df, "euclidean", successive = composite),
    metric_person_vegan(df, "manhattan", successive = composite),
    metric_person_vegan(df, "chord", successive = composite),
    metric_person_vegan(decostand(log1p(df),"norm"), "euclidean", successive = composite),
    metric_person_vegan(df, "chisq", successive = composite),
    metric_person_vegan(decostand(df,"hellinger"), "euclidean", successive = composite),
    metric_person_vegan(df, "jaccard", successive = composite),
    metric_person_vegan(df, "kulczynski", successive = composite),
    metric_person_vegan(df, "bray", successive = composite),
    metric_person_KLdiv(df, successive = composite),
    beta.multi.abund(df)$beta.BRAY

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
