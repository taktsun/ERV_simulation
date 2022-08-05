options(scipen=999)

#======================
# simulation setup
#======================

# study design
siminput <- expand.grid(
  # reps
  rep = 1:2,
  # N
  n = 6,
  # higher meanshift, less rank changes occur, lower ERV expected
  meanshift = c(0.0,2.0),
  # higher auto correlation, lower ERV expected
  autoregressive = c(0.25,0.75),
  #  mean ER endorsement: expects no relationship
  ER_mean = c(0.2,0.3,0.4),
  # higher within-strategy SD, higher ERV expected
  ER_withinSD = c(0.14,0.20,0.26),
  # Number of ER strategies: expects no relationship
  # CJ: Minimum ER to 2 right now, some debugging necessary for 1
  ERn = c(2,3,5),
  # max of scale: expects no relationship
  scalemax = c(100),
  # cross-lagged association: expects no relationship
  cross = c(0),
  # composite metric? (successive comparison included)
  composite = TRUE,
  # measurement correction
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
  # Assign autoregressive regression parameter to diagonal
  Bmat <- diag(autoregressive, ERn)
  # Assign cross-strategy regression coefficient to off-diagonal
  Bmat[Bmat == 0] <- cross
  # Generate data
  out <- VAR.sim(B = Bmat, n = n, lag = 1, include = "none", varcov = diag(ER_withinSD^2, nrow(Bmat)))
  # Center to zero; this is to facilitate your argument ER_mean, but I'm not
  # convinced that this is meaningful
  out <- out - matrix(colMeans(out), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
  # Center to desired mean
  out + matrix(ER_mean, ncol = ncol(out), nrow = nrow(out))
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
  df <- df + get_sign_vector(n = n, ERn = ERn)*meanshift
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
    metric_person_vegan(log1p(df), "euclidean", successive = composite, deco = "norm"),
    metric_person_vegan(df, "chisq", successive = composite),
    metric_person_vegan(df, "euclidean", successive = composite, deco = "hellinger"),
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
