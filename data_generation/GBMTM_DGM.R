# =============================================================================
# GBMTM_DGM.R
# -----------------------------------------------------------------------------
# Data-Generating Model (DGM): Group-Based Mixture Trajectory Model (GBMTM)
#
# Generates multivariate longitudinal data from a K-class GBMTM with two
# outcomes (Y, Z), each measured at m repeated time points. Produces n
# independent datasets, exports them to Mplus .dat format, generates Mplus
# .inp files from MplusAutomation templates, runs all models, and extracts
# fit statistics.
#
# Covariance structure:
#   - Within-outcome: diagonal (no serial correlation)
#   - Variance: constant over time
#   - Between-outcome: contemporaneous correlation (equal across classes)
#
# Requirements: MplusAutomation, MASS, stringr
#   install.packages(c("MplusAutomation", "MASS", "stringr"))
#
# Usage:
#   1. Set wd to your working directory (must contain a /Templates subfolder
#      with the MplusAutomation .txt template files).
#   2. Adjust parameters in the CONFIGURATION section as needed.
#   3. Source or run the script.
# =============================================================================

packages <- c("MplusAutomation", "MASS", "stringr")
invisible(lapply(packages, library, character.only = TRUE))


# =============================================================================
# CONFIGURATION: adjust these for your simulation condition
# =============================================================================

wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"   # Must contain /Templates and /Data subfolders

set.seed(2020)

# Simulation dimensions
n  <- 200    # Number of independent replications (datasets)
N  <- 1000   # Sample size per dataset
m  <- 5      # Number of repeated measurement occasions
K  <- 3      # Number of latent classes (same for both outcomes)

# Class proportions (must sum to 1, length K)
pr <- c(0.25, 0.50, 0.25)

# Time points (quadratic trajectory: 0 to 7 over m occasions)
ltime  <- cbind(seq(0, 7, length = m))
ltime2 <- ltime^2
fdm    <- cbind(1, ltime, ltime2)   # Fixed design matrix [intercept, linear, quadratic]

# Fixed effect coefficients per class, outcome Y1
# Rows: [intercept, linear slope, quadratic]; Columns: classes 1..K
bmeanmat <- matrix(4 * c( 0,  0.438, -0.035,   # class 1: rising then falling
                           0,  0.000,  0.000,   # class 2: flat
                           0, -0.438,  0.035),  # class 3: falling then rising
                   nrow = 3, ncol = K,
                   dimnames = list(c("B0","B1","B2"), paste0("k=",1:K)))

# Fixed effect coefficients per class, outcome Y2
bmeanmat2 <- matrix(4 * c(0.214,  0.125,  0.010,
                           0.000, -0.027,  0.002,
                           1.688,  0.125, -0.010),
                    nrow = 3, ncol = K,
                    dimnames = list(c("B0","B1","B2"), paste0("k=",1:K)))

# Within-outcome residual variance (constant over time)
vary1 <- 3.125   # Outcome Y1
vary2 <- 3.160   # Outcome Y2

# Between-outcome (contemporaneous) correlation, shared across classes
rho_between <- 0.45

# =============================================================================
# DERIVED POPULATION PARAMETERS
# =============================================================================

# Model-implied trajectory means per class
fe  <- lapply(1:K, function(k) fdm %*% bmeanmat[, k])   # Y1 means
fe2 <- lapply(1:K, function(k) fdm %*% bmeanmat2[, k])  # Y2 means

residualmean <- rep(0, m)

# Within-outcome covariance matrices (diagonal, constant variance, no serial corr)
make_within_cov <- function(variance) diag(rep(variance, m))

sig_y1 <- lapply(1:K, function(k) make_within_cov(vary1))
sig_y2 <- lapply(1:K, function(k) make_within_cov(vary2))

# Between-outcome (off-diagonal block) covariance, contemporaneous only
sig_between <- lapply(1:K, function(k)
  rho_between * (sqrt(sig_y1[[k]]) %*% sqrt(sig_y2[[k]])))

# Full joint covariance matrix per class (block structure [Y1,Y2])
sig_joint <- lapply(1:K, function(k)
  rbind(cbind(sig_y1[[k]], sig_between[[k]]),
        cbind(sig_between[[k]], sig_y2[[k]])))

# Joint mean vectors per class [Y1 means stacked on Y2 means]
fe_joint <- lapply(1:K, function(k) rbind(fe[[k]], fe2[[k]]))

# =============================================================================
# DATA GENERATION
# =============================================================================

mult_mixture <- function(rep_id) {
  set.seed(rep_id)

  # Draw class membership counts for this replication
  n_per_class <- as.vector(rmultinom(1, N, prob = pr))

  # Generate data for each class and stack
  class_data <- lapply(1:K, function(k) {
    draws <- do.call("rbind",
      replicate(n_per_class[k],
        t(fe_joint[[k]] + mvrnorm(mu = c(residualmean, residualmean),
                                  Sigma = sig_joint[[k]])),
        simplify = FALSE))
    cbind(draws, class_m = k)
  })

  y <- do.call("rbind", class_data)
  colnames(y) <- c("Y1","Y2","Y3","Y4","Y5","Z1","Z2","Z3","Z4","Z5","class_m")

  # Random row order and add subject IDs
  mydata        <- data.frame(repl = rep_id, y)
  mydata        <- mydata[sample(nrow(mydata)), ]
  mydata$id     <- 1:N
  mydata
}

res <- lapply(1:n, mult_mixture)

# =============================================================================
# EXPORT TO MPLUS
# =============================================================================

dir.create(file.path(wd, "Data"), showWarnings = FALSE)

for (i in 1:n) {
  prepareMplusData(as.data.frame(res[[i]]),
                   file.path(wd, "Data", paste0("Mplus_", i, ".dat")),
                   overwrite = TRUE)
}

# Write list of data file names (used by some template approaches)
filenames <- noquote(list.files(file.path(wd, "Data")))
write.table(filenames, file = file.path(wd, "Data", "Data.txt"),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# =============================================================================
# CREATE AND RUN MPLUS MODELS
# =============================================================================

# Generate .inp files from MplusAutomation templates in /Templates subfolder
template_files <- list.files(file.path(wd, "Templates"))
for (j in template_files) {
  createModels(file.path(wd, "Templates", j))
}

# Detect Mplus executable: searches the default install location on Windows;
# assumes 'mplus' is on PATH on Mac/Linux.
.mplus_cmd <- if (.Platform$OS.type == "windows") {
  candidates <- c("C:/Program Files/Mplus/Mplus.exe",
                  "C:/Program Files (x86)/Mplus/Mplus.exe")
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1L] else stop("Mplus not found. Set .mplus_cmd manually.")
} else {
  "mplus"
}

# Run all models; skip existing output unless input file is newer
runModels(wd, recursive = TRUE, showOutput = FALSE,
          replaceOutfile  = "never",
          logFile         = file.path(wd, "ComparisonLog.txt"),
          Mplus_command   = .mplus_cmd)

# =============================================================================
# EXTRACT FIT STATISTICS
# =============================================================================

library(plyr)
allModelParams  <- readModels(wd, recursive = TRUE, what = "summaries")
justSummaries   <- do.call("rbind.fill", sapply(allModelParams, "[", "summaries"))

# Remove incomplete runs (no NGroups = model did not converge/complete)
justSummaries <- justSummaries[complete.cases(justSummaries$NGroups), ]

# String helper functions
left  <- function(text, n) substr(text, 1, n)
mid   <- function(text, start, n) substr(text, start, start + n - 1)
right <- function(text, n) substr(text, nchar(text) - (n - 1), nchar(text))

# Extract replication number and model type from Mplus title string
justSummaries$rep   <- str_remove_all(right(justSummaries$Title, 3), " ")
justSummaries$model <- str_remove_all(
                         str_remove_all(
                           str_remove_all(mid(justSummaries$Title, 9, 10), "rep"), "re"), " ")
justSummaries$class <- as.numeric(str_remove_all(left(justSummaries$Title, 2), "-"))

# Keep relevant columns and convert to factors
df        <- justSummaries[c(2, 6:8, 10:21)]
df$model  <- as.factor(df$model)
df$rep    <- as.factor(df$rep)

save(df, file = file.path(wd, "df.Rda"))
message("Done. Results saved to: ", file.path(wd, "df.Rda"))
