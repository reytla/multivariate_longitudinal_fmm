# =============================================================================
# MCPGMM_DGM.R
# -----------------------------------------------------------------------------
# Data-Generating Model (DGM): Multivariate Covariance Pattern Growth Mixture
# Model (MCPGMM)
#
# Identical pipeline to GBMTM_DGM.R. Key differences vs. GBMTM/MLCGA:
#
#   1. Within-outcome serial correlation: Toeplitz structure
#      (correlations decay with lag: 0.45, 0.40, 0.35, 0.30)
#   2. Time-varying residual variance (same as MLCGA)
#   3. Between-outcome covariance uses a diagonal-restricted formula:
#      rho * diag(sqrt(sigma_y1)) %*% diag(sqrt(sigma_y2))
#      rather than the full matrix product used in GBMTM/MLCGA
#
# Requirements: MplusAutomation, MASS, stringr
# Usage: see GBMTM_DGM.R header
# =============================================================================

packages <- c("MplusAutomation", "MASS", "stringr")
invisible(lapply(packages, library, character.only = TRUE))


# =============================================================================
# CONFIGURATION: adjust these for your simulation condition
# =============================================================================

wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"

set.seed(2020)

n  <- 200
N  <- 1000
m  <- 5
K  <- 3

pr <- c(0.25, 0.50, 0.25)

ltime  <- cbind(seq(0, 7, length = m))
ltime2 <- ltime^2
fdm    <- cbind(1, ltime, ltime2)

bmeanmat <- matrix(4 * c( 0,  0.438, -0.035,
                           0,  0.000,  0.000,
                           0, -0.438,  0.035),
                   nrow = 3, ncol = K,
                   dimnames = list(c("B0","B1","B2"), paste0("k=",1:K)))

bmeanmat2 <- matrix(4 * c(0.214,  0.125,  0.010,
                           0.000, -0.027,  0.002,
                           1.688,  0.125, -0.010),
                    nrow = 3, ncol = K,
                    dimnames = list(c("B0","B1","B2"), paste0("k=",1:K)))

# Time-varying variance (same as MLCGA)
vary1_tv <- c(2.025, 2.575, 3.125, 3.675, 4.225)
vary2_tv <- c(2.060, 2.610, 3.160, 3.710, 4.260)

# Within-outcome serial (Toeplitz) correlation, shared across classes
# Lag-1: 0.45, Lag-2: 0.40, Lag-3: 0.35, Lag-4: 0.30
rho_within <- toeplitz(c(1, 0.45, 0.40, 0.35, 0.30))

# Between-outcome (contemporaneous) correlation
rho_between <- 0.45

# =============================================================================
# DERIVED POPULATION PARAMETERS
# =============================================================================

fe  <- lapply(1:K, function(k) fdm %*% bmeanmat[, k])
fe2 <- lapply(1:K, function(k) fdm %*% bmeanmat2[, k])
residualmean <- rep(0, m)

# Within-outcome covariance: Toeplitz correlation scaled by time-varying SD
make_cov_toeplitz <- function(variance_vec, rho_mat) {
  # Standard formula: D_sd %*% R %*% D_sd
  d_sd <- diag(sqrt(variance_vec))
  d_sd %*% rho_mat %*% d_sd
}

sig_y1 <- lapply(1:K, function(k) make_cov_toeplitz(vary1_tv, rho_within))
sig_y2 <- lapply(1:K, function(k) make_cov_toeplitz(vary2_tv, rho_within))

# Between-outcome: contemporaneous only; Cov(Y_t, Z_t) = rho * sd_y(t) * sd_z(t)
# sig_y1 is a full Toeplitz matrix, so marginal SDs must be taken from vary1_tv
# directly (not from sqrt(sig_y1), which is a matrix square root, not element-wise)
sig_between <- lapply(1:K, function(k)
  rho_between * diag(sqrt(vary1_tv) * sqrt(vary2_tv)))

sig_joint <- lapply(1:K, function(k)
  rbind(cbind(sig_y1[[k]], sig_between[[k]]),
        cbind(sig_between[[k]], sig_y2[[k]])))

fe_joint <- lapply(1:K, function(k) rbind(fe[[k]], fe2[[k]]))

# =============================================================================
# DATA GENERATION
# =============================================================================

mult_mixture <- function(rep_id) {
  set.seed(rep_id)
  n_per_class <- as.vector(rmultinom(1, N, prob = pr))

  class_data <- lapply(1:K, function(k) {
    draws <- do.call("rbind",
      replicate(n_per_class[k],
        t(fe_joint[[k]] + mvrnorm(mu = c(residualmean, residualmean),
                                  Sigma = sig_joint[[k]])),
        simplify = FALSE))
    cbind(draws, class_m = k)
  })

  y             <- do.call("rbind", class_data)
  colnames(y)   <- c("Y1","Y2","Y3","Y4","Y5","Z1","Z2","Z3","Z4","Z5","class_m")
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

filenames <- noquote(list.files(file.path(wd, "Data")))
write.table(filenames, file = file.path(wd, "Data", "Data.txt"),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# =============================================================================
# CREATE AND RUN MPLUS MODELS
# =============================================================================

template_files <- list.files(file.path(wd, "Templates"))
for (j in template_files) {
  createModels(file.path(wd, "Templates", j))
}

.mplus_cmd <- if (.Platform$OS.type == "windows") {
  candidates <- c("C:/Program Files/Mplus/Mplus.exe",
                  "C:/Program Files (x86)/Mplus/Mplus.exe")
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1L] else stop("Mplus not found. Set .mplus_cmd manually.")
} else {
  "mplus"
}

runModels(wd, recursive = TRUE, showOutput = FALSE,
          replaceOutfile = "never",
          logFile        = file.path(wd, "ComparisonLog.txt"),
          Mplus_command  = .mplus_cmd)

# =============================================================================
# EXTRACT FIT STATISTICS
# =============================================================================

library(plyr)
allModelParams <- readModels(wd, recursive = TRUE, what = "summaries")
justSummaries  <- do.call("rbind.fill", sapply(allModelParams, "[", "summaries"))
justSummaries  <- justSummaries[complete.cases(justSummaries$NGroups), ]

left  <- function(text, n) substr(text, 1, n)
mid   <- function(text, start, n) substr(text, start, start + n - 1)
right <- function(text, n) substr(text, nchar(text) - (n - 1), nchar(text))

justSummaries$rep   <- str_remove_all(right(justSummaries$Title, 3), " ")
justSummaries$model <- str_remove_all(
                         str_remove_all(
                           str_remove_all(mid(justSummaries$Title, 9, 10), "rep"), "re"), " ")
justSummaries$class <- as.numeric(str_remove_all(left(justSummaries$Title, 2), "-"))

df       <- justSummaries[c(2, 6:8, 10:21)]
df$model <- as.factor(df$model)
df$rep   <- as.factor(df$rep)

save(df, file = file.path(wd, "df.Rda"))
message("Done. Results saved to: ", file.path(wd, "df.Rda"))
