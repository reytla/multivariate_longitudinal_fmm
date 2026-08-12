# Set working directory to project root when run via Rscript
.file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg))
  setwd(normalizePath(file.path(dirname(sub("--file=", "", .file_arg[1L])), "..", ".."),
                      winslash = "/"))
source("ovl/ovl_trajectories.R")

pass <- 0L
fail <- 0L
check <- function(label, cond) {
  if (isTRUE(cond)) {
    cat("[PASS]", label, "\n")
    pass <<- pass + 1L
  } else {
    cat("[FAIL]", label, "\n")
    fail <<- fail + 1L
  }
}

# Shared fixtures: from ovl_trajectories.R worked example and GBMTM_DGM.R
t_vec <- seq(0, 7, length = 5)        # DGM time points
X     <- cbind(1, t_vec, t_vec^2)

beta1 <- 4 * c( 0,  0.438, -0.035)
beta2 <- 4 * c( 0,  0.000,  0.000)
beta3 <- 4 * c( 0, -0.438,  0.035)

means_y1 <- list(
  Class_1 = X %*% beta1,
  Class_2 = X %*% beta2,
  Class_3 = X %*% beta3
)
cov_y1   <- 3.125 * diag(5)
covs_y1  <- list(cov_y1, cov_y1, cov_y1)

# -----------------------------------------------------------------------
# T1: Single outcome returns a matrix, not a list
# -----------------------------------------------------------------------
ovl1 <- compute_ovl(means_y1, covs_y1)
check("T1  Single outcome returns a matrix", is.matrix(ovl1))

# -----------------------------------------------------------------------
# T2: Diagonal = 1 (OVL of a class with itself)
# -----------------------------------------------------------------------
check("T2  Diagonal entries are all 1", all(diag(ovl1) == 1))

# -----------------------------------------------------------------------
# T3: Off-diagonal in [0, 1]
# -----------------------------------------------------------------------
off <- ovl1[upper.tri(ovl1)]
check("T3  Off-diagonal OVL values in [0, 1]",
      all(off >= 0) && all(off <= 1))

# -----------------------------------------------------------------------
# T4: Matrix is symmetric
# -----------------------------------------------------------------------
check("T4  OVL matrix is symmetric",
      all(ovl1 == t(ovl1)))

# -----------------------------------------------------------------------
# T5: Identical classes -> OVL = 1
# -----------------------------------------------------------------------
means_id <- list(Class_1 = X %*% beta1,
                 Class_2 = X %*% beta1)   # both classes identical
covs_id  <- list(cov_y1, cov_y1)
ovl_id   <- compute_ovl(means_id, covs_id)
check("T5  Identical classes give OVL = 1",
      abs(ovl_id[1, 2] - 1) < 1e-6)

# -----------------------------------------------------------------------
# T6: Well-separated classes -> OVL close to 0
# -----------------------------------------------------------------------
# Means 30 SDs apart: negligible overlap
means_sep <- list(Class_1 = rep(0,  5),
                  Class_2 = rep(100, 5))
cov_tiny  <- diag(5)    # SD = 1
ovl_sep   <- compute_ovl(means_sep, list(cov_tiny, cov_tiny))
check("T6  Well-separated classes give OVL ~ 0",
      ovl_sep[1, 2] < 1e-4)

# -----------------------------------------------------------------------
# T7: OVL is symmetric in j,k  (already tested via matrix symmetry,
#     but also verify the raw .ovl_pair function directly)
# -----------------------------------------------------------------------
mu1 <- X %*% beta1; mu2 <- X %*% beta2
sd1 <- rep(sqrt(3.125), 5); sd2 <- sd1
ovl_12 <- .ovl_pair(mu1, mu2, sd1, sd2)
ovl_21 <- .ovl_pair(mu2, mu1, sd2, sd1)
check("T7  .ovl_pair is symmetric", abs(ovl_12 - ovl_21) < 1e-10)

# -----------------------------------------------------------------------
# T8: Multiple outcomes return a named list with correct length
# -----------------------------------------------------------------------
beta1_y2 <- 4 * c(0.214,  0.125,  0.010)
beta2_y2 <- 4 * c(0.000, -0.027,  0.002)
beta3_y2 <- 4 * c(1.688,  0.125, -0.010)

means_y2 <- list(Class_1 = X %*% beta1_y2,
                 Class_2 = X %*% beta2_y2,
                 Class_3 = X %*% beta3_y2)
cov_y2   <- 3.16 * diag(5)
covs_y2  <- list(cov_y2, cov_y2, cov_y2)

ovl_multi <- compute_ovl(
  means = list(Y1 = means_y1, Y2 = means_y2),
  covs  = list(Y1 = covs_y1,  Y2 = covs_y2),
  average_outcomes = TRUE
)
check("T8a Multiple outcomes returns a list",  is.list(ovl_multi))
check("T8b List has 3 elements (Y1, Y2, Average)", length(ovl_multi) == 3L)
check("T8c Y1 element is a matrix",  is.matrix(ovl_multi$Y1))
check("T8d Average element is a matrix", is.matrix(ovl_multi$Average))

# -----------------------------------------------------------------------
# T9: Average matrix is the mean of the per-outcome matrices
# -----------------------------------------------------------------------
expected_avg <- (ovl_multi$Y1 + ovl_multi$Y2) / 2
check("T9  Average matrix = mean(Y1, Y2)",
      max(abs(ovl_multi$Average - round(expected_avg, 3))) < 1e-3)

# -----------------------------------------------------------------------
# T10: Shared covariance recycled across classes
# -----------------------------------------------------------------------
# Pass a single matrix instead of a list of K matrices
ovl_recycled <- compute_ovl(means_y1, list(cov_y1))
check("T10 Single shared covariance recycled across classes",
      is.matrix(ovl_recycled) && all(diag(ovl_recycled) == 1))

# -----------------------------------------------------------------------
# T11: average_outcomes = FALSE omits the Average element
# -----------------------------------------------------------------------
ovl_no_avg <- compute_ovl(
  means = list(Y1 = means_y1, Y2 = means_y2),
  covs  = list(Y1 = covs_y1,  Y2 = covs_y2),
  average_outcomes = FALSE
)
check("T11 average_outcomes=FALSE: no Average element",
      is.list(ovl_no_avg) && is.null(ovl_no_avg$Average) && length(ovl_no_avg) == 2L)

# -----------------------------------------------------------------------
# T12: digits parameter controls rounding
# -----------------------------------------------------------------------
ovl_d2 <- compute_ovl(means_y1, covs_y1, digits = 2)
ovl_d4 <- compute_ovl(means_y1, covs_y1, digits = 4)
# Values rounded to 2 dp should differ from those rounded to 4 dp
check("T12 digits parameter affects rounding precision",
      !identical(ovl_d2, ovl_d4))
# More specific: all off-diagonal values have at most 2 decimal places
off_d2 <- ovl_d2[upper.tri(ovl_d2)]
check("T12b digits=2 gives values with at most 2 decimal places",
      all(off_d2 == round(off_d2, 2)))

# -----------------------------------------------------------------------
# T13: K=2 case works correctly
# -----------------------------------------------------------------------
means_k2 <- list(Class_1 = X %*% beta1, Class_2 = X %*% beta2)
covs_k2  <- list(cov_y1, cov_y1)
ovl_k2   <- compute_ovl(means_k2, covs_k2)
check("T13 K=2 returns 2x2 matrix", all(dim(ovl_k2) == c(2L, 2L)))
check("T13b K=2 diagonal = 1", all(diag(ovl_k2) == 1))

# -----------------------------------------------------------------------
# T14: OVL increases when classes move closer together
# -----------------------------------------------------------------------
# Move class 2 toward class 1: OVL should increase
means_closer <- list(Class_1 = X %*% beta1,
                     Class_2 = X %*% (beta1 * 0.5))  # halfway between
ovl_far    <- compute_ovl(list(Class_1 = X %*% beta1,
                               Class_2 = X %*% beta2), list(cov_y1, cov_y1))
ovl_close  <- compute_ovl(means_closer, list(cov_y1, cov_y1))
check("T14 OVL increases as classes move closer",
      ovl_close[1, 2] > ovl_far[1, 2])

cat("\n============================\n")
cat(sprintf("Results: %d passed, %d failed\n", pass, fail))
