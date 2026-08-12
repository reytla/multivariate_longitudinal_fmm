# =============================================================================
# ovl_trajectories.R
# -----------------------------------------------------------------------------
# Overlap Coefficient (OVL) for longitudinal finite mixture models.
#
# The OVL quantifies the degree of distributional overlap between pairs of
# latent classes at each measurement occasion. Values near 0 indicate
# well-separated classes; values near 1 indicate near-complete overlap.
#
# Method:
#   For each pair of classes (j, k) and each time point t, the OVL is the
#   integral of the minimum of two univariate normal densities:
#
#     OVL_t(j,k) = integral_{-Inf}^{Inf} min( N(x; mu_j(t), sigma_j(t)),
#                                              N(x; mu_k(t), sigma_k(t)) ) dx
#
#   The overall pairwise OVL is the mean of OVL_t across all T time points.
#   When multiple outcomes are provided, OVL matrices are returned per outcome
#   and optionally averaged across outcomes.
#
# Reference:
#   Nowakowska, E., Koronacki, J., & Lipovetsky, S. (2014). Tractable measure
#   of component overlap for Gaussian mixture models. arXiv:1407.7172 [math.ST].
#   https://arxiv.org/abs/1407.7172
#
# Usage:
#   See compute_ovl() and the worked example at the bottom of this file.
# =============================================================================


# -----------------------------------------------------------------------------
# Internal: OVL for a single pair of classes, single outcome
# -----------------------------------------------------------------------------
# mu1, mu2 : numeric vectors of length T (class means at each time point)
# sd1, sd2 : numeric vectors of length T (class SDs at each time point)
# Returns  : scalar OVL averaged across time points
.ovl_pair <- function(mu1, mu2, sd1, sd2) {
  stopifnot(length(mu1) == length(mu2),
            length(sd1) == length(sd2),
            length(mu1) == length(sd1))

  integrand <- function(x, m1, m2, s1, s2) {
    pmin(dnorm(x, mean = m1, sd = s1),
         dnorm(x, mean = m2, sd = s2))
  }

  T <- length(mu1)
  ovl_t <- numeric(T)
  for (t in seq_len(T)) {
    ovl_t[t] <- integrate(
      integrand, lower = -Inf, upper = Inf,
      m1 = mu1[t], m2 = mu2[t], s1 = sd1[t], s2 = sd2[t],
      subdivisions = 200L
    )$value
  }
  mean(ovl_t)
}


# -----------------------------------------------------------------------------
# Internal: full K x K OVL matrix for a single outcome
# -----------------------------------------------------------------------------
# means_k : list of K numeric vectors (one per class, length T)
# covs_k  : list of K numeric matrices (T x T covariance matrix per class);
#           only the diagonal (variances) is used
# Returns : K x K symmetric matrix with 1s on the diagonal
.ovl_matrix <- function(means_k, covs_k, digits) {
  K <- length(means_k)

  if (length(covs_k) != K) {
    stop("'means' and 'covs' must have the same number of classes.")
  }

  class_labels <- names(means_k)
  if (is.null(class_labels)) class_labels <- paste0("Class_", seq_len(K))

  mat <- matrix(1, nrow = K, ncol = K,
                dimnames = list(class_labels, class_labels))

  for (i in seq_len(K - 1)) {
    for (j in seq(i + 1, K)) {
      sd_i <- sqrt(diag(as.matrix(covs_k[[i]])))
      sd_j <- sqrt(diag(as.matrix(covs_k[[j]])))
      val  <- .ovl_pair(means_k[[i]], means_k[[j]], sd_i, sd_j)
      mat[i, j] <- mat[j, i] <- round(val, digits)
    }
  }
  mat
}


# =============================================================================
#' Compute pairwise Overlap Coefficients (OVL) for a longitudinal mixture model
#'
#' @description
#' Quantifies class separation in a finite mixture model with repeated measures
#' by computing the pairwise OVL between all latent classes. The OVL is the
#' area under the minimum of two probability density functions; it equals 0 for
#' perfectly separated classes and 1 for identical distributions.
#'
#' @param means
#'   \strong{Single outcome:} a named or unnamed list of K numeric vectors, one
#'   per class, each of length T (number of time points), containing the
#'   class-specific trajectory means.
#'
#'   \strong{Multiple outcomes:} a named or unnamed list of P such lists (one
#'   per outcome). Names are used as outcome labels in the output.
#'
#' @param covs
#'   \strong{Single outcome:} a list of K numeric matrices (T x T), one per
#'   class. Only the diagonal elements (variances) are used; off-diagonal
#'   elements are ignored. A single shared T x T matrix may be provided and
#'   will be recycled across all K classes.
#'
#'   \strong{Multiple outcomes:} a list of P such lists, matching the structure
#'   of \code{means}. A single shared list of K matrices will be recycled
#'   across all P outcomes.
#'
#' @param outcome_names Optional character vector of length P. Used to label
#'   output matrices. Defaults to names of \code{means} or
#'   \code{"Outcome_1", "Outcome_2", ...}.
#'
#' @param average_outcomes Logical (default \code{TRUE}). When multiple outcomes
#'   are provided, also return a matrix averaged across outcomes.
#'
#' @param digits Integer (default \code{3}). Decimal places to round OVL values.
#'
#' @return
#'   \strong{Single outcome:} a K x K symmetric numeric matrix. Diagonal = 1.
#'   Off-diagonal [j, k] = OVL between classes j and k, averaged across T time
#'   points.
#'
#'   \strong{Multiple outcomes:} a named list of K x K matrices (one per
#'   outcome), plus an \code{"Average"} matrix if
#'   \code{average_outcomes = TRUE}.
#'
#' @details
#'   The OVL at time t between classes j and k is:
#'   \deqn{OVL_t(j,k) = \int_{-\infty}^{\infty}
#'         \min\bigl(f_j(x,t),\, f_k(x,t)\bigr)\, dx}
#'   where \eqn{f_c(x,t) = \mathcal{N}(x;\, \mu_c(t),\, \sigma_c(t))}.
#'   The reported OVL is \eqn{\bar{OVL}(j,k) = T^{-1} \sum_t OVL_t(j,k)}.
#'
#'   Trajectory means \eqn{\mu_c(t)} are typically the model-implied fixed
#'   effects for class c at time t (e.g. \code{X \%*\% beta_c} for a polynomial
#'   time trend). Standard deviations are the square roots of the diagonal
#'   elements of the class-specific residual covariance matrix.
#'
#' @seealso
#'   \code{\link{print_ovl}} for formatted console output.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Example 1: single outcome, 3 classes, 5 time points
#' # -------------------------------------------------------------------------
#' # Quadratic time trend design matrix (t = 0, 1, 2, 3, 4)
#' t_vec  <- 0:4
#' X      <- cbind(1, t_vec, t_vec^2)
#'
#' # Class-specific fixed effect coefficients [intercept, slope, quadratic]
#' beta1  <- c(0,  1.752, -0.140)   # rising then falling
#' beta2  <- c(0,  0.000,  0.000)   # flat
#' beta3  <- c(0, -1.752,  0.140)   # falling then rising
#'
#' # Model-implied trajectory means per class (length-5 vectors)
#' means <- list(
#'   Class_1 = X %*% beta1,
#'   Class_2 = X %*% beta2,
#'   Class_3 = X %*% beta3
#' )
#'
#' # Shared diagonal residual covariance (variance = 3.125 at all time points)
#' cov_mat <- 3.125 * diag(5)
#' covs <- list(cov_mat, cov_mat, cov_mat)
#'
#' ovl <- compute_ovl(means, covs)
#' print_ovl(ovl)
#'
#' # -------------------------------------------------------------------------
#' # Example 2: two outcomes, 3 classes, 5 time points
#' # -------------------------------------------------------------------------
#' beta1_y2 <- c(0.856,  0.500,  0.040)
#' beta2_y2 <- c(0.000, -0.108,  0.008)
#' beta3_y2 <- c(6.752,  0.500, -0.040)
#'
#' means_y2 <- list(
#'   Class_1 = X %*% beta1_y2,
#'   Class_2 = X %*% beta2_y2,
#'   Class_3 = X %*% beta3_y2
#' )
#'
#' cov_mat2 <- 3.16 * diag(5)
#' covs_y2  <- list(cov_mat2, cov_mat2, cov_mat2)
#'
#' ovl_multi <- compute_ovl(
#'   means = list(Y1 = means, Y2 = means_y2),
#'   covs  = list(Y1 = covs,  Y2 = covs_y2),
#'   average_outcomes = TRUE
#' )
#' print_ovl(ovl_multi)
# =============================================================================
compute_ovl <- function(means,
                        covs,
                        outcome_names    = NULL,
                        average_outcomes = TRUE,
                        digits           = 3) {

  # --- Detect single vs multiple outcomes ------------------------------------
  # Single outcome: means[[1]] is a numeric vector
  # Multiple outcomes: means[[1]] is itself a list
  single_outcome <- is.numeric(means[[1]])

  if (single_outcome) {
    means_list <- list(means)
    covs_list  <- list(covs)
  } else {
    means_list <- means
    covs_list  <- covs
  }

  P <- length(means_list)

  # --- Recycle a single shared covariance list across all outcomes -----------
  if (!is.list(covs_list[[1]])) {
    # covs was supplied as a flat list of K matrices: share across all outcomes
    covs_list <- rep(list(covs_list), P)
  }

  # --- Outcome labels --------------------------------------------------------
  if (is.null(outcome_names)) {
    outcome_names <- names(means_list)
  }
  if (is.null(outcome_names) || length(outcome_names) != P) {
    outcome_names <- if (P == 1) "OVL" else paste0("Outcome_", seq_len(P))
  }
  names(means_list) <- names(covs_list) <- outcome_names

  # --- Recycle a single shared covariance matrix across all classes ----------
  for (p in seq_len(P)) {
    K <- length(means_list[[p]])
    if (length(covs_list[[p]]) == 1 && K > 1) {
      covs_list[[p]] <- rep(covs_list[[p]], K)
    }
  }

  # --- Compute OVL matrix per outcome ----------------------------------------
  results <- lapply(outcome_names, function(p) {
    .ovl_matrix(means_list[[p]], covs_list[[p]], digits)
  })
  names(results) <- outcome_names

  # --- Average across outcomes -----------------------------------------------
  if (P > 1 && average_outcomes) {
    avg_mat        <- Reduce("+", results) / P
    avg_mat        <- round(avg_mat, digits)
    diag(avg_mat)  <- 1
    results[["Average"]] <- avg_mat
  }

  # --- Return ----------------------------------------------------------------
  if (single_outcome) results[[1]] else results
}


# =============================================================================
#' Pretty-print OVL output
#'
#' @description
#' Prints OVL matrices from \code{\link{compute_ovl}} with a clear header,
#' highlighting that lower values indicate better class separation.
#'
#' @param ovl Output from \code{compute_ovl}: either a single matrix or a
#'   named list of matrices.
#' @param title Optional string printed above the output.
#'
#' @return Invisibly returns \code{ovl}.
#'
#' @examples
#' # See compute_ovl() examples.
# =============================================================================
print_ovl <- function(ovl, title = "Overlap Coefficient (OVL) Matrix") {
  cat(rep("=", nchar(title) + 4), "\n", sep = "")
  cat(" ", title, "\n")
  cat(rep("=", nchar(title) + 4), "\n", sep = "")
  cat(" Interpretation: 0 = perfect separation, 1 = complete overlap\n\n")

  if (is.matrix(ovl)) {
    print(ovl)
  } else {
    for (nm in names(ovl)) {
      cat("-- ", nm, " --\n", sep = "")
      print(ovl[[nm]])
      cat("\n")
    }
  }
  invisible(ovl)
}


# =============================================================================
# WORKED EXAMPLE
# (Remove or comment out before sourcing as a library file)
# =============================================================================
if (FALSE) {

  # --- Setup -----------------------------------------------------------------
  t_vec <- 0:4
  X     <- cbind(1, t_vec, t_vec^2)

  # --- Example 1: single outcome, 3 classes ----------------------------------
  beta1 <- 4 * c(0,  0.438, -0.035)
  beta2 <- 4 * c(0,  0.000,  0.000)
  beta3 <- 4 * c(0, -0.438,  0.035)

  means_y1 <- list(
    Class_1 = X %*% beta1,
    Class_2 = X %*% beta2,
    Class_3 = X %*% beta3
  )

  cov_y1 <- 3.125 * diag(5)           # shared diagonal covariance
  covs_y1 <- list(cov_y1, cov_y1, cov_y1)

  ovl_y1 <- compute_ovl(means_y1, covs_y1)
  print_ovl(ovl_y1, title = "OVL - Outcome Y1 (3 classes)")

  # --- Example 2: two outcomes, 3 classes ------------------------------------
  beta1_y2 <- 4 * c(0.214,  0.125,  0.010)
  beta2_y2 <- 4 * c(0.000, -0.027,  0.002)
  beta3_y2 <- 4 * c(1.688,  0.125, -0.010)

  means_y2 <- list(
    Class_1 = X %*% beta1_y2,
    Class_2 = X %*% beta2_y2,
    Class_3 = X %*% beta3_y2
  )

  cov_y2   <- 3.16 * diag(5)
  covs_y2  <- list(cov_y2, cov_y2, cov_y2)

  ovl_multi <- compute_ovl(
    means = list(Y1 = means_y1, Y2 = means_y2),
    covs  = list(Y1 = covs_y1,  Y2 = covs_y2),
    average_outcomes = TRUE
  )
  print_ovl(ovl_multi, title = "OVL - Two outcomes (3 classes)")

  # --- Example 3: class-specific (heterogeneous) covariance matrices ---------
  # Useful when classes have different residual variances
  cov_c1 <- diag(c(2.5, 2.8, 3.0, 3.2, 3.1))   # class 1: variance increases
  cov_c2 <- diag(c(3.0, 3.0, 3.0, 3.0, 3.0))   # class 2: constant variance
  cov_c3 <- diag(c(4.0, 3.5, 3.0, 2.5, 2.0))   # class 3: variance decreases

  ovl_het <- compute_ovl(
    means = means_y1,
    covs  = list(cov_c1, cov_c2, cov_c3)
  )
  print_ovl(ovl_het, title = "OVL - Heterogeneous covariance")
}
