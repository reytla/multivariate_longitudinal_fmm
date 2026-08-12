# =============================================================================
# abc_trajectories.R
# -----------------------------------------------------------------------------
# Area Between Curves (ABC) for label switching in longitudinal mixture models.
#
# In finite mixture models, class labels are unidentified: class "1" in one
# model fit may correspond to class "3" in another, or to the population's
# class "2". The ABC-based algorithm resolves this by trying all K! possible
# relabellings of the estimated classes, computing the average ABC between
# estimated and population trajectories for each, and selecting the relabelling
# that minimises this criterion.
#
# The ABC between two polynomial trajectories is computed analytically (not
# numerically). For trajectories of degree d, the difference d(t) is itself a
# polynomial of degree d whose absolute value integrates in closed form:
#
#   d(t)  = sum_{i=0}^{d} (b_i - eb_i) * t^i
#   F(t)  = sum_{i=1}^{d} (b_i - eb_i) * t^i / i        [antiderivative]
#   ABC   = sum over sub-intervals of |F(right) - F(left)|
#
# Sign-change points of d(t) within (t_min, t_max) are found via polyroot()
# (base R), which handles any polynomial degree. Coefficients may differ in
# length between population and estimated trajectories; the shorter vector is
# zero-padded to the degree of the longer before differencing.
# This gives machine-precision results: ABC = 0 exactly when estimated equals
# population parameters.
#
# Usage:
#   See compute_abc(), print_abc(), plot_abc(), and the worked example at the
#   bottom of this file.
# =============================================================================


# -----------------------------------------------------------------------------
# Internal: ABC for one class, one outcome
# -----------------------------------------------------------------------------
# beta_pop : numeric vector c(b0, b1, ..., bd): population trajectory coefficients
# beta_est : numeric vector of coefficients for the estimated trajectory;
#            may differ in length from beta_pop (shorter is zero-padded)
# t_min    : scalar integration lower bound
# t_max    : scalar integration upper bound
# Returns  : scalar ABC >= 0
.poly_abc_pair <- function(beta_pop, beta_est, t_min, t_max) {
  # Pad to equal length, then difference
  d  <- max(length(beta_pop), length(beta_est))
  dc <- c(beta_pop, rep(0, d - length(beta_pop))) -
        c(beta_est, rep(0, d - length(beta_est)))

  # Antiderivative: F(t) = sum_{i=1}^{d} dc[i] * t^i / i
  F <- function(t) {
    i <- seq_along(dc)
    sum(dc * t^i / i)
  }

  # Real roots of d(t) inside (t_min, t_max) via polyroot()
  # Trim trailing near-zero coefficients so polyroot receives the correct degree
  eps     <- .Machine$double.eps ^ 0.5
  nonzero <- which(abs(dc) > eps)

  internal_roots <- if (length(nonzero) > 0L && max(nonzero) >= 2L) {
    r <- polyroot(dc[seq_len(max(nonzero))])
    Re(r[abs(Im(r)) < eps & Re(r) > t_min + eps & Re(r) < t_max - eps])
  } else {
    numeric(0)
  }

  bp <- c(t_min, sort(internal_roots), t_max)
  sum(abs(diff(vapply(bp, F, numeric(1)))))
}


# -----------------------------------------------------------------------------
# Internal: generate all K! permutations of 1:K
# -----------------------------------------------------------------------------
# K       : positive integer
# Returns : list of K! integer vectors, each of length K
.all_permutations <- function(K) {
  if (K == 1L) return(list(1L))
  result <- vector("list", factorial(K))
  idx    <- 1L
  for (first in seq_len(K)) {
    rest <- setdiff(seq_len(K), first)
    for (sp in .all_permutations(length(rest))) {
      result[[idx]] <- c(first, rest[sp])
      idx           <- idx + 1L
    }
  }
  result
}


# -----------------------------------------------------------------------------
# Internal: mean ABC across outcomes and classes for one permutation
# -----------------------------------------------------------------------------
# perm            : integer vector length K; perm[k] is the estimated class
#                   that is matched to population class k
# pop_betas_list  : list of P matrices, each K x d (row = class, col = coeff)
# est_betas_list  : same structure as pop_betas_list; d may differ from pop
# t_min, t_max    : integration bounds
# Returns         : scalar mean ABC
.abc_for_permutation <- function(perm, pop_betas_list, est_betas_list,
                                 t_min, t_max) {
  P   <- length(pop_betas_list)
  K   <- nrow(pop_betas_list[[1]])
  acc <- 0
  for (p in seq_len(P)) {
    for (k in seq_len(K)) {
      acc <- acc + .poly_abc_pair(
        beta_pop = pop_betas_list[[p]][k, ],
        beta_est = est_betas_list[[p]][perm[k], ],
        t_min    = t_min,
        t_max    = t_max
      )
    }
  }
  acc / (P * K)
}


# -----------------------------------------------------------------------------
# Internal: apply a class-label permutation to a class membership vector
# -----------------------------------------------------------------------------
# class_vec : integer vector of estimated class labels (values in 1:K)
# perm      : integer vector length K; perm[k] = estimated label -> pop class k
# Returns   : integer vector of same length with relabelled classes
.apply_permutation <- function(class_vec, perm) {
  # Build inverse map: inv_perm[perm[k]] = k
  inv_perm          <- integer(length(perm))
  inv_perm[perm]    <- seq_along(perm)
  inv_perm[class_vec]
}


# =============================================================================
#' Compute ABC-based optimal class-label alignment for a mixture model
#'
#' @description
#' Solves the label-switching problem for longitudinal mixture models with
#' polynomial trajectory means of any degree. For each of the K! possible
#' relabellings of the estimated classes, the algorithm computes the average
#' Area Between Curves (ABC) between estimated and population trajectories
#' across all classes and outcomes. The relabelling that minimises this
#' criterion is returned as the optimal permutation.
#'
#' ABC values are computed analytically (not via numerical quadrature): for
#' trajectories of degree d, the difference polynomial has degree d and its
#' antiderivative is \eqn{F(t) = \sum_{i=1}^{d} \Delta b_i \, t^i / i}.
#' Sign-change points within \code{[t_min, t_max]} are found via
#' \code{polyroot()} (base R), yielding machine-precision results (ABC = 0
#' exactly when estimated equals population). Population and estimated
#' trajectories may have different numbers of coefficients; the shorter
#' is zero-padded to the degree of the longer.
#'
#' @param pop_betas
#'   \strong{Single outcome:} a K x d numeric matrix, where row k contains the
#'   population fixed-effect coefficients \code{c(b0, b1, ..., b_{d-1})} for
#'   class k. Any polynomial degree d >= 1 is supported.
#'
#'   \strong{Multiple outcomes:} a named or unnamed list of P such matrices (one
#'   per outcome). Matrices across outcomes may have different numbers of
#'   columns (degrees).
#'
#'   \emph{Note on DGM convention:} the data-generating scripts store
#'   coefficients as a d x K matrix \code{bmeanmat} (rows = coefficients,
#'   columns = classes). Pass \code{t(bmeanmat)} to obtain the K x d matrix
#'   expected here.
#'
#' @param est_betas Estimated trajectory coefficients. Must have the same
#'   structure as \code{pop_betas}: K x d matrix for a single outcome, or a
#'   list of P such matrices for multiple outcomes. The number of columns
#'   (degree) may differ from \code{pop_betas}; the shorter is zero-padded.
#'
#' @param t_min Scalar lower bound of the integration domain. Default \code{0}.
#'
#' @param t_max Scalar upper bound of the integration domain. Default \code{7}
#'   (matching the DGM time range \code{seq(0, 7, length = 5)}).
#'
#' @param class_assignments Optional integer vector of length N containing
#'   estimated class memberships (values in \code{1:K}). When supplied, the
#'   function also returns a relabelled version using the optimal permutation.
#'   Default \code{NULL}.
#'
#' @param outcome_names Optional character vector of length P used to label
#'   output. Defaults to names of \code{pop_betas} (if a named list), or
#'   \code{"Outcome_1"}, \code{"Outcome_2"}, etc.
#'
#' @param digits Integer. Decimal places to round reported ABC values.
#'   Default \code{4}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{abc_by_class}}{Length-K numeric vector: ABC per class under the
#'     optimal permutation, averaged across outcomes.}
#'   \item{\code{abc_by_outcome}}{Length-P numeric vector: ABC per outcome under
#'     the optimal permutation, averaged across classes.}
#'   \item{\code{abc_mean}}{Scalar grand mean ABC (equals
#'     \code{mean(abc_by_class)} and \code{mean(abc_by_outcome)}).}
#'   \item{\code{optimal_permutation}}{Integer vector of length K.
#'     \code{optimal_permutation[k]} is the estimated class label that best
#'     corresponds to population class k.}
#'   \item{\code{all_permutations_abc}}{Named numeric vector of length K!:
#'     mean ABC for every permutation tried. Useful for diagnostics (e.g. to
#'     check whether the minimum is well-separated from alternatives).}
#'   \item{\code{relabeled_assignments}}{Integer vector (same length as
#'     \code{class_assignments}) with optimal relabelling applied, or
#'     \code{NULL} if \code{class_assignments} was not supplied.}
#'   \item{\code{t_min}, \code{t_max}}{Integration bounds, echoed back.}
#' }
#'
#' @details
#'   The integration domain \code{[t_min, t_max]} should span the observed
#'   time range. For the DGM scripts in this repository the time values are
#'   \code{seq(0, 7, length = 5)} = {0, 1.75, 3.5, 5.25, 7}, so the default
#'   \code{t_min = 0}, \code{t_max = 7} is appropriate.
#'
#'   For K = 3 (the default simulation design), all 3! = 6 permutations are
#'   evaluated. Computation time is negligible even inside simulation loops
#'   thanks to the analytical integration.
#'
#' @seealso \code{\link{print_abc}} for formatted console output,
#'   \code{\link{plot_abc}} for visualisation.
#'
#' @examples
#' # -------------------------------------------------------------------------
#' # Example 1: single outcome, perfect recovery (ABC should be 0)
#' # -------------------------------------------------------------------------
#' pop_betas <- matrix(4 * c( 0,  0.438, -0.035,
#'                             0,  0.000,  0.000,
#'                             0, -0.438,  0.035),
#'                     nrow = 3, ncol = 3, byrow = TRUE,
#'                     dimnames = list(paste0("k=", 1:3), c("B0","B1","B2")))
#'
#' abc_zero <- compute_abc(pop_betas, pop_betas, t_min = 0, t_max = 7)
#' print_abc(abc_zero)   # abc_mean = 0
#'
#' # -------------------------------------------------------------------------
#' # Example 2: two outcomes, classes 1 and 3 swapped in the estimated model
#' # -------------------------------------------------------------------------
#' pop_y2 <- matrix(4 * c(0.214,  0.125,  0.010,
#'                         0.000, -0.027,  0.002,
#'                         1.688,  0.125, -0.010),
#'                  nrow = 3, ncol = 3, byrow = TRUE)
#' est_y1 <- pop_betas[c(3,2,1), ]   # swap classes 1 and 3
#' est_y2 <- pop_y2[c(3,2,1), ]
#'
#' abc_res <- compute_abc(
#'   pop_betas        = list(Y1 = pop_betas, Y2 = pop_y2),
#'   est_betas        = list(Y1 = est_y1,    Y2 = est_y2),
#'   class_assignments = c(1L, 2L, 3L, 1L, 3L)
#' )
#' print_abc(abc_res)   # optimal_permutation should be c(3L, 2L, 1L)
# =============================================================================
compute_abc <- function(pop_betas,
                        est_betas,
                        t_min             = 0,
                        t_max             = 7,
                        class_assignments = NULL,
                        outcome_names     = NULL,
                        digits            = 4) {

  # --- Validate integration bounds -------------------------------------------
  if (!is.numeric(t_min) || !is.numeric(t_max) || length(t_min) != 1 ||
      length(t_max) != 1 || t_min >= t_max) {
    stop("'t_min' and 't_max' must be scalars with t_min < t_max.")
  }

  # --- Detect single vs multiple outcomes ------------------------------------
  # Single outcome: pop_betas is a matrix; multiple: it is a list of matrices
  single_outcome <- is.matrix(pop_betas)

  if (single_outcome) {
    pop_list <- list(pop_betas)
    est_list <- list(est_betas)
  } else {
    if (!is.list(pop_betas) || !is.list(est_betas)) {
      stop("'pop_betas' and 'est_betas' must both be matrices (single outcome) ",
           "or lists of matrices (multiple outcomes).")
    }
    pop_list <- pop_betas
    est_list <- est_betas
  }

  P <- length(pop_list)

  # --- Validate dimensions ---------------------------------------------------
  if (length(est_list) != P) {
    stop("'pop_betas' and 'est_betas' must have the same number of outcomes.")
  }
  for (p in seq_len(P)) {
    if (!is.matrix(pop_list[[p]]) || !is.matrix(est_list[[p]])) {
      stop("Each element of 'pop_betas' and 'est_betas' must be a matrix.")
    }
    if (ncol(pop_list[[p]]) < 1L || ncol(est_list[[p]]) < 1L) {
      stop("Each beta matrix must have at least 1 column.")
    }
    if (nrow(pop_list[[p]]) != nrow(est_list[[p]])) {
      stop("'pop_betas' and 'est_betas' must have the same number of rows (K).")
    }
  }

  K <- nrow(pop_list[[1]])
  if (K < 1) stop("'K' (number of classes) must be at least 1.")

  # --- Outcome labels --------------------------------------------------------
  if (is.null(outcome_names)) outcome_names <- names(pop_list)
  if (is.null(outcome_names) || length(outcome_names) != P) {
    outcome_names <- if (P == 1) "Outcome_1" else paste0("Outcome_", seq_len(P))
  }
  names(pop_list) <- names(est_list) <- outcome_names

  # --- Generate all K! permutations ------------------------------------------
  perms     <- .all_permutations(K)
  n_perms   <- length(perms)

  perm_abc  <- vapply(perms,
                      FUN       = .abc_for_permutation,
                      FUN.VALUE = numeric(1),
                      pop_betas_list = pop_list,
                      est_betas_list = est_list,
                      t_min          = t_min,
                      t_max          = t_max)

  # --- Name each permutation for diagnostics ---------------------------------
  perm_labels      <- vapply(perms,
                             function(p) paste0("(", paste(p, collapse = ","), ")"),
                             character(1))
  names(perm_abc)  <- perm_labels

  # --- Select optimal permutation --------------------------------------------
  opt_idx  <- which.min(perm_abc)
  opt_perm <- perms[[opt_idx]]

  # --- Per-class and per-outcome breakdown under optimal permutation ----------
  abc_matrix <- matrix(0, nrow = K, ncol = P,
                       dimnames = list(paste0("Class_", seq_len(K)),
                                       outcome_names))
  for (p in seq_len(P)) {
    for (k in seq_len(K)) {
      abc_matrix[k, p] <- .poly_abc_pair(
        beta_pop = pop_list[[p]][k, ],
        beta_est = est_list[[p]][opt_perm[k], ],
        t_min    = t_min,
        t_max    = t_max
      )
    }
  }

  abc_by_class   <- round(rowMeans(abc_matrix), digits)
  abc_by_outcome <- round(colMeans(abc_matrix), digits)
  abc_mean       <- round(mean(abc_matrix), digits)

  # --- Relabel class assignments if provided ---------------------------------
  relabeled <- NULL
  if (!is.null(class_assignments)) {
    if (!is.integer(class_assignments)) {
      class_assignments <- as.integer(class_assignments)
    }
    if (any(class_assignments < 1L | class_assignments > K)) {
      stop("'class_assignments' contains values outside 1:K.")
    }
    relabeled <- .apply_permutation(class_assignments, opt_perm)
  }

  list(
    abc_by_class          = abc_by_class,
    abc_by_outcome        = abc_by_outcome,
    abc_mean              = abc_mean,
    optimal_permutation   = opt_perm,
    all_permutations_abc  = round(perm_abc, digits),
    relabeled_assignments = relabeled,
    t_min                 = t_min,
    t_max                 = t_max
  )
}


# =============================================================================
#' Pretty-print ABC output
#'
#' @description
#' Prints a formatted summary of \code{\link{compute_abc}} results, showing the
#' optimal class-label permutation, grand mean ABC, and per-class and
#' per-outcome breakdowns.
#'
#' @param abc Output from \code{compute_abc}.
#' @param title Optional string printed as a header.
#'
#' @return Invisibly returns \code{abc}.
#'
#' @examples
#' # See compute_abc() examples.
# =============================================================================
print_abc <- function(abc, title = "Area Between Curves (ABC): Label Switching") {
  cat(rep("=", nchar(title) + 4), "\n", sep = "")
  cat(" ", title, "\n")
  cat(rep("=", nchar(title) + 4), "\n", sep = "")
  cat(" Integration domain: [", abc$t_min, ", ", abc$t_max, "]\n", sep = "")

  K <- length(abc$optimal_permutation)
  perm_str <- paste(sprintf("est.%d -> pop.%d", abc$optimal_permutation,
                            seq_len(K)),
                    collapse = ";  ")
  cat(" Optimal permutation: ", perm_str, "\n", sep = "")
  cat(" Grand mean ABC:      ", abc$abc_mean, "\n\n", sep = "")

  cat(" ABC by class (averaged across outcomes):\n")
  print(abc$abc_by_class)

  if (length(abc$abc_by_outcome) > 1) {
    cat("\n ABC by outcome (averaged across classes):\n")
    print(abc$abc_by_outcome)
  }

  cat("\n All permutations (mean ABC):\n")
  print(abc$all_permutations_abc)

  invisible(abc)
}


# =============================================================================
#' Plot estimated and population trajectories with shaded ABC regions
#'
#' @description
#' Produces a panel grid (P rows x K columns) of trajectory plots. Each panel
#' shows the population trajectory (solid black) and the optimally-relabelled
#' estimated trajectory (dashed red) for one class and one outcome. The area
#' between the two curves is shaded to visualise the ABC.
#'
#' Uses base R graphics only. No external package dependencies.
#'
#' @param abc_result Output from \code{\link{compute_abc}}.
#'
#' @param pop_betas Population trajectory coefficients (same format as passed
#'   to \code{compute_abc}): K x 3 matrix or list of P such matrices.
#'
#' @param est_betas Estimated trajectory coefficients (same format).
#'
#' @param t_min Scalar lower bound of the plot/integration domain.
#'   Default \code{0}.
#'
#' @param t_max Scalar upper bound. Default \code{7}.
#'
#' @param outcome_names Optional character vector of length P for panel row
#'   labels. Defaults to names from \code{abc_result$abc_by_outcome}.
#'
#' @param class_names Optional character vector of length K for panel column
#'   labels. Defaults to \code{"Class 1"}, \code{"Class 2"}, etc.
#'
#' @param n_points Integer. Number of points used to draw each curve.
#'   Default \code{200}.
#'
#' @param col_pop Colour for the population trajectory line. Default
#'   \code{"black"}.
#'
#' @param col_est Colour for the estimated trajectory line. Default
#'   \code{"#E41A1C"} (red).
#'
#' @param shade_col Fill colour for the shaded ABC region. Default: a
#'   semi-transparent red.
#'
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' # See compute_abc() examples.
# =============================================================================
plot_abc <- function(abc_result,
                     pop_betas,
                     est_betas,
                     t_min         = 0,
                     t_max         = 7,
                     outcome_names = NULL,
                     class_names   = NULL,
                     n_points      = 200,
                     col_pop       = "black",
                     col_est       = "#E41A1C",
                     shade_col     = adjustcolor("#E41A1C", alpha.f = 0.15),
                     ...) {

  # --- Normalise inputs to list-of-P-matrices --------------------------------
  single_outcome <- is.matrix(pop_betas)
  if (single_outcome) {
    pop_list <- list(pop_betas)
    est_list <- list(est_betas)
  } else {
    pop_list <- pop_betas
    est_list <- est_betas
  }

  P    <- length(pop_list)
  K    <- nrow(pop_list[[1]])
  perm <- abc_result$optimal_permutation

  # --- Labels ----------------------------------------------------------------
  if (is.null(outcome_names)) {
    outcome_names <- names(abc_result$abc_by_outcome)
    if (is.null(outcome_names)) outcome_names <- paste0("Outcome ", seq_len(P))
  }
  if (is.null(class_names)) class_names <- paste0("Class ", seq_len(K))

  # --- Plot layout -----------------------------------------------------------
  old_par <- par(mfrow = c(P, K), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)

  t_seq <- seq(t_min, t_max, length.out = n_points)

  for (p in seq_len(P)) {
    for (k in seq_len(K)) {
      beta_pop <- pop_list[[p]][k, ]
      beta_est <- est_list[[p]][perm[k], ]

      y_pop <- drop(outer(t_seq, seq_along(beta_pop) - 1L, `^`) %*% beta_pop)
      y_est <- drop(outer(t_seq, seq_along(beta_est) - 1L, `^`) %*% beta_est)

      y_range <- range(c(y_pop, y_est))
      y_pad   <- diff(y_range) * 0.1
      if (y_pad == 0) y_pad <- 1

      plot(t_seq, y_pop,
           type  = "n",
           ylim  = c(y_range[1] - y_pad, y_range[2] + y_pad),
           xlab  = "Time",
           ylab  = "Mean trajectory",
           main  = paste0(outcome_names[p], ": ", class_names[k]),
           ...)

      # Shaded area between curves
      polygon(c(t_seq, rev(t_seq)),
              c(y_pop, rev(y_est)),
              col    = shade_col,
              border = NA)

      lines(t_seq, y_pop, col = col_pop, lwd = 2, lty = 1)
      lines(t_seq, y_est, col = col_est, lwd = 2, lty = 2)

      mtext(sprintf("ABC = %.4f", abc_result$abc_by_class[k]),
            side = 3, line = 0.1, cex = 0.8, col = "grey40")
    }
  }

  # --- Legend ----------------------------------------------------------------
  # Draw on the last panel using a small inset legend
  legend("topright",
         legend = c("Population", "Estimated (relabelled)"),
         col    = c(col_pop, col_est),
         lty    = c(1, 2),
         lwd    = 2,
         bty    = "n",
         cex    = 0.85)

  invisible(NULL)
}


# =============================================================================
# WORKED EXAMPLE
# (Wrapped in if (FALSE) so it does not run when this file is sourced.
#  Remove the if (FALSE) { ... } guard to execute interactively.)
# =============================================================================
if (FALSE) {

  # --- Coefficient matrices from GBMTM_DGM.R (transposed to K x 3) ----------
  # bmeanmat is 3 x K in the DGMs; t(bmeanmat) gives K x 3 as required here.
  bmeanmat <- matrix(4 * c( 0,  0.438, -0.035,
                             0,  0.000,  0.000,
                             0, -0.438,  0.035),
                     nrow = 3, ncol = 3,
                     dimnames = list(c("B0","B1","B2"), paste0("k=",1:3)))

  bmeanmat2 <- matrix(4 * c(0.214,  0.125,  0.010,
                             0.000, -0.027,  0.002,
                             1.688,  0.125, -0.010),
                      nrow = 3, ncol = 3,
                      dimnames = list(c("B0","B1","B2"), paste0("k=",1:3)))

  pop_y1 <- t(bmeanmat)    # 3 x 3  ->  K x 3
  pop_y2 <- t(bmeanmat2)

  # ---------------------------------------------------------------------------
  # Example 1: single outcome, estimated = population (ABC must be exactly 0)
  # ---------------------------------------------------------------------------
  abc_zero <- compute_abc(pop_betas = pop_y1,
                          est_betas = pop_y1,
                          t_min = 0, t_max = 7)
  print_abc(abc_zero, title = "Example 1: Perfect recovery (ABC = 0)")
  # Expected: abc_mean = 0, optimal_permutation = c(1, 2, 3) or any identity

  # ---------------------------------------------------------------------------
  # Example 2: two outcomes, classes 1 and 3 swapped in the estimated model
  # ---------------------------------------------------------------------------
  set.seed(42)
  est_y1 <- pop_y1[c(3, 2, 1), ]    # swap classes 1 and 3
  est_y2 <- pop_y2[c(3, 2, 1), ]

  abc_res <- compute_abc(
    pop_betas         = list(Y1 = pop_y1, Y2 = pop_y2),
    est_betas         = list(Y1 = est_y1, Y2 = est_y2),
    t_min             = 0,
    t_max             = 7,
    class_assignments = sample(1L:3L, 30L, replace = TRUE)
  )
  print_abc(abc_res, title = "Example 2: Labels 1 & 3 swapped")
  # Expected: optimal_permutation = c(3, 2, 1)

  # Visualise
  plot_abc(abc_res,
           pop_betas = list(Y1 = pop_y1, Y2 = pop_y2),
           est_betas = list(Y1 = est_y1, Y2 = est_y2),
           t_min = 0, t_max = 7)

  # ---------------------------------------------------------------------------
  # Example 3: Simulation over 100 replications
  #
  # Notation:
  #   true_perm[k] = estimated class label that corresponds to population class k
  #                  (i.e. true_perm IS the optimal permutation we expect to recover)
  #
  # Each replication:
  #   1. Draws true_perm at random (the "true" optimal permutation).
  #   2. Constructs estimated betas as pop[ inv(true_perm), ] + noise, so that
  #      estimated row j ≈ population row true_perm^{-1}[j], and therefore
  #      compute_abc() should return optimal_permutation == true_perm.
  #   3. Scrambles N = 90 class assignments: subjects in pop class k are labelled
  #      as estimated class true_perm[k].
  #   4. Runs compute_abc() to recover the optimal permutation and relabelled
  #      assignments.
  #   5. Records whether the correct permutation was identified and whether
  #      the relabelled assignments match the true population classes.
  # ---------------------------------------------------------------------------
  set.seed(2024)

  n_rep    <- 100L
  N        <- 90L          # subjects per replication
  K        <- 3L
  noise_sd <- 0.05         # std dev of per-coefficient Gaussian noise

  # Helper: invert a permutation (inv_perm[ perm[k] ] = k)
  invert_perm <- function(p) { inv <- integer(length(p)); inv[p] <- seq_along(p); inv }

  # Storage
  correct_perm  <- logical(n_rep)   # was the optimal permutation correct?
  abc_means     <- numeric(n_rep)   # grand mean ABC per replication
  relabels_list <- vector("list", n_rep)  # relabelled assignments (length N each)
  saved_perms   <- vector("list", n_rep)  # true_perm per replication
  saved_est_y1  <- vector("list", n_rep)  # estimated Y1 betas per replication
  saved_est_y2  <- vector("list", n_rep)  # estimated Y2 betas per replication
  saved_abc     <- vector("list", n_rep)  # full compute_abc output per replication

  # True class assignments (balanced: 30 per class, fixed across reps)
  true_classes <- rep(1L:K, each = N %/% K)

  for (rep in seq_len(n_rep)) {

    # 1. Draw the "true" optimal permutation: true_perm[k] = est class for pop class k
    true_perm <- sample.int(K)

    # 2. Build estimated betas using inv(true_perm) as the row index, so that
    #    est[ inv(true_perm)[k], ] ≈ pop[k, ]  =>  optimal_permutation = true_perm
    row_idx    <- invert_perm(true_perm)
    est_y1_sim <- pop_y1[row_idx, ] + matrix(rnorm(K * 3, sd = noise_sd), K, 3)
    est_y2_sim <- pop_y2[row_idx, ] + matrix(rnorm(K * 3, sd = noise_sd), K, 3)

    # 3. Scramble class assignments: pop class k -> est class true_perm[k]
    scrambled_classes <- true_perm[true_classes]

    # 4. ABC-based label alignment
    abc_rep <- compute_abc(
      pop_betas         = list(Y1 = pop_y1, Y2 = pop_y2),
      est_betas         = list(Y1 = est_y1_sim, Y2 = est_y2_sim),
      t_min             = 0,
      t_max             = 7,
      class_assignments = scrambled_classes
    )

    # 5. Record results
    correct_perm[rep]    <- identical(abc_rep$optimal_permutation, true_perm)
    abc_means[rep]       <- abc_rep$abc_mean
    relabels_list[[rep]] <- abc_rep$relabeled_assignments
    saved_perms[[rep]]   <- true_perm
    saved_est_y1[[rep]]  <- est_y1_sim
    saved_est_y2[[rep]]  <- est_y2_sim
    saved_abc[[rep]]     <- abc_rep
  }

  # --- Summary ----------------------------------------------------------------
  cat("\n===== Simulation Summary (", n_rep, " replications) =====\n", sep = "")
  cat(sprintf(" Correct permutation recovered : %d / %d  (%.1f%%)\n",
              sum(correct_perm), n_rep, 100 * mean(correct_perm)))
  cat(sprintf(" Grand mean ABC  mean (SD)      : %.4f (%.4f)\n",
              mean(abc_means), sd(abc_means)))
  cat(sprintf("                range           : [%.4f, %.4f]\n",
              min(abc_means), max(abc_means)))

  # Check relabelling accuracy: after relabelling, do recovered labels == truth?
  relabel_accuracy <- vapply(relabels_list,
                             function(rl) mean(rl == true_classes),
                             numeric(1))
  cat(sprintf(" Mean relabelling accuracy      : %.4f (SD %.4f)\n",
              mean(relabel_accuracy), sd(relabel_accuracy)))

  # Quick histogram of ABC values
  hist(abc_means,
       breaks = 20,
       main   = "Distribution of grand mean ABC across 100 replications",
       xlab   = "Grand mean ABC",
       col    = "steelblue",
       border = "white")
  abline(v = mean(abc_means), col = "firebrick", lwd = 2, lty = 2)
  legend("topright", legend = sprintf("Mean = %.4f", mean(abc_means)),
         col = "firebrick", lty = 2, lwd = 2, bty = "n")

  # --- Plot A: single-replication trajectory plot ----------------------------
  # Use plot_abc() for the first replication that has a non-identity permutation,
  # so that the relabelling does real work (classes are actually swapped).
  nontrivial_rep <- which(
    vapply(saved_perms, function(p) !identical(p, seq_len(K)), logical(1))
  )[1]
  cat(sprintf("\nPlot A: replication %d, injected permutation (%s)\n",
              nontrivial_rep,
              paste(saved_perms[[nontrivial_rep]], collapse = ",")))
  plot_abc(
    abc_result    = saved_abc[[nontrivial_rep]],
    pop_betas     = list(Y1 = pop_y1, Y2 = pop_y2),
    est_betas     = list(Y1 = saved_est_y1[[nontrivial_rep]],
                         Y2 = saved_est_y2[[nontrivial_rep]]),
    t_min         = 0,
    t_max         = 7
  )

  # --- Plot B: fan plot, population vs all relabelled estimated trajectories ---
  # One panel per class, outcome Y1. Thin blue lines = estimated trajectory from
  # each replication after ABC relabelling; thick black = population.
  t_seq      <- seq(0, 7, length.out = 200)
  eval_quad  <- function(beta, t) beta[1] + beta[2] * t + beta[3] * t^2

  old_par <- par(mfrow = c(1, K), mar = c(4, 4, 3, 1))
  for (k in seq_len(K)) {

    # Estimated trajectory for class k in each replication, after relabelling:
    # optimal_permutation[k] gives the estimated row that maps to population class k
    est_curves <- lapply(seq_len(n_rep), function(r) {
      est_row <- saved_abc[[r]]$optimal_permutation[k]
      eval_quad(saved_est_y1[[r]][est_row, ], t_seq)
    })
    pop_curve <- eval_quad(pop_y1[k, ], t_seq)

    ylim_k <- range(c(unlist(est_curves), pop_curve))
    plot(t_seq, pop_curve,
         type  = "n",
         ylim  = ylim_k,
         xlab  = "Time",
         ylab  = "Mean trajectory",
         main  = sprintf("Class %d (Y1), %d reps", k, n_rep))

    for (r in seq_len(n_rep))
      lines(t_seq, est_curves[[r]],
            col = adjustcolor("steelblue", alpha.f = 0.25), lwd = 1)

    lines(t_seq, pop_curve, col = "black", lwd = 2.5)
  }
  # Shared legend on the last panel
  legend("topright",
         legend = c("Population", "Estimated (relabelled)"),
         col    = c("black", adjustcolor("steelblue", alpha.f = 0.7)),
         lwd    = c(2.5, 1),
         bty    = "n",
         cex    = 0.85)
  par(old_par)

  # ---------------------------------------------------------------------------
  # Example 4: Mixed polynomial degrees (linear, quadratic, cubic)
  #
  # Population trajectories on [0, 7]:
  #   Class 1 (linear):    f(t) = 2 + 0.4*t                  [2.0 -> 4.8]
  #   Class 2 (quadratic): f(t) = 4 + 0.5*t - 0.1*t^2       [4.0 -> 2.6]
  #   Class 3 (cubic):     f(t) = 6 - 0.3*t^2 + 0.03*t^3    [6.0 -> 2.7]
  #
  # All stored in a K x 4 matrix, zero-padded to degree 3 (the highest present).
  # Estimated model has classes 1 and 3 swapped.
  # Expected: optimal_permutation = c(3, 2, 1), abc_mean = 0 (exact recovery).
  # ---------------------------------------------------------------------------
  pop_mixed <- matrix(c(2,  0.4,  0,    0,
                        4,  0.5, -0.1,  0,
                        6,  0,   -0.3,  0.03),
                      nrow = 3, ncol = 4, byrow = TRUE,
                      dimnames = list(paste0("k=", 1:3),
                                      c("B0", "B1", "B2", "B3")))

  # Swap classes 1 and 3 in the estimated model
  est_mixed <- pop_mixed[c(3, 2, 1), ]

  set.seed(1)
  abc_mixed <- compute_abc(
    pop_betas         = pop_mixed,
    est_betas         = est_mixed,
    t_min             = 0,
    t_max             = 7,
    class_assignments = sample(1L:3L, 30L, replace = TRUE)
  )
  print_abc(abc_mixed, title = "Example 4: Mixed degrees (linear, quadratic, cubic)")
  # Expected: optimal_permutation = c(3, 2, 1), abc_mean = 0

  # Trajectory plot: panel grid (1 outcome x 3 classes)
  # Solid black = population, dashed red = estimated after optimal relabelling
  plot_abc(abc_mixed,
           pop_betas = pop_mixed,
           est_betas = est_mixed,
           t_min     = 0,
           t_max     = 7)

  # Overlay all three population shapes on one panel for comparison
  t_seq    <- seq(0, 7, length.out = 200)
  eval_poly <- function(beta, t) drop(outer(t, seq_along(beta) - 1L, `^`) %*% beta)

  cols   <- c("steelblue", "firebrick", "forestgreen")
  ltypes <- c(1, 2, 4)
  labels <- c("Class 1 (linear)", "Class 2 (quadratic)", "Class 3 (cubic)")

  ylim_all <- range(vapply(seq_len(3), function(k)
    range(eval_poly(pop_mixed[k, ], t_seq)), numeric(2)))

  plot(t_seq, eval_poly(pop_mixed[1, ], t_seq),
       type = "l", col = cols[1], lty = ltypes[1], lwd = 2,
       ylim = ylim_all, xlab = "Time", ylab = "f(t)",
       main = "Population trajectories: mixed degrees")
  for (k in 2:3)
    lines(t_seq, eval_poly(pop_mixed[k, ], t_seq),
          col = cols[k], lty = ltypes[k], lwd = 2)
  legend("topright", legend = labels, col = cols, lty = ltypes,
         lwd = 2, bty = "n", cex = 0.85)
}
