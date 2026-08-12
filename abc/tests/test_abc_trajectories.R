# =============================================================================
# test_abc_trajectories.R
# -----------------------------------------------------------------------------
# Unit tests for abc/abc_trajectories.R using the testthat framework.
#
# Run via:
#   testthat::test_file("abc/tests/test_abc_trajectories.R")
#
# Or, from within the tests/ directory:
#   testthat::test_file("test_abc_trajectories.R")
# =============================================================================

library(testthat)

# Locate abc_trajectories.R relative to this test file.
# Works when called via test_file(), source(), or interactively.
.this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) getwd()
)
.module_path <- file.path(.this_dir, "..", "abc_trajectories.R")
if (!file.exists(.module_path)) {
  # Fallback: assume working directory is simulation_code root
  .module_path <- "abc/abc_trajectories.R"
}
source(.module_path, local = FALSE)


# =============================================================================
# Shared fixtures
# =============================================================================

# Population coefficient matrices from GBMTM_DGM.R  (K x 3, rows = classes)
pop_y1 <- matrix(4 * c( 0,  0.438, -0.035,
                         0,  0.000,  0.000,
                         0, -0.438,  0.035),
                 nrow = 3, ncol = 3, byrow = TRUE,
                 dimnames = list(paste0("k=", 1:3), c("B0","B1","B2")))

pop_y2 <- matrix(4 * c(0.214,  0.125,  0.010,
                        0.000, -0.027,  0.002,
                        1.688,  0.125, -0.010),
                 nrow = 3, ncol = 3, byrow = TRUE,
                 dimnames = list(paste0("k=", 1:3), c("B0","B1","B2")))


# =============================================================================
# Test 1: ABC = 0 when estimated equals population (analytical exactness)
# =============================================================================
test_that("ABC is exactly zero when estimated equals population parameters", {
  result <- compute_abc(pop_betas = pop_y1,
                        est_betas = pop_y1,
                        t_min = 0, t_max = 7)

  expect_equal(result$abc_mean, 0, tolerance = 1e-10)
  expect_true(all(result$abc_by_class == 0))
})


# =============================================================================
# Test 2: Known constant offset, analytical closed-form verification
# =============================================================================
# When d(t) = delta (constant), ABC = |delta| * (t_max - t_min).
# Here delta = 2, t_range = 7, so expected ABC = 14 per class.
# With K = 3 classes and all offsets equal, grand mean ABC = 14.
test_that("Constant-offset ABC matches closed-form solution delta * (t_max - t_min)", {
  delta  <- 2
  t_min  <- 0
  t_max  <- 7

  # All three classes shifted by delta on b0; identity permutation is optimal
  pop_const <- matrix(c(0, 0, 0,
                        0, 0, 0,
                        0, 0, 0),
                      nrow = 3, ncol = 3, byrow = TRUE)
  est_const <- pop_const
  est_const[, 1] <- est_const[, 1] + delta   # shift b0 for every class

  result <- compute_abc(pop_betas = pop_const,
                        est_betas = est_const,
                        t_min = t_min, t_max = t_max)

  expect_equal(result$abc_mean, delta * (t_max - t_min), tolerance = 1e-12)
  expect_true(all(abs(result$abc_by_class - delta * (t_max - t_min)) < 1e-10))
})


# =============================================================================
# Test 3: Correct permutation is recovered when classes 1 and 3 are swapped
# =============================================================================
test_that("Optimal permutation is correctly identified after a 1-3 class swap", {
  est_y1_swapped <- pop_y1[c(3, 2, 1), ]   # classes 1 and 3 exchanged

  result <- compute_abc(pop_betas = pop_y1,
                        est_betas = est_y1_swapped,
                        t_min = 0, t_max = 7)

  expect_equal(result$optimal_permutation, c(3L, 2L, 1L))
})


# =============================================================================
# Test 4: Optimal permutation has strictly lower ABC than identity when swapped
# =============================================================================
test_that("Optimal permutation ABC is lower than identity permutation ABC", {
  est_y1_swapped <- pop_y1[c(3, 2, 1), ]

  result <- compute_abc(pop_betas = pop_y1,
                        est_betas = est_y1_swapped,
                        t_min = 0, t_max = 7)

  identity_label <- "(1,2,3)"
  expect_true(result$abc_mean < result$all_permutations_abc[identity_label])
})


# =============================================================================
# Test 5: ABC is symmetric: swapping pop and est gives identical abc_mean
# =============================================================================
test_that("ABC is symmetric: compute_abc(A, B)$abc_mean == compute_abc(B, A)$abc_mean", {
  est_y1_swapped <- pop_y1[c(3, 2, 1), ]

  abc_forward  <- compute_abc(pop_betas = pop_y1,
                              est_betas = est_y1_swapped,
                              t_min = 0, t_max = 7)
  abc_backward <- compute_abc(pop_betas = est_y1_swapped,
                              est_betas = pop_y1,
                              t_min = 0, t_max = 7)

  expect_equal(abc_forward$abc_mean, abc_backward$abc_mean, tolerance = 1e-10)
})


# =============================================================================
# Test 6: Class relabeling is applied correctly to class_assignments
# =============================================================================
test_that("relabeled_assignments correctly maps estimated labels to population labels", {
  # Swap classes 1 and 3 in the estimated model
  est_y1_swapped <- pop_y1[c(3, 2, 1), ]

  # A known assignment vector: pop labels 1,2,3 map from est labels 3,2,1
  # So an estimated assignment of c(3,1,2,1,3) should relabel to c(1,3,2,3,1)
  est_assignments <- c(3L, 1L, 2L, 1L, 3L)
  expected_relabeled <- c(1L, 3L, 2L, 3L, 1L)

  result <- compute_abc(pop_betas         = pop_y1,
                        est_betas         = est_y1_swapped,
                        t_min             = 0,
                        t_max             = 7,
                        class_assignments = est_assignments)

  expect_equal(result$relabeled_assignments, expected_relabeled)
})


# =============================================================================
# Test 7: .all_permutations generates the correct number of unique permutations
# =============================================================================
test_that(".all_permutations generates K! unique permutations for K=3 and K=4", {
  perms3 <- .all_permutations(3L)
  expect_equal(length(perms3), 6L)

  # All are distinct
  perm3_strs <- vapply(perms3, paste, character(1), collapse = ",")
  expect_equal(length(unique(perm3_strs)), 6L)

  perms4 <- .all_permutations(4L)
  expect_equal(length(perms4), 24L)

  perm4_strs <- vapply(perms4, paste, character(1), collapse = ",")
  expect_equal(length(unique(perm4_strs)), 24L)

  # Each permutation is a complete set of 1:K
  for (p in perms3) expect_equal(sort(p), 1:3)
  for (p in perms4) expect_equal(sort(p), 1:4)
})


# =============================================================================
# Test 8: Two outcomes return correct structure and grand mean consistency
# =============================================================================
test_that("Two-outcome input returns correct structure and consistent grand mean", {
  result <- compute_abc(
    pop_betas = list(Y1 = pop_y1, Y2 = pop_y2),
    est_betas = list(Y1 = pop_y1[c(3,2,1), ], Y2 = pop_y2[c(3,2,1), ]),
    t_min = 0, t_max = 7
  )

  expect_equal(length(result$abc_by_outcome), 2L)
  expect_equal(length(result$abc_by_class), 3L)

  # Grand mean must equal mean of per-outcome means (within rounding)
  expect_equal(result$abc_mean,
               round(mean(result$abc_by_outcome), 4),
               tolerance = 1e-4)
  expect_equal(result$abc_mean,
               round(mean(result$abc_by_class), 4),
               tolerance = 1e-4)
})


# =============================================================================
# Test 9: Input validation catches malformed inputs
# =============================================================================
test_that("compute_abc stops with informative errors on malformed inputs", {

  # t_min >= t_max
  expect_error(
    compute_abc(pop_y1, pop_y1, t_min = 7, t_max = 0),
    regexp = "t_min.*t_max|t_max.*t_min",
    ignore.case = TRUE
  )

  # Mismatched number of classes between pop and est
  est_wrong_K <- pop_y1[1:2, ]   # only 2 rows instead of 3
  expect_error(
    compute_abc(pop_y1, est_wrong_K),
    regexp = "same number of rows|K"
  )

  # Wrong number of columns (not 3)
  est_wrong_cols <- cbind(pop_y1, extra = 1)   # 4 columns
  expect_error(
    compute_abc(pop_y1, est_wrong_cols),
    regexp = "3 columns"
  )

  # Mismatched outcome counts (list lengths differ)
  expect_error(
    compute_abc(
      pop_betas = list(Y1 = pop_y1, Y2 = pop_y2),
      est_betas = list(Y1 = pop_y1)
    ),
    regexp = "same number of outcomes"
  )
})


# =============================================================================
# Test 10: Linear difference (b2 = 0) is handled correctly
# =============================================================================
# When the quadratic term cancels, d(t) is linear with one crossing point.
# For d(t) = -t + 3.5 on [0,7]: root at t=3.5, ABC = 2 * integral_0^3.5 (3.5-t) dt
# = 2 * [3.5*3.5 - 3.5^2/2] = 2 * [12.25 - 6.125] = 2 * 6.125 = 12.25
test_that("Linear difference (b2 = 0) gives correct analytical ABC", {
  pop_linear <- matrix(c(0, 0, 0), nrow = 1, ncol = 3)   # K=1, flat at 0
  est_linear <- matrix(c(3.5, -1, 0), nrow = 1, ncol = 3) # d(t) = -3.5 + t -> wait

  # d(t) = pop - est = (0 - 3.5) + (0 - (-1))*t + 0*t^2 = -3.5 + t
  # root: -3.5 + t = 0 -> t = 3.5 (inside [0,7])
  # ABC = integral_0^3.5 |(-3.5+t)| dt + integral_3.5^7 |(-3.5+t)| dt
  #     = integral_0^3.5 (3.5-t) dt + integral_3.5^7 (t-3.5) dt
  #     = [3.5t - t^2/2]_0^3.5 + [t^2/2 - 3.5t]_3.5^7
  #     = (12.25 - 6.125) + ((24.5 - 24.5) - (6.125 - 12.25))
  #     = 6.125 + 6.125 = 12.25

  result <- compute_abc(pop_betas = pop_linear,
                        est_betas = est_linear,
                        t_min = 0, t_max = 7)

  expect_equal(result$abc_mean, 12.25, tolerance = 1e-12)
})
