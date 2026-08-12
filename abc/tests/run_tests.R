# Set working directory to project root when run via Rscript
.file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.file_arg))
  setwd(normalizePath(file.path(dirname(sub("--file=", "", .file_arg[1L])), "..", ".."),
                      winslash = "/"))
source("abc/abc_trajectories.R")

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

pop_y1 <- matrix(4 * c( 0,  0.438, -0.035,
                         0,  0.000,  0.000,
                         0, -0.438,  0.035),
                 nrow = 3, ncol = 3, byrow = TRUE)

pop_y2 <- matrix(4 * c(0.214,  0.125,  0.010,
                        0.000, -0.027,  0.002,
                        1.688,  0.125, -0.010),
                 nrow = 3, ncol = 3, byrow = TRUE)

# Test 1: ABC = 0 exactly when pop == est
r1 <- compute_abc(pop_y1, pop_y1, t_min = 0, t_max = 7)
check("T1  ABC = 0 when pop == est", abs(r1$abc_mean) < 1e-10)
check("T1b all abc_by_class = 0", all(r1$abc_by_class == 0))

# Test 2: Known constant offset  ABC = delta * (t_max - t_min)
pop_c <- matrix(0, 3, 3)
est_c <- pop_c
est_c[, 1] <- 2          # shift intercept by delta = 2
r2 <- compute_abc(pop_c, est_c, t_min = 0, t_max = 7)
check("T2  Constant offset ABC = 14", abs(r2$abc_mean - 14) < 1e-10)

# Test 3: Permutation recovery after 1 <-> 3 swap
est_swap <- pop_y1[c(3, 2, 1), ]
r3 <- compute_abc(pop_y1, est_swap, t_min = 0, t_max = 7)
check("T3  Optimal permutation = c(3,2,1)", all(r3$optimal_permutation == c(3L, 2L, 1L)))

# Test 4: Optimal < identity
check("T4  Optimal ABC < identity ABC",
      r3$abc_mean < r3$all_permutations_abc["(1,2,3)"])

# Test 5: Symmetry
r5a <- compute_abc(pop_y1, est_swap)
r5b <- compute_abc(est_swap, pop_y1)
check("T5  Symmetry: ABC(A,B) == ABC(B,A)", abs(r5a$abc_mean - r5b$abc_mean) < 1e-10)

# Test 6: Relabeling of class_assignments
est_assignments  <- c(3L, 1L, 2L, 1L, 3L)
expected_relabel <- c(1L, 3L, 2L, 3L, 1L)
r6 <- compute_abc(pop_y1, est_swap, class_assignments = est_assignments)
check("T6  Class relabeling correct", all(r6$relabeled_assignments == expected_relabel))

# Test 7: K! permutation counts
p3 <- .all_permutations(3L)
p4 <- .all_permutations(4L)
check("T7a .all_permutations(3) has 6 entries", length(p3) == 6L)
check("T7b .all_permutations(4) has 24 entries", length(p4) == 24L)
check("T7c K=3 permutations are all unique",
      length(unique(vapply(p3, paste, character(1), collapse = ","))) == 6L)
check("T7d K=3 permutations cover 1:3",
      all(vapply(p3, function(p) all(sort(p) == 1:3), logical(1))))

# Test 8: Two-outcome structure
r8 <- compute_abc(
  pop_betas = list(Y1 = pop_y1, Y2 = pop_y2),
  est_betas = list(Y1 = pop_y1[c(3,2,1), ], Y2 = pop_y2[c(3,2,1), ]),
  t_min = 0, t_max = 7
)
check("T8a Two outcomes: length(abc_by_outcome) == 2", length(r8$abc_by_outcome) == 2L)
check("T8b Two outcomes: length(abc_by_class) == 3",  length(r8$abc_by_class) == 3L)
check("T8c Grand mean == mean(abc_by_outcome)",
      abs(r8$abc_mean - round(mean(r8$abc_by_outcome), 4)) < 1e-4)

# Test 9: Input validation
err_caught <- function(expr) tryCatch({expr; FALSE}, error = function(e) TRUE)
check("T9a t_min >= t_max raises error",
      err_caught(compute_abc(pop_y1, pop_y1, t_min = 7, t_max = 0)))
check("T9b Mismatched K raises error",
      err_caught(compute_abc(pop_y1, pop_y1[1:2, ])))
check("T9c Zero-column matrix raises error",
      err_caught(compute_abc(matrix(numeric(0), nrow = 3L, ncol = 0L),
                             matrix(numeric(0), nrow = 3L, ncol = 0L))))
check("T9d Mismatched outcome list length raises error",
      err_caught(compute_abc(list(Y1=pop_y1, Y2=pop_y2), list(Y1=pop_y1))))

# Test 10: Linear difference  d(t) = -3.5 + t  on [0,7]  => ABC = 12.25
pop_lin <- matrix(c(0, 0, 0), nrow = 1, ncol = 3)
est_lin <- matrix(c(3.5, -1, 0), nrow = 1, ncol = 3)
r10 <- compute_abc(pop_lin, est_lin, t_min = 0, t_max = 7)
check("T10 Linear difference ABC = 12.25", abs(r10$abc_mean - 12.25) < 1e-12)

# Test 11: Cubic trajectory  d(t) = -t^3  on [0,7]  => ABC = 7^4/4 = 600.25
# pop = 0 (flat), est has pure cubic coefficient = 1; no interior sign change
pop_cub <- matrix(c(0, 0, 0, 0), nrow = 1, ncol = 4)
est_cub <- matrix(c(0, 0, 0, 1), nrow = 1, ncol = 4)
r11 <- compute_abc(pop_cub, est_cub, t_min = 0, t_max = 7)
check("T11 Cubic difference ABC = 600.25", abs(r11$abc_mean - 600.25) < 1e-10)

# Test 12: Zero-padding equivalence, adding a zero column should not change ABC
pop_y1_4  <- cbind(pop_y1,  0)
est_swap_4 <- cbind(pop_y1[c(3, 2, 1), ], 0)
r12 <- compute_abc(pop_y1_4, est_swap_4, t_min = 0, t_max = 7)
check("T12 Zero-padded K×4 gives same ABC as K×3",
      abs(r12$abc_mean - r3$abc_mean) < 1e-10)

cat("\n============================\n")
cat(sprintf("Results: %d passed, %d failed\n", pass, fail))
