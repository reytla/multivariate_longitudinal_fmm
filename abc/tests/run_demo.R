setwd("h:/.../simulation_code")
source("abc/abc_trajectories.R")

pop_y1 <- matrix(4 * c( 0,  0.438, -0.035,
                         0,  0.000,  0.000,
                         0, -0.438,  0.035),
                 nrow = 3, ncol = 3, byrow = TRUE)
pop_y2 <- matrix(4 * c(0.214,  0.125,  0.010,
                        0.000, -0.027,  0.002,
                        1.688,  0.125, -0.010),
                 nrow = 3, ncol = 3, byrow = TRUE)

est_y1 <- pop_y1[c(3, 2, 1), ]
est_y2 <- pop_y2[c(3, 2, 1), ]

abc_res <- compute_abc(
  pop_betas         = list(Y1 = pop_y1, Y2 = pop_y2),
  est_betas         = list(Y1 = est_y1, Y2 = est_y2),
  class_assignments = c(1L, 3L, 2L, 1L, 2L, 3L)
)
print_abc(abc_res)
