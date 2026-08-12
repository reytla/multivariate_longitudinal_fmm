# =============================================================================
# requirements.R
# -----------------------------------------------------------------------------
# Installs all R packages required by this repository.
# Run once before using any scripts:
#
#   source("requirements.R")
#
# Tested with R 4.5.2 on Windows. Package versions used during development
# are noted below; newer versions should generally work.
#
# For fully reproducible environments (pinned transitive dependencies) use
# renv instead:
#   install.packages("renv")
#   renv::restore()          # restores from renv.lock if present
#   renv::init()             # or initialise a new lockfile
# =============================================================================

pkgs <- c(
  # Data generation and Mplus interfacing
  "MplusAutomation",   # 1.1.1: generate .inp files, run Mplus, extract output
  "MASS",              # 7.3.64: mvrnorm() for multivariate normal simulation
  "stringr",           # 1.5.1:  string helpers for parsing Mplus output titles
  "plyr",              # 1.8.9:  ldply() for collating results across replications

  # Testing (abc/tests/ only, not needed for simulation pipeline)
  "testthat"           # 3.2.3
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(to_install) == 0) {
  message("All packages already installed.")
} else {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
}

# Report installed versions
cat("\nInstalled package versions:\n")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT FOUND")
  cat(sprintf("  %-20s %s\n", p, v))
}
cat(sprintf("  %-20s %s\n", "R", R.version.string))
