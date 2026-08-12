# =============================================================================
# run_models.R
# -----------------------------------------------------------------------------
# Generates Mplus .inp files from MplusAutomation templates and runs them.
#
# Use this script for the MC (Monte Carlo classification) stage, where the
# correct number of classes K is known and the focus is on parameter recovery,
# class assignment accuracy, and trajectory recovery.
#
# The DGM scripts (GBMTM_DGM.R, MLCGA_DGM.R, MCPGMM_DGM.R) already call
# createModels() and runModels() internally. This standalone script is useful
# when you want to re-run or extend a model set without re-generating data.
#
# Requirements: MplusAutomation
# =============================================================================

library(MplusAutomation)

# =============================================================================
# CONFIGURATION
# =============================================================================

wd <- "PATH/TO/YOUR/WORKING/DIRECTORY"   # Must contain a /Templates subfolder

# replaceOutfile options:
#   "never"        : skip any model that already has a .out file
#   "modifiedDate" : re-run only if the .inp is newer than the .out
#   "always"       : always re-run all models
replace_rule <- "modifiedDate"

# =============================================================================
# CREATE .inp FILES FROM TEMPLATES
# =============================================================================

template_files <- list.files(file.path(wd, "Templates"), full.names = FALSE)

for (j in template_files) {
  createModels(file.path(wd, "Templates", j))
}

# =============================================================================
# RUN ALL MODELS
# =============================================================================

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

runModels(wd,
          recursive      = TRUE,
          showOutput     = TRUE,
          replaceOutfile = replace_rule,
          logFile        = file.path(wd, "ComparisonLog.txt"),
          Mplus_command  = .mplus_cmd)

message("All models complete. Log: ", file.path(wd, "ComparisonLog.txt"))
