# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# run_all.R -- reproduce every table and figure in the manuscript
# =============================================================================
# Usage:
#   1. Place the 3 required data files in input/ (see input/README.md).
#   2. From the repo root: Rscript run_all.R
#   3. Every output lands under output/ -- see FIGURE_TABLE_SOURCES.md for
#      which file corresponds to which manuscript table/figure.
#
# HOW THIS WORKS
#   Each script in scripts/ was written for its own original location in the
#   analysis working folder and resolves its data/output paths relative to
#   itself using one of two conventions:
#     (a) 01, 03, 08            -> reads  ../00-DATA/<file>
#                                -> writes ../03-OUTPUTS/
#     (b) 16,17,18,19,20,23,24,25,27 (August scripts)
#                                -> reads  ../00-DATA/<file>
#                                -> writes ./outputs/<subfolder>/
#   Rather than edit 12 scripts by hand (two path dialects, easy to typo one),
#   this file creates 3 symlinks so every script finds exactly the sibling
#   directory it already expects, pointing at input/ and output/. No script
#   in scripts/ is modified.
# =============================================================================

suppressWarnings(suppressMessages({

REPO_ROOT   <- normalizePath(dirname(sub("--file=", "",
                  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))),
                  mustWork = FALSE)
if (length(REPO_ROOT) == 0 || !nzchar(REPO_ROOT)) REPO_ROOT <- getwd()

SCRIPTS_DIR <- file.path(REPO_ROOT, "scripts")
INPUT_DIR   <- file.path(REPO_ROOT, "input")
OUTPUT_DIR  <- file.path(REPO_ROOT, "output")

}))

cat("=============================================================\n")
cat("  Golden Hour -- run_all.R\n")
cat("=============================================================\n")
cat(sprintf("Repo root : %s\n", REPO_ROOT))

# -----------------------------------------------------------------------------
# 0. Packages -- install anything missing, then load
# -----------------------------------------------------------------------------
REQUIRED_PACKAGES <- c("tidyverse", "lubridate", "scales", "broom",
                        "officer", "flextable")

cat("\n-- Checking packages --------------------------------------\n")
missing_pkgs <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace,
                                           logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("Installing missing package(s): %s\n",
              paste(missing_pkgs, collapse = ", ")))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(REQUIRED_PACKAGES, function(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))))
cat("All packages loaded.\n")

# -----------------------------------------------------------------------------
# 1. Check input/ has the 3 required data files
# -----------------------------------------------------------------------------
REQUIRED_DATA_FILES <- c(
  "ZIM_db_master_joined_to_20260401.csv",
  "ZIM_db_master_joined_to_20260525.csv",
  "zim_db_maternal_outcomes_20260501_cleaned.csv"
)

cat("\n-- Checking input/ -------------------------------------------\n")
present <- file.exists(file.path(INPUT_DIR, REQUIRED_DATA_FILES))
if (!all(present)) {
  missing_files <- REQUIRED_DATA_FILES[!present]
  stop(
    "\n\n",
    "=============================================================\n",
    "  MISSING DATA FILES\n",
    "=============================================================\n",
    "run_all.R cannot proceed -- the following file(s) are not in\n",
    "  ", INPUT_DIR, "\n\n",
    paste0("  - ", missing_files, collapse = "\n"), "\n\n",
    "These files are confidential and are not included in this\n",
    "repository. See input/README.md for where to obtain them.\n",
    "=============================================================\n",
    call. = FALSE
  )
}
cat("All 3 required data files found.\n")

# -----------------------------------------------------------------------------
# 2. Wire up the path redirection (see header comment) -- symlinks, with a
#    directory-copy fallback if the filesystem/OS doesn't support them
#    (e.g. Windows without Developer Mode / admin rights).
# -----------------------------------------------------------------------------
cat("\n-- Setting up input/output redirection ------------------------\n")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "03-OUTPUTS"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "outputs"),    showWarnings = FALSE, recursive = TRUE)

link_or_copy <- function(target, link_path) {
  # remove any stale link/dir/file from a previous run before recreating
  existing_link <- Sys.readlink(link_path)
  if (!is.na(existing_link) && nzchar(existing_link)) {
    unlink(link_path)
  } else if (dir.exists(link_path) || file.exists(link_path)) {
    unlink(link_path, recursive = TRUE)
  }
  ok <- tryCatch(isTRUE(file.symlink(target, link_path)), error = function(e) FALSE)
  if (!isTRUE(ok)) {
    cat(sprintf("  (symlink unavailable, copying instead: %s)\n", basename(link_path)))
    dir.create(link_path, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(target, full.names = TRUE), link_path, recursive = TRUE)
  }
  ok
}

invisible(link_or_copy(INPUT_DIR,                          file.path(REPO_ROOT, "00-DATA")))
invisible(link_or_copy(file.path(OUTPUT_DIR, "03-OUTPUTS"), file.path(REPO_ROOT, "03-OUTPUTS")))
invisible(link_or_copy(file.path(OUTPUT_DIR, "outputs"),    file.path(SCRIPTS_DIR, "outputs")))

cat("Redirection ready: 00-DATA/, 03-OUTPUTS/ -> output/03-OUTPUTS/, scripts/outputs/ -> output/outputs/\n")

# -----------------------------------------------------------------------------
# 3. Run every script, in dependency order (20 depends on 18's output; all
#    others are independent). Continue past a failing script; report at end.
# -----------------------------------------------------------------------------
SCRIPTS_IN_ORDER <- c(
  "01_descriptive_overview.R",
  "03_temperature_trends.R",
  "08_maternal_script_uptake.R",
  "16_data_readiness_checks.R",
  "17_eligible_cohort_tables.R",
  "18_outcome_or_and_sensitivity.R",
  "19_reference_group_sensitivity.R",
  "20_forest_table4.R",
  "23_table1_maternal_uptake_12mo.R",
  "24_csection_monthly_uptake.R",
  "25_documentation_all_babies.R",
  "27_lbw_monthly_uptake.R"
)

cat("\n=============================================================\n")
cat("  Running scripts\n")
cat("=============================================================\n")

results <- vector("list", length(SCRIPTS_IN_ORDER))
names(results) <- SCRIPTS_IN_ORDER

for (script in SCRIPTS_IN_ORDER) {
  cat(sprintf("\n--- Running %s ---\n", script))
  script_path <- file.path(SCRIPTS_DIR, script)
  status <- tryCatch({
    # each script is self-contained; run in a fresh child R process so one
    # script's environment/options can't leak into the next
    res <- system2("Rscript", shQuote(script_path), stdout = TRUE, stderr = TRUE)
    exit_code <- attr(res, "status")
    if (!is.null(exit_code) && exit_code != 0) {
      cat(paste(tail(res, 20), collapse = "\n"), "\n")
      stop(sprintf("exited with non-zero status (%s)", exit_code))
    }
    "OK"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  results[[script]] <- status
  cat(sprintf("--- %s: %s ---\n", script,
              if (identical(status, "OK")) "OK" else status))
}

# -----------------------------------------------------------------------------
# 4. Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("  SUMMARY\n")
cat("=============================================================\n")
ok_scripts   <- names(results)[vapply(results, identical, logical(1), "OK")]
fail_scripts <- setdiff(names(results), ok_scripts)

cat(sprintf("Succeeded (%d/%d):\n", length(ok_scripts), length(results)))
for (s in ok_scripts) cat(sprintf("  [OK]   %s\n", s))

if (length(fail_scripts) > 0) {
  cat(sprintf("\nFailed (%d/%d):\n", length(fail_scripts), length(results)))
  for (s in fail_scripts) cat(sprintf("  [FAIL] %s -- %s\n", s, results[[s]]))
  cat("\nSee each script's own *_log.txt under output/ for full detail.\n")
} else {
  cat("\nAll scripts completed successfully.\n")
}
cat(sprintf("\nOutputs written to: %s\n", OUTPUT_DIR))
cat("=============================================================\n")
