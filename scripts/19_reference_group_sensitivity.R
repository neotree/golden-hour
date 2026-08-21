# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# JULY 2026 REQUESTS -- Script 19: REFERENCE-GROUP SENSITIVITY (cross-file)
# =============================================================================
# Backs the reference-group figures quoted in the manuscript Limitations
# paragraph, so every number in the outcomes document traces to an R script
# (Tables 2/3 come from Script 17; Tables 4/5 from Script 18; cohort flow from
# Scripts 16/17; this script supplies the reference-group sensitivity numbers).
#
# It reproduces, for BOTH cleaned master extracts:
#   - eligible cohort (apgar1 < 7 excluded, Aug 2026 rule), delivinter-known N and %,
#   - the Neither (no DCC, no S2S) group among matched infants: n, deaths, rate.
# and then the "recovered" reference infants = new-file Neither minus old-file
# Neither (valid because the exposed groups are identical across the two files),
# with their mortality. These give: delivinter-known 37.6% -> 57.7%; Neither
# crude mortality 9.3% -> 8.2%; recovered infants ~7.6%.
#
# Eligibility / exposure / outcome definitions are IDENTICAL to Scripts 17/18
# (team decision 21 Jul: exclude 5-minute Apgar < 7; resus limb dropped).
#
# DSH note: ASCII-only.
# =============================================================================

library(tidyverse)
library(lubridate)

args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
SCRIPT_DIR <- if (length(script_flag) > 0)
  dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE)) else getwd()

# Both cleaned master extracts (new first, then the earlier one).
DATA_FILES <- c("ZIM_db_master_joined_to_20260525.csv",
                "ZIM_db_master_joined_to_20260401.csv")
DATA_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", "00-DATA"), mustWork = FALSE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", "reference_group_apgar1")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "19_reference_group_sensitivity_log.txt")
log_con  <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")
APGAR_CUTOFF <- 7

cat("=============================================================\n")
cat("  Script 19: Reference-group sensitivity (cross-file)\n")
cat("=============================================================\n\n")

# Process one master file -> one-row summary of eligible/known/Neither counts.
summarise_file <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  # DATE-PARSING FIX (Aug 2026): read_csv's type guesser converts
# "datetimeadmission" to POSIXct; ymd_hms() then re-coerces it via as.character(),
# which for a midnight timestamp yields a date-only string ("2025-08-08") that
# ymd_hms cannot parse -> NA -> the record is silently dropped by the window
# filter. Four SMCH inborn admissions were lost this way in the July run
# (4,856 instead of 4,860). Read the column as character and take the date part.
raw <- read_csv(path, show_col_types = FALSE, guess_max = 100000,
                  col_types = cols(datetimeadmission = col_character()))

  df <- raw %>%
    mutate(
      adm_date = as.Date(substr(datetimeadmission, 1, 10)),
      inborn = inorout %in% c("Yes","true","True","In"),
      is_readmission = readmission == "Y",
      matched = match_type == "direct_match",
      apgar1_num = suppressWarnings(as.numeric(apgar1)),
      apgar1_num = if_else(!is.na(apgar1_num) & (apgar1_num < 0 | apgar1_num > 10),
                           NA_real_, apgar1_num),
      excl_apgar = !is.na(apgar1_num) & apgar1_num < APGAR_CUTOFF,
      di_known = !(is.na(delivinter) | trimws(delivinter) == "" | trimws(delivinter) == "U"),
      has_dcc  = di_known & grepl("DCC", delivinter),
      has_s2s  = di_known & grepl("S2S", delivinter),
      is_neither = di_known & !has_dcc & !has_s2s,
      nnd = neotreeoutcome == "NND")

  pop <- df %>%
    filter(facility == "SMCH", inborn == TRUE,
           is_readmission == FALSE | is.na(is_readmission),
           adm_date >= INT_START, adm_date <= INT_END)
  elig <- pop %>% filter(excl_apgar == FALSE)
  elig_known <- elig %>% filter(di_known == TRUE)
  neither_m  <- elig_known %>% filter(is_neither == TRUE, matched == TRUE)

  tibble(
    file = fname,
    eligible = nrow(elig),
    delivinter_known = nrow(elig_known),
    delivinter_known_pct = round(100 * nrow(elig_known) / nrow(elig), 1),
    neither_matched_n = nrow(neither_m),
    neither_deaths = sum(neither_m$nnd, na.rm = TRUE),
    neither_rate_pct = round(100 * sum(neither_m$nnd, na.rm = TRUE) / max(nrow(neither_m), 1), 1))
}

res <- map_dfr(DATA_FILES, summarise_file)
cat("\nPer-file summary:\n")
print(as.data.frame(res), row.names = FALSE)

new <- res[1, ]; old <- res[2, ]
rec_n <- new$neither_matched_n - old$neither_matched_n
rec_d <- new$neither_deaths    - old$neither_deaths
rec_rate <- round(100 * rec_d / max(rec_n, 1), 1)

cat("\n-------------------------------------------------------------\n")
cat("  REFERENCE-GROUP SENSITIVITY (manuscript Limitations)\n")
cat("-------------------------------------------------------------\n")
cat(sprintf("  delivinter-known, %% of eligible : %.1f%% (old) -> %.1f%% (new)\n",
            old$delivinter_known_pct, new$delivinter_known_pct))
cat(sprintf("  Neither group (matched)         : %d (old) -> %d (new)  [x%.1f]\n",
            old$neither_matched_n, new$neither_matched_n,
            new$neither_matched_n / max(old$neither_matched_n, 1)))
cat(sprintf("  Neither crude mortality         : %.1f%% (old) -> %.1f%% (new)\n",
            old$neither_rate_pct, new$neither_rate_pct))
cat(sprintf("  Recovered reference infants     : n=%d, deaths=%d, mortality=%.1f%%\n",
            rec_n, rec_d, rec_rate))
cat("  (Recovered = new-file Neither minus old-file Neither; valid because the\n")
cat("   exposed groups are identical across the two files.)\n")

out <- res %>% mutate(label = c("new (20260525)", "old (20260401)")) %>%
  bind_rows(tibble(file = "recovered (new - old)", label = "recovered",
                   eligible = NA, delivinter_known = NA, delivinter_known_pct = NA,
                   neither_matched_n = rec_n, neither_deaths = rec_d,
                   neither_rate_pct = rec_rate))
write_csv(out, file.path(OUTPUT_DIR, "19a_reference_group_sensitivity.csv"))

cat("\nSaved: 19a_reference_group_sensitivity.csv\n")
cat("\n=== Script 19 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
