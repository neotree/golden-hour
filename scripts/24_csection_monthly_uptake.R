# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# AUGUST 2026 REQUESTS -- Script 24: DCC / ESSC UPTAKE BY MONTH, CAESAREAN ONLY
# =============================================================================
# Rachel, 28 July 2026 (point 2): "a table of by-month uptake number and % of
#   eligible babies born by C-section (emergency and elective can be combined)
#   for DCC and ESSC."
#
# This is the Caesarean-restricted analogue of Tables 9c/9d (Script 09), so the
# eligibility definition, status coding and percentage definitions are kept
# IDENTICAL to Script 09 for comparability with the tables Rachel already has.
#
# POPULATION : maternal outcome script (Prisca's delivery log), facility SMCH,
#              live births (neotreeoutcome LB or ENND), delivered by Caesarean
#              section (modedelivery 4 or 5 -- elective and emergency combined,
#              as requested).
# ELIGIBLE   : apgar1 > 6 OR resus == "N"   (Script 09 eligibility proxy).
#              NOTE: this is the Script 08/09 definition, NOT the resus-only
#              definition used for Table 1 (Script 23). Both are reported in the
#              log so the difference is visible.
# STATUS     : dcc / s2s coded Y = received, N = not received, U = unknown,
#              "" / NA = not recorded.
# PERCENTAGES (per month, Caesarean eligible babies):
#   pct_yes          = Y / n_eligible_cs           <- the coverage figure
#   pct_yes_of_known = Y / (Y + N)                 <- the "of documented" figure
#   pct_unknown      = (U + not recorded) / n_eligible_cs
#
# WINDOWS: both are produced --
#   "12mo"     1 Nov 2024 - 31 Oct 2025  (the key study period; use this one)
#   "extended" 1 Nov 2024 - 31 Mar 2026  (as in the earlier 9c/9d tables)
#
# DSH note: ASCII-only.  Data file: zim_db_maternal_outcomes_20260501_cleaned.csv
# =============================================================================

library(tidyverse)
library(lubridate)

args <- commandArgs(trailingOnly = FALSE)
sf <- args[grep("--file=", args)]
SCRIPT_DIR <- if (length(sf) > 0)
  dirname(normalizePath(sub("--file=", "", sf), mustWork = FALSE)) else getwd()

MATERNAL <- normalizePath(file.path(SCRIPT_DIR, "..", "00-DATA",
  "zim_db_maternal_outcomes_20260501_cleaned.csv"), mustWork = FALSE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", "maternal_august")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "24_csection_monthly_uptake_log.txt")
log_con <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
END_12MO  <- as.Date("2025-10-31")
END_EXT   <- as.Date("2026-03-31")

# Read every column as character (see Script 22/23 header: readr's type guesser
# mangles dateadmission for the Nov 2025 - Mar 2026 records).
raw <- read_csv(MATERNAL, show_col_types = FALSE,
                col_types = cols(.default = col_character()))

df <- raw %>%
  mutate(
    adm_date  = as.Date(dateadmission),
    month_label = format(adm_date, "%Y-%m"),
    live_birth = neotreeoutcome %in% c("LB", "ENND"),
    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    # Script 09 eligibility proxy
    eligible = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus == "N"                        ~ TRUE,
      TRUE                                ~ FALSE),
    eligible_resus_only = resus == "N" & !is.na(resus),   # Table 1 definition
    caesarean = modedelivery %in% c("4", "5"),
    dcc_Y = dcc == "Y", dcc_N = dcc == "N", dcc_U = dcc == "U",
    dcc_NR = is.na(dcc) | trimws(dcc) == "",
    s2s_Y = s2s == "Y", s2s_N = s2s == "N", s2s_U = s2s == "U",
    s2s_NR = is.na(s2s) | trimws(s2s) == "")

sum_na <- function(x) sum(x, na.rm = TRUE)

monthly_cs <- function(end_date, elig_var) {
  d <- df %>% filter(facility == "SMCH", live_birth == TRUE, caesarean == TRUE,
                     .data[[elig_var]] == TRUE,
                     adm_date >= INT_START, adm_date <= end_date)
  body <- d %>% group_by(month_label) %>% summarise(
    n_eligible_cs = n(),
    dcc_Y = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR),
    .groups = "drop")
  tot <- d %>% summarise(month_label = "TOTAL",
    n_eligible_cs = n(),
    dcc_Y = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR))
  bind_rows(arrange(body, month_label), tot) %>%
    mutate(
      dcc_pct_yes          = round(100 * dcc_Y / n_eligible_cs, 1),
      dcc_pct_yes_of_known = if_else(dcc_Y + dcc_N > 0,
                                     round(100 * dcc_Y / (dcc_Y + dcc_N), 1), NA_real_),
      dcc_pct_unknown      = round(100 * (dcc_U + dcc_not_recorded) / n_eligible_cs, 1),
      essc_pct_yes          = round(100 * essc_Y / n_eligible_cs, 1),
      essc_pct_yes_of_known = if_else(essc_Y + essc_N > 0,
                                      round(100 * essc_Y / (essc_Y + essc_N), 1), NA_real_),
      essc_pct_unknown      = round(100 * (essc_U + essc_not_recorded) / n_eligible_cs, 1)) %>%
    relocate(month_label, n_eligible_cs,
             dcc_Y, dcc_pct_yes, dcc_pct_yes_of_known, dcc_N, dcc_U, dcc_not_recorded, dcc_pct_unknown,
             essc_Y, essc_pct_yes, essc_pct_yes_of_known, essc_N, essc_U, essc_not_recorded, essc_pct_unknown)
}

for (w in list(list(tag = "12mo", end = END_12MO), list(tag = "extended", end = END_EXT))) {
  tb <- monthly_cs(w$end, "eligible")
  cat(sprintf("\n=========== CAESAREAN DCC/ESSC UPTAKE BY MONTH (%s: %s to %s) ===========\n",
              w$tag, INT_START, w$end))
  print(as.data.frame(tb), row.names = FALSE)
  write_csv(tb, file.path(OUTPUT_DIR, sprintf("24a_csection_monthly_uptake_%s.csv", w$tag)))
}

# ---- Context for the footnotes ---------------------------------------------
cat("\n-------------------------------------------------------------\n")
cat("  CONTEXT / SENSITIVITY\n")
cat("-------------------------------------------------------------\n")
ctx <- df %>% filter(facility == "SMCH", live_birth == TRUE,
                     adm_date >= INT_START, adm_date <= END_12MO)
cat(sprintf("12-month window, SMCH live births                  : %d\n", nrow(ctx)))
cat(sprintf("  of which Caesarean (modedelivery 4 or 5)         : %d\n", sum_na(ctx$caesarean)))
cat(sprintf("    modedelivery == 4                              : %d\n",
            sum(ctx$modedelivery == "4", na.rm = TRUE)))
cat(sprintf("    modedelivery == 5                              : %d\n",
            sum(ctx$modedelivery == "5", na.rm = TRUE)))
cat(sprintf("  Caesarean AND eligible (apgar1>6 OR resus=N)     : %d\n",
            sum(ctx$caesarean & ctx$eligible, na.rm = TRUE)))
cat(sprintf("  Caesarean AND eligible (resus=N only, Table 1 def): %d\n",
            sum(ctx$caesarean & ctx$eligible_resus_only, na.rm = TRUE)))

# Same table under the Table 1 (resus-only) eligibility, for cross-checking
# against Table 1's Caesarean row.
tb_alt <- monthly_cs(END_12MO, "eligible_resus_only")
write_csv(tb_alt, file.path(OUTPUT_DIR, "24b_csection_monthly_uptake_12mo_resusonly.csv"))
cat("\nAlso saved 24b (same table, resus-only eligibility) to reconcile with Table 1.\n")
cat("TOTAL row of 24b:\n")
print(as.data.frame(tb_alt %>% filter(month_label == "TOTAL")), row.names = FALSE)

cat("\n=== Script 24 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
