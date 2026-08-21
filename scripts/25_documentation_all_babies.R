# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# AUGUST 2026 REQUESTS -- Script 25: DCC / ESSC DOCUMENTATION, ALL LIVE BIRTHS
# =============================================================================
# Rachel, 28 July 2026 (point 3): table 8b gives, per month, live births,
#   eligible babies, and Y / N / Unknown for DCC and ESSC -- but the Y/N/U
#   counts add up only to the ELIGIBLE babies. She would also like the Y / N /
#   Unknown documentation for ALL babies: "Every baby should have a status
#   recorded even if they were not eligible. Then the midwives should have put
#   N - if that makes sense."
#
# So this script reproduces the 8b layout with the status breakdown given TWICE
# per month: once for ALL live births and once for ELIGIBLE live births (the
# original 8b columns), so the two are directly comparable and the ineligible
# remainder is visible.
#
# POPULATION : maternal outcome script (Prisca's delivery log), facility SMCH,
#              live births (neotreeoutcome LB or ENND).
# ELIGIBLE   : apgar1 > 6 OR resus == "N"  -- the Script 08/09 proxy, kept
#              identical to the original table 8b.
# STATUS     : Y = received, N = not received, U = unknown,
#              "" / NA = not recorded (reported separately from U, because
#              "midwife wrote nothing" and "midwife wrote unknown" are different
#              documentation failures).
#
# WINDOWS: both are produced --
#   "12mo"     1 Nov 2024 - 31 Oct 2025  (the key study period; use this one)
#   "extended" 1 Nov 2024 - 31 Mar 2026  (as in the earlier 8b table)
#
# OUTPUTS: 25a_documentation_all_vs_eligible_<window>.csv  (wide, 8b-style)
#          25b_documentation_ineligible_<window>.csv       (the remainder only)
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
LOG_FILE <- file.path(OUTPUT_DIR, "25_documentation_all_babies_log.txt")
log_con <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
END_12MO  <- as.Date("2025-10-31")
END_EXT   <- as.Date("2026-03-31")

raw <- read_csv(MATERNAL, show_col_types = FALSE,
                col_types = cols(.default = col_character()))

df <- raw %>%
  mutate(
    adm_date = as.Date(dateadmission),
    month_label = format(adm_date, "%Y-%m"),
    live_birth = neotreeoutcome %in% c("LB", "ENND"),
    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    eligible = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus == "N"                        ~ TRUE,
      TRUE                                ~ FALSE),
    dcc_Y = dcc == "Y", dcc_N = dcc == "N", dcc_U = dcc == "U",
    dcc_NR = is.na(dcc) | trimws(dcc) == "",
    s2s_Y = s2s == "Y", s2s_N = s2s == "N", s2s_U = s2s == "U",
    s2s_NR = is.na(s2s) | trimws(s2s) == "")

sum_na <- function(x) sum(x, na.rm = TRUE)

# status breakdown for an arbitrary subset, with a column-name suffix
status_block <- function(d, suffix) {
  out <- d %>% group_by(month_label) %>% summarise(
    n            = n(),
    dcc_Y        = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U        = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y       = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U       = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR),
    .groups = "drop")
  tot <- d %>% summarise(month_label = "TOTAL", n = n(),
    dcc_Y = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR))
  out <- bind_rows(arrange(out, month_label), tot) %>%
    mutate(dcc_pct_documented  = round(100 * (dcc_Y + dcc_N) / n, 1),
           essc_pct_documented = round(100 * (essc_Y + essc_N) / n, 1))
  names(out)[names(out) != "month_label"] <-
    paste0(names(out)[names(out) != "month_label"], "_", suffix)
  out
}

build <- function(end_date, tag) {
  base <- df %>% filter(facility == "SMCH", live_birth == TRUE,
                        adm_date >= INT_START, adm_date <= end_date)
  all_b  <- status_block(base, "all")
  elig_b <- status_block(base %>% filter(eligible == TRUE), "eligible")
  inel_b <- status_block(base %>% filter(eligible == FALSE), "ineligible")

  wide <- all_b %>%
    left_join(elig_b, by = "month_label") %>%
    left_join(inel_b %>% select(month_label, n_ineligible), by = "month_label") %>%
    rename(n_live_births = n_all, n_eligible = n_eligible) %>%
    relocate(month_label, n_live_births, n_eligible, n_ineligible)

  cat(sprintf("\n===== DCC/ESSC DOCUMENTATION, ALL vs ELIGIBLE LIVE BIRTHS (%s: %s to %s) =====\n",
              tag, INT_START, end_date))
  cat("\n-- All live births --\n")
  print(as.data.frame(wide %>% select(month_label, n_live_births, n_eligible, n_ineligible,
                                      ends_with("_all"))), row.names = FALSE)
  cat("\n-- Eligible live births (original table 8b columns) --\n")
  print(as.data.frame(wide %>% select(month_label, n_eligible, ends_with("_eligible"))),
        row.names = FALSE)
  write_csv(wide, file.path(OUTPUT_DIR,
            sprintf("25a_documentation_all_vs_eligible_%s.csv", tag)))
  write_csv(inel_b, file.path(OUTPUT_DIR,
            sprintf("25b_documentation_ineligible_%s.csv", tag)))
  wide
}

w12 <- build(END_12MO, "12mo")
wex <- build(END_EXT,  "extended")

# ---- Headline numbers for the covering note --------------------------------
tot12 <- w12 %>% filter(month_label == "TOTAL")
cat("\n-------------------------------------------------------------\n")
cat("  HEADLINE (12-month window)\n")
cat("-------------------------------------------------------------\n")
cat(sprintf("  Live births                         : %d\n", tot12$n_live_births))
cat(sprintf("  Eligible                            : %d\n", tot12$n_eligible))
cat(sprintf("  Ineligible                          : %d\n", tot12$n_ineligible))
cat(sprintf("  DCC  documented (Y or N), all babies: %d (%.1f%%)\n",
            tot12$dcc_Y_all + tot12$dcc_N_all, tot12$dcc_pct_documented_all))
cat(sprintf("  DCC  documented (Y or N), eligible  : %d (%.1f%%)\n",
            tot12$dcc_Y_eligible + tot12$dcc_N_eligible, tot12$dcc_pct_documented_eligible))
cat(sprintf("  ESSC documented (Y or N), all babies: %d (%.1f%%)\n",
            tot12$essc_Y_all + tot12$essc_N_all, tot12$essc_pct_documented_all))
cat(sprintf("  ESSC documented (Y or N), eligible  : %d (%.1f%%)\n",
            tot12$essc_Y_eligible + tot12$essc_N_eligible, tot12$essc_pct_documented_eligible))

cat("\n=== Script 25 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
