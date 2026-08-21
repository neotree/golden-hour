# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# AUGUST 2026 REQUESTS -- Script 23: TABLE 1 restricted to the 12-month period
# =============================================================================
# Rachel, 31 July 2026 (point 4): Table 1 must cover the KEY STUDY PERIOD ONLY,
#   1 November 2024 - 31 October 2025, to match the rest of the paper. This
#   supersedes Script 22 (which ran 1 Nov 2024 - 31 Mar 2026).
#
# CHANGE FROM SCRIPT 22:
#   - INT_END moved from 2026-03-31 to 2025-10-31.
#   - The "Maternal script >12 months" Yes/No block is DROPPED: with the window
#     restricted, every birth falls in "No", so the row carries no information.
#     (Agreed with David, 6 Aug 2026.)
#
# ELIGIBILITY (Rachel's definition for this table): "Did this baby require
#   immediate resuscitation = No", i.e. resus == "N" in the maternal script.
# POPULATION: SMCH live births (neotreeoutcome LB or ENND) in the maternal
#   outcome script (Prisca's delivery log).
# COLUMNS: Eligible infants n; Received DCC n (%); Received ESSC n (%);
#   Received both n (%). Percentages are of the row's eligible n (coverage,
#   i.e. Unknown/not-recorded records count in the denominator).
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
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", "table1_12mo")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "23_table1_maternal_uptake_12mo_log.txt")
log_con <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")   # <-- 31 Jul 2026 change (was 2026-03-31)

# NOTE (carried over from Script 22): read EVERY column as character. readr's
# type guesser infers "dateadmission" as datetime from the first block of rows,
# then silently mangles the Nov 2025 - Mar 2026 records (all timestamped
# 22:00:00). Reading as character and parsing the date ourselves reproduces
# base R / pandas exactly. Keep this even though the window now stops in Oct
# 2025 -- the mangling happens at read time, before the filter.
raw <- read_csv(MATERNAL, show_col_types = FALSE,
                col_types = cols(.default = col_character()))

df <- raw %>%
  mutate(
    adm_date = as.Date(dateadmission),
    live_birth = neotreeoutcome %in% c("LB", "ENND"),
    eligible = resus == "N",
    dcc_y = dcc == "Y", s2s_y = s2s == "Y",
    both_y = (dcc == "Y") & (s2s == "Y"),
    bw = suppressWarnings(as.numeric(bwtdis)),
    bw = if_else(!is.na(bw) & (bw < 300 | bw > 7000), NA_real_, bw),
    ga = suppressWarnings(as.numeric(gestation)),
    ga = if_else(!is.na(ga) & (ga < 20 | ga > 44), NA_real_, ga),
    mode_deliv = case_when(
      modedelivery %in% c("1","2","3","6") ~ "Vaginal (incl. instrumental)",
      modedelivery %in% c("4","5")         ~ "Caesarean section",
      TRUE ~ NA_character_),
    bw_band = case_when(is.na(bw) ~ NA_character_, bw < 1000 ~ "<1000 g",
      bw < 1500 ~ "1000-1499 g", bw < 2500 ~ "1500-2499 g",
      bw < 4000 ~ "2500-3999 g", TRUE ~ ">=4000 g"),
    ga_band = case_when(is.na(ga) ~ NA_character_, ga < 28 ~ "<28 weeks",
      ga < 32 ~ "28-31+6 weeks", ga < 37 ~ "32-36+6 weeks", TRUE ~ ">=37 weeks"))

elig <- df %>% filter(facility == "SMCH", live_birth == TRUE, eligible == TRUE,
                      adm_date >= INT_START, adm_date <= INT_END)
cat(sprintf("Eligible SMCH live births (resus = No), %s to %s: %d\n\n",
            INT_START, INT_END, nrow(elig)))

trow <- function(d, characteristic, level) {
  n <- nrow(d)
  f <- function(k) if (n > 0) sprintf("%d (%.1f%%)", k, 100*k/n) else "0"
  tibble(characteristic = characteristic, level = level, eligible_n = n,
         received_dcc = f(sum(d$dcc_y, na.rm = TRUE)),
         received_essc = f(sum(d$s2s_y, na.rm = TRUE)),
         received_both = f(sum(d$both_y, na.rm = TRUE)))
}
blk <- function(var, label, lvls)
  map_dfr(lvls, ~trow(elig %>% filter(.data[[var]] == .x), label, .x))

table1 <- bind_rows(
  trow(elig, "Overall", ""),
  blk("mode_deliv", "Mode of delivery",
      c("Vaginal (incl. instrumental)", "Caesarean section")),
  blk("bw_band", "Birth weight",
      c("<1000 g","1000-1499 g","1500-2499 g","2500-3999 g",">=4000 g")),
  blk("ga_band", "Gestational age",
      c("<28 weeks","28-31+6 weeks","32-36+6 weeks",">=37 weeks")))

print(as.data.frame(table1), row.names = FALSE)
write_csv(table1, file.path(OUTPUT_DIR, "23a_table1_maternal_uptake_12mo.csv"))

# For reference: uptake among KNOWN (Y or N) records, overall.
kn <- function(v) { k <- sum(elig[[v]] %in% c("Y","N")); y <- sum(elig[[v]] == "Y")
  sprintf("%d/%d (%.1f%%)", y, k, 100*y/max(k,1)) }
cat(sprintf("\nFor reference, uptake of KNOWN (Y or N) records: DCC %s; ESSC %s\n",
            kn("dcc"), kn("s2s")))
cat(sprintf("Missing: mode %d, birthweight %d, gestation %d\n",
            sum(is.na(elig$mode_deliv)), sum(is.na(elig$bw_band)), sum(is.na(elig$ga_band))))
cat("\n=== Script 23 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
