# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# JULY 2026 REQUESTS -- Script 16: DATA READINESS CHECKS
# =============================================================================
# Purpose (plan section 5):
#   Before running the new eligible-cohort outcome analysis, confirm that the
#   variables the new specification depends on are actually present and usable
#   in the master NeoTree file. Specifically:
#     - gender / sex           (new adjustment covariate, Michelle's request)
#     - resus / resustype      (needed for the eligibility exclusion BVM / CPR)
#     - apgar1                 (needed for the eligibility exclusion apgar1 < 6)
#     - delivinter             (four-group exposure: Neither/DCC/ESSC/Both)
#     - birthweight, gender    (adjusted-model covariates)
#     - neotreeoutcome         (outcome: death before discharge = NND)
#
#   This script does NOT model anything. It prints and saves completeness and
#   value distributions so the team can confirm the eligibility definition is
#   runnable (and pick the resus fallback if resustype is too incomplete).
#
# KEY DEFINITIONS BEING TESTED (from Rachel's 14 July 15:20 email):
#   Population  : all inborn neonates admitted to SMCH NNU, Nov 2024 - Oct 2025.
#   Eligibility : EXCLUDE apgar1 < 6 OR resustype contains BVM or CPR.
#   Exposure    : Neither / DCC only / ESSC(S2S) only / Both, from delivinter.
#
# CONFIRMED VARIABLE CODES (from dictionary_zim_admissions.xlsx ValueMaps):
#   gender     : F, M, U (Unsure), AmbG (Ambiguous genitalia)
#   resus      : Y / N   (were any resuscitation efforts made)
#   resustype  : Stim, BVM, CPR, RESUSMED, Suc, O2   (multi-select)
#   delivinter : DCC, S2S, BF, Norm (=NONE), U (=Unknown)   (multi-select)
#
# DSH note: ASCII-only.
# =============================================================================

library(tidyverse)
library(lubridate)

# -----------------------------------------------------------------------------
# 0. SETUP -- locate own directory; data one level up in 00-DATA
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

# Data file: the 25-May-2026 master (latest cleaning pipeline). Verified against
# the Apr file: the intervention cohort is identical (4,860 babies; matched and
# deaths differ by 1) but it recovers ~1,000 more babies with KNOWN DCC/ESSC
# exposure (1,866 -> 2,868), which directly strengthens the exposure analysis.
# To cross-check against the previous file, set this to ...20260401.csv.
DATA_FILE <- "ZIM_db_master_joined_to_20260525.csv"
DATA_PATH <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA", DATA_FILE), mustWork = FALSE)
# Outputs tagged by data-file date (see Script 18) so a cross-check run against
# the older master lands in its own subfolder.
OUTPUT_TAG <- sub("^.*_(\\d{8})\\.csv$", "\\1", DATA_FILE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", paste0(OUTPUT_TAG, "_apgar1"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUTPUT_DIR, "16_data_readiness_checks_log.txt")
log_con  <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

# Analysis window -- QI paper period (Marcia: keep same period as QI).
# Change these two dates if the team confirms a different window.
INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")

# ELIGIBILITY (team decision, 21 Jul 2026 -- Marcia + Rachel):
#   Exclude 5-minute Apgar < 7 (apgar5); DROP the resus/resustype limb because the
#   resus fields were incompletely filled. This supersedes the earlier
#   "apgar1 < 6 OR resustype BVM/CPR" draft. Keep in step with Scripts 17/18.
APGAR_FIELD    <- "apgar1"        # AUG 2026: 1-minute Apgar (Rachel 28 Jul, confirmed with Marcia)
APGAR_CUTOFF   <- 7               # exclude if Apgar < this (apgar1 < 7 = 0-6)
EXCLUSION_MODE <- "apgar_only"    # "apgar_only" | "apgar_resustype" | "apgar_both"

cat("=============================================================\n")
cat("  Script 16: DATA READINESS CHECKS (July requests)\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat(sprintf("  Data   : %s\n", DATA_PATH))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD + DERIVE POPULATION
# -----------------------------------------------------------------------------
# DATE-PARSING FIX (Aug 2026): read_csv's type guesser converts
# "datetimeadmission" to POSIXct; ymd_hms() then re-coerces it via as.character(),
# which for a midnight timestamp yields a date-only string ("2025-08-08") that
# ymd_hms cannot parse -> NA -> the record is silently dropped by the window
# filter. Four SMCH inborn admissions were lost this way in the July run
# (4,856 instead of 4,860). Read the column as character and take the date part.
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000,
                 col_types = cols(datetimeadmission = col_character()))
cat(sprintf("Rows loaded: %d\n", nrow(raw)))

needed <- c("facility","inorout","readmission","match_type","datetimeadmission",
            "delivinter","birthweight","gestation","apgar1","apgar5","temperature",
            "gender","resus","resustype","neotreeoutcome","modedelivery")
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  cat("\n*** WARNING: expected columns not found in the CSV: ***\n")
  cat("   ", paste(missing_cols, collapse = ", "), "\n")
  cat("    (Check the master file version. Everything else will still run.)\n\n")
}

# Prefer the derived canonical grams column birthweight_g if present (see
# variable-issues catalogue E); else fall back to raw birthweight.
bw_src <- if ("birthweight_g" %in% names(raw)) raw$birthweight_g else raw$birthweight

df <- raw %>%
  mutate(
    adm_date = as.Date(substr(datetimeadmission, 1, 10)),
    inborn = case_when(
      inorout %in% c("Yes","true","True","In")   ~ TRUE,
      inorout %in% c("No","false","False","Out")  ~ FALSE,
      TRUE ~ NA),
    is_readmission = readmission == "Y",
    birthweight_num = suppressWarnings(as.numeric(bw_src)),
    # rescue residual kg (<=20 -> *1000; catalogue B1), then 300-7000 g filter
    birthweight_num = if_else(!is.na(birthweight_num) & birthweight_num > 0 &
                                birthweight_num <= 20, birthweight_num * 1000, birthweight_num),
    birthweight_num = if_else(!is.na(birthweight_num) &
                                (birthweight_num < 300 | birthweight_num > 7000),
                              NA_real_, birthweight_num),
    apgar1_num      = suppressWarnings(as.numeric(apgar1)),
    apgar1_num      = if_else(!is.na(apgar1_num) & (apgar1_num < 0 | apgar1_num > 10),
                              NA_real_, apgar1_num),
    apgar5_num      = suppressWarnings(as.numeric(apgar5)),
    apgar5_num      = if_else(!is.na(apgar5_num) & (apgar5_num < 0 | apgar5_num > 10),
                              NA_real_, apgar5_num),
    apgar_excl      = if (APGAR_FIELD == "apgar5") apgar5_num else apgar1_num,
    # sex normalised to M/F (legacy Male/Female variants; U/AmbG -> NA)
    sex = case_when(
      toupper(trimws(as.character(gender))) %in% c("F","FEMALE") ~ "F",
      toupper(trimws(as.character(gender))) %in% c("M","MALE")   ~ "M",
      TRUE ~ NA_character_))

pop <- df %>%
  filter(facility == "SMCH", inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission),
         adm_date >= INT_START, adm_date <= INT_END)

N <- nrow(pop)
cat(sprintf("\nSMCH inborn non-readmission, %s to %s : %d\n", INT_START, INT_END, N))
cat(sprintf("  direct_match (has discharge/outcome)     : %d\n",
            sum(pop$match_type == "direct_match", na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 2. COMPLETENESS OF KEY VARIABLES
# -----------------------------------------------------------------------------
cat("\n-------------------------------------------------------------\n")
cat("  COMPLETENESS in the population (n non-missing / %)\n")
cat("-------------------------------------------------------------\n")
check_vars <- c("delivinter","birthweight","gestation","apgar1","apgar5","temperature",
                "gender","resus","resustype","neotreeoutcome","modedelivery")
completeness <- map_dfr(check_vars, function(v) {
  x <- pop[[v]]
  nn <- sum(!is.na(x) & trimws(as.character(x)) != "")
  tibble(variable = v, n_present = nn, pct_present = round(100 * nn / N, 1))
})
print(as.data.frame(completeness), row.names = FALSE)
write_csv(completeness, file.path(OUTPUT_DIR, "16a_completeness.csv"))

# -----------------------------------------------------------------------------
# 3. VALUE DISTRIBUTIONS FOR THE NEW VARIABLES
# -----------------------------------------------------------------------------
show_counts <- function(v, label) {
  cat(sprintf("\n%s (%s):\n", label, v))
  tab <- pop %>% count(value = .data[[v]], sort = TRUE)
  print(as.data.frame(head(tab, 20)), row.names = FALSE)
}
show_counts("gender", "Sex (raw; legacy variants revealed here)")
cat(sprintf("\nSex normalised to M/F (usable for adjustment): %d / %d (%.1f%%)\n",
            sum(!is.na(pop$sex)), N, 100*sum(!is.na(pop$sex))/N))
print(as.data.frame(pop %>% count(sex)), row.names = FALSE)
show_counts("resus", "Resuscitation Y/N")
show_counts("resustype", "Resuscitation type (multi-select, raw)")
show_counts("delivinter", "Delivery interventions (multi-select, raw)")
show_counts("modedelivery", "Mode of delivery (raw)")

# BVM/CPR ascertainment: resustype (post-Mar-2025 field) vs resus (100% complete,
# carries the same multi-select tokens) vs either. Also split by era, because the
# resustype field is structurally absent before ~Mar 2025 -- so the resustype-only
# exclusion is applied at different strength across the study window.
pop <- pop %>%
  mutate(
    rt_bvmcpr = !is.na(resustype) & grepl("BVM|CPR", resustype),
    rs_bvmcpr = !is.na(resus)     & grepl("BVM|CPR", resus),
    either_bvmcpr = rt_bvmcpr | rs_bvmcpr,
    era = if_else(adm_date < as.Date("2025-03-01"), "pre-Mar2025", "Mar2025-on"))
cat("\n-------------------------------------------------------------\n")
cat("  BVM/CPR ASCERTAINMENT: resustype vs resus, by era\n")
cat("-------------------------------------------------------------\n")
asc <- pop %>% group_by(era) %>%
  summarise(n = n(),
            via_resustype = sum(rt_bvmcpr),
            via_resus     = sum(rs_bvmcpr),
            via_either    = sum(either_bvmcpr), .groups = "drop")
print(as.data.frame(asc), row.names = FALSE)
cat(sprintf("\n  Whole window: resustype %d, resus %d, either %d BVM/CPR babies.\n",
            sum(pop$rt_bvmcpr), sum(pop$rs_bvmcpr), sum(pop$either_bvmcpr)))
cat("  Context only: the team dropped the resus/resustype limb (this field was\n")
cat("  incompletely filled and structurally absent pre-Mar2025). Eligibility is\n")
cat("  now Apgar-based -- see section 4.\n")
write_csv(asc, file.path(OUTPUT_DIR, "16c_bvmcpr_ascertainment_by_era.csv"))

# -----------------------------------------------------------------------------
# 4. ELIGIBILITY EXCLUSION -- final rule: exclude Apgar (apgar5) < cutoff
# -----------------------------------------------------------------------------
# Team decision 21 Jul: exclude 5-minute Apgar < 7; resus/resustype limb dropped.
# "Unassessable" = the chosen Apgar field is missing -> eligibility not confirmable.
pop <- pop %>%
  mutate(
    excl_apgar5 = !is.na(apgar5_num) & apgar5_num < 7,   # July rule, reference only
    excl_apgar1 = !is.na(apgar1_num) & apgar1_num < 7,   # AUG 2026 selected rule
    excl_selected = !is.na(apgar_excl) & apgar_excl < APGAR_CUTOFF,
    excluded    = excl_selected,
    unassessable = is.na(apgar_excl),
    eligible    = !excluded)

cat("\n-------------------------------------------------------------\n")
cat(sprintf("  ELIGIBILITY -- FINAL RULE: exclude %s < %d\n", APGAR_FIELD, APGAR_CUTOFF))
cat("-------------------------------------------------------------\n")
cat(sprintf("  Population                                   : %d\n", N))
cat(sprintf("  Excluded by apgar1 < 7 (SELECTED, Aug 2026)  : %d\n", sum(pop$excl_apgar1)))
cat(sprintf("  Excluded by apgar5 < 7 (July rule, reference): %d\n", sum(pop$excl_apgar5)))
cat(sprintf("  ELIGIBLE cohort (%s < %d)                : %d\n",
            APGAR_FIELD, APGAR_CUTOFF, sum(pop$eligible)))
cat(sprintf("  Of eligible, %s missing (eligibility uncertain): %d\n",
            APGAR_FIELD, sum(pop$eligible & pop$unassessable)))

# -----------------------------------------------------------------------------
# 5. FOUR-GROUP EXPOSURE AMONG ELIGIBLE (known delivinter only)
# -----------------------------------------------------------------------------
elig <- pop %>% filter(eligible)
elig <- elig %>%
  mutate(
    di_known = !(is.na(delivinter) | trimws(delivinter) == "" | trimws(delivinter) == "U"),
    has_dcc  = di_known & grepl("DCC", delivinter),
    has_s2s  = di_known & grepl("S2S", delivinter),
    exposure4 = case_when(
      !di_known           ~ NA_character_,
      has_dcc & has_s2s   ~ "Both DCC and ESSC",
      has_dcc & !has_s2s  ~ "DCC only",
      !has_dcc & has_s2s  ~ "ESSC only",
      TRUE                ~ "Neither"))
cat("\n-------------------------------------------------------------\n")
cat("  FOUR-GROUP EXPOSURE among ELIGIBLE (delivinter known)\n")
cat("-------------------------------------------------------------\n")
cat(sprintf("  Eligible with delivinter KNOWN          : %d / %d (%.1f%%)\n",
            sum(elig$di_known), nrow(elig), 100*sum(elig$di_known)/nrow(elig)))
print(as.data.frame(elig %>% filter(di_known) %>% count(exposure4)), row.names = FALSE)

# outcome availability among eligible-known
elig_known <- elig %>% filter(di_known)
cat(sprintf("\n  Eligible-known with direct_match (outcome available): %d\n",
            sum(elig_known$match_type == "direct_match", na.rm = TRUE)))
nnd_ok <- elig_known %>% filter(match_type == "direct_match")
cat(sprintf("  Deaths before discharge (NND) among those : %d\n",
            sum(nnd_ok$neotreeoutcome == "NND", na.rm = TRUE)))

readiness <- tibble(
  metric = c("population","excluded_apgar5_lt7_ref","excluded_apgar1_lt7_selected",
             "eligible_selected","eligible_unassessable",
             "eligible_delivinter_known","eligible_known_matched",
             "nnd_events_in_analysis_set"),
  value  = c(N, sum(pop$excl_apgar5), sum(pop$excl_apgar1),
             sum(pop$eligible), sum(pop$eligible & pop$unassessable),
             sum(elig$di_known), sum(nnd_ok$match_type == "direct_match", na.rm=TRUE),
             sum(nnd_ok$neotreeoutcome == "NND", na.rm = TRUE)))
write_csv(readiness, file.path(OUTPUT_DIR, "16b_readiness_summary.csv"))

cat("\nSaved: 16a_completeness.csv, 16b_readiness_summary.csv, 16c_...csv\n")
cat("\n=== Script 16 complete ===\n")
cat("REVIEW BEFORE MODELLING (final rule = exclude apgar5 < 7, resus limb dropped):\n")
cat(sprintf("  - Is %s complete enough? (see 16a_completeness.csv; ~90%% expected)\n", APGAR_FIELD))
cat("  - Is gender complete enough to adjust for? (see 16a_completeness.csv)\n")
cat("  - Is the delivinter-known denominator acceptable for Tables 2-4?\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
