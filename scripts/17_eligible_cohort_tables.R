# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# JULY 2026 REQUESTS -- Script 17: ELIGIBLE-COHORT CHARACTERISTICS + OUTCOMES
# =============================================================================
# Produces the descriptive tables Felicity and Hannah asked for, on the
# ELIGIBLE inborn NNU cohort, broken down by the four-group exposure:
#   Neither / DCC only / ESSC(S2S) only / Both DCC and ESSC.
#
#   Table 2  -- baseline characteristics by group
#               (mode of delivery, birthweight bands, gestational-age bands)
#   Table 3  -- outcome by group: death before discharge (NND), n (%)
#
# SPEC: population = inborn NNU admissions, Nov 2024-Oct 2025; eligibility EXCLUDE
#   5-minute Apgar < 7 (team decision 21 Jul; resus limb dropped); exposure four
#   groups (ESSC-only reported as counts here, no OR -- see Script 18).
#
# DATA-QUALITY HANDLING (see Neotree_variable_issues_catalogue.md):
#   - birthweight : prefer derived birthweight_g; rescue kg; 300-7000 g filter (B1,E)
#   - modedelivery: collapsed to Vaginal/Caesarean, handling BOTH numeric (1-7)
#                   and any legacy text codes defensively (A2,C1)
#   - gender      : normalise legacy Male/Female variants to M/F (A1)
#   - delivinter  : multi-select substring match, robust to legacy coding
#
# OUTCOME NOTE: death before discharge = neotreeoutcome == "NND", computed on
#   direct_match records only (a discharge outcome is required). The full
#   outcome distribution is printed so the team can see the denominator.
#
# DSH note: ASCII-only.
# =============================================================================

library(tidyverse)
library(lubridate)

# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
SCRIPT_DIR <- if (length(script_flag) > 0)
  dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE)) else getwd()

# Data file: 25-May-2026 master (latest cleaning pipeline). See Script 16 header
# note. Switch to ...20260401.csv to cross-check against the previous file.
DATA_FILE <- "ZIM_db_master_joined_to_20260525.csv"
DATA_PATH <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA", DATA_FILE), mustWork = FALSE)
# Outputs tagged by data-file date (see Script 18) so a cross-check run against
# the older master lands in its own subfolder.
OUTPUT_TAG <- sub("^.*_(\\d{8})\\.csv$", "\\1", DATA_FILE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", paste0(OUTPUT_TAG, "_apgar1"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "17_eligible_cohort_tables_log.txt")
log_con  <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")

# ELIGIBILITY -- keep in step with Script 18. Team decision (21 Jul, Marcia+Rachel):
# exclude 5-minute Apgar < 7; drop the resus limb. See Script 18 header.
APGAR_FIELD    <- "apgar1"   # AUG 2026: 1-minute Apgar (Rachel 28 Jul)
APGAR_CUTOFF   <- 7
EXCLUSION_MODE <- "apgar_only"   # "apgar_only" | "apgar_resustype" | "apgar_both"

cat("=============================================================\n")
cat("  Script 17: Eligible-cohort characteristics + outcomes\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD + DERIVE (shared derivation, mirrors Scripts 16/18)
# -----------------------------------------------------------------------------
# DATE-PARSING FIX (Aug 2026): read_csv's type guesser converts
# "datetimeadmission" to POSIXct; ymd_hms() then re-coerces it via as.character(),
# which for a midnight timestamp yields a date-only string ("2025-08-08") that
# ymd_hms cannot parse -> NA -> the record is silently dropped by the window
# filter. Four SMCH inborn admissions were lost this way in the July run
# (4,856 instead of 4,860). Read the column as character and take the date part.
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000,
                 col_types = cols(datetimeadmission = col_character()))
bw_src <- if ("birthweight_g" %in% names(raw)) raw$birthweight_g else raw$birthweight

collapse_mode <- function(x) {
  u <- toupper(trimws(as.character(x)))
  case_when(
    u %in% c("4","5","ECS","ELCS","EMCS","CSPRLAB","CSPOLAB") | grepl("AESAR|CESAR", u) ~ "Caesarean section",
    u %in% c("1","2","3","6","7","SVD","IVD","VENT","FOR","ID","VACUUM","FORCEPS","BREECH") |
      grepl("VAGIN", u) ~ "Vaginal (incl. instrumental)",
    is.na(x) | u == "" ~ NA_character_,
    TRUE ~ "Other/Unknown")
}

df <- raw %>%
  mutate(
    adm_date = as.Date(substr(datetimeadmission, 1, 10)),
    inborn = case_when(
      inorout %in% c("Yes","true","True","In")  ~ TRUE,
      inorout %in% c("No","false","False","Out") ~ FALSE, TRUE ~ NA),
    is_readmission = readmission == "Y",
    matched = match_type == "direct_match",

    birthweight_num = suppressWarnings(as.numeric(bw_src)),
    birthweight_num = if_else(!is.na(birthweight_num) & birthweight_num > 0 &
                                birthweight_num <= 20, birthweight_num*1000, birthweight_num),
    birthweight_num = if_else(!is.na(birthweight_num) &
                                (birthweight_num < 300 | birthweight_num > 7000), NA_real_, birthweight_num),
    gestation_num = suppressWarnings(as.numeric(gestation)),
    gestation_num = if_else(!is.na(gestation_num) &
                              (gestation_num < 20 | gestation_num > 44), NA_real_, gestation_num),
    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    apgar1_num = if_else(!is.na(apgar1_num) & (apgar1_num < 0 | apgar1_num > 10), NA_real_, apgar1_num),
    apgar5_num = suppressWarnings(as.numeric(apgar5)),
    apgar5_num = if_else(!is.na(apgar5_num) & (apgar5_num < 0 | apgar5_num > 10), NA_real_, apgar5_num),
    apgar_excl = if (APGAR_FIELD == "apgar5") apgar5_num else apgar1_num,

    sex = case_when(
      toupper(trimws(as.character(gender))) %in% c("F","FEMALE") ~ "F",
      toupper(trimws(as.character(gender))) %in% c("M","MALE")   ~ "M", TRUE ~ NA_character_),

    mode_deliv = collapse_mode(modedelivery),

    bw_band = case_when(
      is.na(birthweight_num)        ~ NA_character_,
      birthweight_num < 1000        ~ "<1000 g",
      birthweight_num < 1500        ~ "1000-1499 g",
      birthweight_num < 2500        ~ "1500-2499 g",
      birthweight_num < 4000        ~ "2500-3999 g",
      TRUE                          ~ ">=4000 g"),
    ga_band = case_when(
      is.na(gestation_num)  ~ NA_character_,
      gestation_num < 28    ~ "<28 weeks",
      gestation_num < 32    ~ "28-31+6 weeks",
      gestation_num < 37    ~ "32-36+6 weeks",
      TRUE                  ~ ">=37 weeks"),

    # BVM/CPR ascertainment from both fields (resustype is post-Mar-2025 only)
    rt_bvmcpr = !is.na(resustype) & grepl("BVM|CPR", resustype),
    rs_bvmcpr = !is.na(resus)     & grepl("BVM|CPR", resus),
    excl_apgar = !is.na(apgar_excl) & apgar_excl < APGAR_CUTOFF,
    excl_resus = if (EXCLUSION_MODE == "apgar_only") FALSE
                 else if (EXCLUSION_MODE == "apgar_both") (rt_bvmcpr | rs_bvmcpr)
                 else rt_bvmcpr,
    excluded = excl_apgar | excl_resus,

    # four-group exposure (known delivinter only)
    di_known = !(is.na(delivinter) | trimws(delivinter) == "" | trimws(delivinter) == "U"),
    has_dcc  = di_known & grepl("DCC", delivinter),
    has_s2s  = di_known & grepl("S2S", delivinter),
    exposure4 = factor(case_when(
      !di_known          ~ NA_character_,
      has_dcc & has_s2s  ~ "Both DCC and ESSC",
      has_dcc & !has_s2s ~ "DCC only",
      !has_dcc & has_s2s ~ "ESSC only",
      TRUE               ~ "Neither"),
      levels = c("Neither","DCC only","ESSC only","Both DCC and ESSC")),

    nnd = neotreeoutcome == "NND")

# Eligible cohort with known exposure
elig <- df %>%
  filter(facility == "SMCH", inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission),
         adm_date >= INT_START, adm_date <= INT_END,
         excluded == FALSE, di_known == TRUE)

cat(sprintf("Eligible cohort with known exposure: %d\n", nrow(elig)))
cat("Exposure group sizes:\n")
print(as.data.frame(elig %>% count(exposure4)), row.names = FALSE)
cat("\nFull outcome distribution (eligible cohort):\n")
print(as.data.frame(elig %>% count(neotreeoutcome, sort = TRUE)), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. TABLE 2 -- CHARACTERISTICS BY GROUP  (n and column %)
# -----------------------------------------------------------------------------
# Helper: for a categorical variable, n (col%) within each exposure group + Overall.
char_block <- function(data, var, var_label) {
  grps <- c("Overall", levels(data$exposure4))
  lvls <- data %>% filter(!is.na(.data[[var]])) %>% pull(.data[[var]]) %>% unique()
  # keep a sensible order for banded variables
  ord <- c("Vaginal (incl. instrumental)","Caesarean section","Other/Unknown",
           "<1000 g","1000-1499 g","1500-2499 g","2500-3999 g",">=4000 g",
           "<28 weeks","28-31+6 weeks","32-36+6 weeks",">=37 weeks")
  lvls <- c(intersect(ord, lvls), setdiff(lvls, ord))
  out <- list()
  for (lv in lvls) {
    row <- tibble(characteristic = var_label, level = lv)
    for (g in grps) {
      d <- if (g == "Overall") data else data %>% filter(exposure4 == g)
      denom <- sum(!is.na(d[[var]]))
      n <- sum(d[[var]] == lv, na.rm = TRUE)
      row[[g]] <- sprintf("%d (%.1f%%)", n, ifelse(denom > 0, 100*n/denom, 0))
    }
    out[[lv]] <- row
  }
  bind_rows(out)
}

# header row with group Ns
n_row <- tibble(characteristic = "N (eligible, known exposure)", level = "")
n_row[["Overall"]] <- as.character(nrow(elig))
for (g in levels(elig$exposure4)) n_row[[g]] <- as.character(sum(elig$exposure4 == g))

table2 <- bind_rows(
  n_row,
  char_block(elig, "mode_deliv", "Mode of delivery"),
  char_block(elig, "bw_band",    "Birth weight"),
  char_block(elig, "ga_band",    "Gestational age"))

cat("\n================= TABLE 2: Characteristics by group =================\n")
print(as.data.frame(table2), row.names = FALSE)
write_csv(table2, file.path(OUTPUT_DIR, "17a_table2_characteristics.csv"))

# missingness note for the three characteristics
miss <- tibble(
  variable = c("mode_deliv","bw_band","ga_band"),
  n_missing = c(sum(is.na(elig$mode_deliv)), sum(is.na(elig$bw_band)), sum(is.na(elig$ga_band))),
  pct_missing = round(100*c(sum(is.na(elig$mode_deliv)), sum(is.na(elig$bw_band)),
                            sum(is.na(elig$ga_band)))/nrow(elig),1))
cat("\nCharacteristic missingness (eligible cohort):\n")
print(as.data.frame(miss), row.names = FALSE)
write_csv(miss, file.path(OUTPUT_DIR, "17b_table2_missingness.csv"))

# -----------------------------------------------------------------------------
# 3. TABLE 3 -- OUTCOME BY GROUP (death before discharge, matched records)
# -----------------------------------------------------------------------------
elig_m <- elig %>% filter(matched == TRUE)
cat(sprintf("\nEligible cohort with a discharge outcome (direct_match): %d\n", nrow(elig_m)))

by_grp <- elig_m %>%
  group_by(group = as.character(exposure4)) %>%
  summarise(total_infants = n(), deaths = sum(nnd, na.rm = TRUE), .groups = "drop")
overall <- tibble(group = "Overall", total_infants = nrow(elig_m),
                  deaths = sum(elig_m$nnd, na.rm = TRUE))
table3 <- bind_rows(overall, by_grp) %>%
  mutate(group = factor(group, levels = c("Overall", levels(elig_m$exposure4)))) %>%
  arrange(group) %>%
  mutate(death_pct = ifelse(total_infants > 0, round(100*deaths/total_infants, 1), NA_real_),
         death_n_pct = sprintf("%d (%.1f%%)", deaths, death_pct))

cat("\n================= TABLE 3: Death before discharge by group =================\n")
print(as.data.frame(table3), row.names = FALSE)
write_csv(table3, file.path(OUTPUT_DIR, "17c_table3_outcomes.csv"))

cat("\nSaved: 17a_table2_characteristics.csv, 17b_table2_missingness.csv,\n")
cat("       17c_table3_outcomes.csv\n")
cat("\n=== Script 17 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
