# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 01: Descriptive Overview of the Dataset
# =============================================================================
# Purpose:
#   Produce a comprehensive descriptive summary of the master dataset as it
#   applies to the golden hour analysis. Covers date ranges, completeness of
#   all key variables, distributions, and flags important data gaps (notably
#   the absence of DCC/ESSC variables at SMCH).
#
# Data source:
#   ZIM_db_master_joined_to_20260525.csv
#   (see brief.md -> Data and infrastructure -> Primary data file)
#
# Filters applied in this script:
#   - Facility: SMCH only
#   - Inborn only (harmonised across script versions)
#   - Readmissions excluded
#
# Population labels used throughout:
#   all_smch     : all SMCH records (no date filter, no inborn filter)
#   smch_inborn  : SMCH inborn non-readmission records (full history)
#   int_pop      : intervention period (Nov 2024 - Oct 2025), inborn
#   hist_pop     : historical period (Jan 2022 - Oct 2024), inborn
#   q10_pop      : int_pop + hist_pop combined (for Q10 temperature trends)
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)

# -- Locate this script's directory -------------------------------------------
# Works when run via: Rscript 01_descriptive_overview.R
# Falls back to getwd() when sourced interactively.
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

# -- Path configuration -------------------------------------------------------
# Paths are relative to the script location (Golden-hour/02-CODE/).
DATA_PATH <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA",
            "ZIM_db_master_joined_to_20260525.csv"),
  mustWork = FALSE
)

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -- Redirect all output to a log file ----------------------------------------
# Uses an explicit file connection so the log is written correctly whether the
# script is run via Rscript or via source() in RStudio.
# split = TRUE sends output to both the file and the console simultaneously.
LOG_FILE <- file.path(OUTPUT_DIR, "01_descriptive_overview_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")
# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")   # intervention period start
INT_END    <- as.Date("2025-10-31")   # intervention period end (inclusive)
HIST_START <- as.Date("2022-01-01")   # historical period start (for Q10)
HIST_END   <- as.Date("2024-10-31")   # historical period end


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
# guess_max = 100000: sample up to 100k rows for type-guessing so that
# mixed-type columns whose problem values fall after row 1000 are handled
# correctly rather than silently coerced to NA with a parse warning.
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)

# Report any residual parse problems so they are visible in the log
parse_probs <- problems(raw)
if (nrow(parse_probs) > 0) {
  cat(sprintf("  WARNING: %d parse problem(s) remain after guess_max fix:\n",
              nrow(parse_probs)))
  print(head(parse_probs, 20))
} else {
  cat("  No parse problems detected.\n")
}

cat(sprintf("  Rows loaded      : %d\n", nrow(raw)))
cat(sprintf("  Columns          : %d\n", ncol(raw)))
cat(sprintf("  Facilities       : %s\n",
            paste(sort(unique(raw$facility)), collapse = ", ")))


# -----------------------------------------------------------------------------
# 2. INITIAL CLEANING AND DERIVED VARIABLES
# -----------------------------------------------------------------------------

df <- raw %>%
  mutate(
    # -- Parse admission datetime ---------------------------------------------
    adm_dt = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date = as.Date(adm_dt),
    adm_year  = year(adm_dt),
    adm_month = floor_date(adm_dt, "month"),

    # -- Harmonise inorout coding across script versions ----------------------
    # Values observed: "Yes"/"true"/"True"/"In" -> inborn
    #                  "No"/"false"/"False"/"Out" -> outborn
    #                  "" / NA -> unknown
    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")      ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out")   ~ FALSE,
      TRUE                                               ~ NA
    ),

    # -- Numeric conversions for key clinical variables -----------------------
    temperature_num  = suppressWarnings(as.numeric(temperature)),
    apgar1_num       = suppressWarnings(as.numeric(apgar1)),
    apgar5_num       = suppressWarnings(as.numeric(apgar5)),
    apgar10_num      = suppressWarnings(as.numeric(apgar10)),
    birthweight_num  = suppressWarnings(as.numeric(birthweight)),
    gestation_num    = suppressWarnings(as.numeric(gestation)),
    hr_num           = suppressWarnings(as.numeric(hr)),
    rr_num           = suppressWarnings(as.numeric(rr)),
    satsair_num      = suppressWarnings(as.numeric(satsair)),
    age_raw          = suppressWarnings(as.numeric(age)),   # age in hours (raw)
    # Cap age at 28 days (672 hours). Values above this are data entry errors
    # (confirmed: the pipeline stores age in hours throughout; the outliers are
    # not days-stored-as-hours -- they are erroneous values that should be NA).
    # See Neotree/brief.md -> Known variable issues -> age.
    age_num          = if_else(age_raw > 672, NA_real_, age_raw),
    bloodsugarmmol_num = suppressWarnings(as.numeric(bloodsugarmmol)),

    # -- Resuscitation: any non-none intervention ------------------------------
    resus_any = case_when(
      resus %in% c("NONE", "Norm")  ~ FALSE,
      is.na(resus) | resus == ""    ~ NA,
      TRUE                          ~ TRUE
    ),

    # -- DCC eligibility proxy: inborn + (no resus required OR apgar1 > 6) ---
    # Note: DCC/ESSC direct variables (corddelay, ikmc) are NOT present at SMCH.
    dcc_eligible_proxy = case_when(
      is.na(apgar1_num) & is.na(resus_any) ~ NA,
      !is.na(apgar1_num) & apgar1_num > 6  ~ TRUE,
      resus_any == FALSE                    ~ TRUE,
      TRUE                                  ~ FALSE
    ),

    # -- Normal temperature at admission (36.5-37.5 C, WHO definition) --------
    temp_normal = case_when(
      is.na(temperature_num)                               ~ NA,
      temperature_num >= 36.5 & temperature_num <= 37.5   ~ TRUE,
      TRUE                                                 ~ FALSE
    ),

    # -- Hypothermia categories (WHO) -----------------------------------------
    temp_cat = case_when(
      is.na(temperature_num)              ~ NA_character_,
      temperature_num < 32.0             ~ "Severe hypothermia (<32)",
      temperature_num < 36.0             ~ "Moderate hypothermia (32-35.9)",
      temperature_num < 36.5             ~ "Mild hypothermia (36-36.4)",
      temperature_num <= 37.5            ~ "Normal (36.5-37.5)",
      temperature_num <= 38.5            ~ "Mild hyperthermia (37.6-38.5)",
      TRUE                               ~ "Hyperthermia (>38.5)"
    ),

    # -- Study period flag ----------------------------------------------------
    period = case_when(
      adm_date >= INT_START  & adm_date <= INT_END   ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END  ~ "Historical",
      TRUE                                            ~ "Other"
    ),

    # -- Matched status -------------------------------------------------------
    matched = match_type == "direct_match",

    # -- Readmission flag (already in data, coerce to logical) ----------------
    is_readmission = readmission == "Y",

    # -- Delivery interventions (delivinter): DCC, S2S, BF -------------------
    # Introduced Oct 2024; zero non-missing values in the historical period.
    # Multi-select stored in two encoding formats: {DCC,S2S} or ["DCC","S2S"].
    # grepl() handles both. "U" = Unknown (excluded from uptake denominators).
    dcc_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                     ~ NA,
      grepl("DCC", delivinter, fixed = TRUE)           ~ TRUE,
      TRUE                                             ~ FALSE
    ),
    s2s_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                     ~ NA,
      grepl("S2S", delivinter, fixed = TRUE)           ~ TRUE,
      TRUE                                             ~ FALSE
    ),
    bf_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                     ~ NA,
      grepl("BF", delivinter, fixed = TRUE)            ~ TRUE,
      TRUE                                             ~ FALSE
    ),
    no_intervention = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                     ~ NA,
      grepl("DCC|S2S|BF", delivinter)                  ~ FALSE,
      TRUE                                             ~ TRUE
    ),
    any_gh_intervention = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                     ~ NA,
      grepl("DCC|S2S", delivinter)                     ~ TRUE,
      TRUE                                             ~ FALSE
    ),
    # TRUE = delivinter recorded and not Unknown -- use as uptake denominator
    delivinter_known = (
      !is.na(delivinter) &
      trimws(delivinter) != "" &
      trimws(delivinter) != "U"
    )
  )

cat("\nDerived variables created.\n")


# -----------------------------------------------------------------------------
# 3. DEFINE ANALYSIS POPULATIONS
# -----------------------------------------------------------------------------

all_smch <- df %>% filter(facility == "SMCH")

smch_inborn <- all_smch %>%
  filter(inborn == TRUE, is_readmission == FALSE | is.na(is_readmission))

int_pop  <- smch_inborn %>% filter(period == "Intervention")
hist_pop <- smch_inborn %>% filter(period == "Historical")
q10_pop  <- smch_inborn %>% filter(period %in% c("Intervention", "Historical"))

cat("\n=== POPULATION SIZES ===\n")
cat(sprintf("  Full master dataset         : %d records\n", nrow(df)))
cat(sprintf("  All SMCH records            : %d records\n", nrow(all_smch)))
cat(sprintf("  SMCH inborn, no readmissions: %d records\n", nrow(smch_inborn)))
cat(sprintf("  Intervention period (inborn): %d records\n", nrow(int_pop)))
cat(sprintf("    of which matched          : %d records (%.1f%%)\n",
            sum(int_pop$matched, na.rm = TRUE),
            100 * mean(int_pop$matched, na.rm = TRUE)))
cat(sprintf("  Historical period (inborn)  : %d records\n", nrow(hist_pop)))
cat(sprintf("  Q10 combined (inborn)       : %d records\n", nrow(q10_pop)))


# -----------------------------------------------------------------------------
# 4. DATE RANGE OVERVIEW
# -----------------------------------------------------------------------------

cat("\n=== DATE RANGES (SMCH inborn) ===\n")
cat(sprintf("  Earliest admission: %s\n", min(smch_inborn$adm_date, na.rm = TRUE)))
cat(sprintf("  Latest admission  : %s\n", max(smch_inborn$adm_date, na.rm = TRUE)))
cat(sprintf("  Intervention      : %s to %s\n", INT_START, INT_END))
cat(sprintf("  Historical (Q10)  : %s to %s\n", HIST_START, HIST_END))

# Admissions by year
cat("\nAdmissions by year (SMCH inborn):\n")
smch_inborn %>%
  count(adm_year) %>%
  print(n = Inf)

# Admissions by month (intervention period)
cat("\nAdmissions by month (intervention period, SMCH inborn):\n")
int_pop %>%
  mutate(month = format(adm_month, "%Y-%m")) %>%
  count(month) %>%
  print(n = Inf)


# -----------------------------------------------------------------------------
# 5.  DCC / ESSC UPTAKE AT SMCH (delivinter variable)
# -----------------------------------------------------------------------------
# DCC and ESSC are captured in the multi-select field 'delivinter', which was
# added to the SMCH Neotree form in October 2024. Codes:
#   DCC  = Delayed Cord Clamping (>60 secs)
#   S2S  = Skin to Skin (within 1st hour of life)  [ESSC]
#   BF   = Breastfeeding (within 1st hour of life)
#   NONE/Norm = no Golden Hour intervention received
#   U    = Unknown (excluded from uptake denominators)
#
# Historical period has zero delivinter records (expected -- field not in form).
# Note: corddelay and ikmc remain absent at SMCH; delivinter is the correct key.

cat("\n")
cat("================================================================\n")
cat("  DCC / ESSC UPTAKE AT SMCH (delivinter)\n")
cat("================================================================\n")

# Raw value distribution (intervention period)
cat("\ndelivinter raw values -- intervention period (SMCH inborn):\n")
int_pop %>%
  mutate(dv = ifelse(is.na(delivinter) | delivinter == "",
                     "_MISSING_", delivinter)) %>%
  count(dv, sort = TRUE) %>%
  mutate(pct = round(100 * n / nrow(int_pop), 1)) %>%
  print(n = Inf)

# Recording status summary
n_int_total  <- nrow(int_pop)
n_dv_missing <- sum(is.na(int_pop$delivinter) |
                    trimws(int_pop$delivinter) == "")
n_dv_unknown <- sum(trimws(int_pop$delivinter) == "U", na.rm = TRUE)
n_dv_known   <- sum(int_pop$delivinter_known, na.rm = TRUE)

cat(sprintf("\ndelivinter recording status (intervention period, n=%d):\n",
            n_int_total))
cat(sprintf("  Recorded (any value)    : %d (%.1f%%)\n",
            n_int_total - n_dv_missing,
            100 * (n_int_total - n_dv_missing) / n_int_total))
cat(sprintf("  of which Unknown (U)    : %d (%.1f%% of total)\n",
            n_dv_unknown, 100 * n_dv_unknown / n_int_total))
cat(sprintf("  Missing (not recorded)  : %d (%.1f%%)\n",
            n_dv_missing, 100 * n_dv_missing / n_int_total))
cat(sprintf("  Known (excl. U/missing) : %d (%.1f%% of total)\n",
            n_dv_known, 100 * n_dv_known / n_int_total))

# Uptake rates (denominator = delivinter_known)
cat(sprintf("\n--- Uptake (denominator = known delivinter, n=%d) ---\n",
            n_dv_known))
int_known <- int_pop %>% filter(delivinter_known == TRUE)

cat(sprintf("  DCC received          : %d / %d (%.1f%%)\n",
            sum(int_known$dcc_received, na.rm = TRUE), n_dv_known,
            100 * mean(int_known$dcc_received, na.rm = TRUE)))
cat(sprintf("  S2S received (ESSC)   : %d / %d (%.1f%%)\n",
            sum(int_known$s2s_received, na.rm = TRUE), n_dv_known,
            100 * mean(int_known$s2s_received, na.rm = TRUE)))
cat(sprintf("  BF within 1st hour    : %d / %d (%.1f%%)\n",
            sum(int_known$bf_received, na.rm = TRUE), n_dv_known,
            100 * mean(int_known$bf_received, na.rm = TRUE)))
cat(sprintf("  DCC + S2S (both)      : %d / %d (%.1f%%)\n",
            sum(int_known$dcc_received & int_known$s2s_received, na.rm=TRUE),
            n_dv_known,
            100 * mean(int_known$dcc_received & int_known$s2s_received,
                       na.rm = TRUE)))
cat(sprintf("  Any GH (DCC or S2S)   : %d / %d (%.1f%%)\n",
            sum(int_known$any_gh_intervention, na.rm = TRUE), n_dv_known,
            100 * mean(int_known$any_gh_intervention, na.rm = TRUE)))
cat(sprintf("  No intervention       : %d / %d (%.1f%%)\n",
            sum(int_known$no_intervention, na.rm = TRUE), n_dv_known,
            100 * mean(int_known$no_intervention, na.rm = TRUE)))

# Monthly uptake (Q1-Q2)
cat("\n--- Monthly DCC and S2S uptake (known delivinter only) ---\n")
monthly_uptake <- int_pop %>%
  filter(delivinter_known == TRUE) %>%
  mutate(month = format(adm_month, "%Y-%m")) %>%
  group_by(month) %>%
  summarise(
    n_known  = n(),
    n_dcc    = sum(dcc_received,  na.rm = TRUE),
    n_s2s    = sum(s2s_received,  na.rm = TRUE),
    n_bf     = sum(bf_received,   na.rm = TRUE),
    n_both   = sum(dcc_received & s2s_received, na.rm = TRUE),
    n_none   = sum(no_intervention, na.rm = TRUE),
    pct_dcc  = round(100 * n_dcc / n_known, 1),
    pct_s2s  = round(100 * n_s2s / n_known, 1),
    pct_both = round(100 * n_both / n_known, 1),
    pct_none = round(100 * n_none / n_known, 1),
    .groups  = "drop"
  )

print(monthly_uptake, n = Inf)
write_csv(monthly_uptake,
          file.path(OUTPUT_DIR, "01d_monthly_dcc_s2s_uptake.csv"))
cat("\nSaved: 01d_monthly_dcc_s2s_uptake.csv\n")

# Uptake by birthweight group (Q3)
cat("\n--- DCC/S2S uptake by birthweight group (known delivinter only) ---\n")
bw_uptake <- int_pop %>%
  filter(delivinter_known == TRUE) %>%
  group_by(bwgroup) %>%
  summarise(
    n_known  = n(),
    n_dcc    = sum(dcc_received,  na.rm = TRUE),
    n_s2s    = sum(s2s_received,  na.rm = TRUE),
    n_both   = sum(dcc_received & s2s_received, na.rm = TRUE),
    n_none   = sum(no_intervention, na.rm = TRUE),
    pct_dcc  = round(100 * n_dcc / n_known, 1),
    pct_s2s  = round(100 * n_s2s / n_known, 1),
    pct_both = round(100 * n_both / n_known, 1),
    pct_none = round(100 * n_none / n_known, 1),
    .groups  = "drop"
  ) %>%
  arrange(factor(bwgroup,
                 levels = c("ELBW", "VLBW", "LBW", "NBW", "HBW", "Unknown")))

print(bw_uptake, n = Inf)
write_csv(bw_uptake,
          file.path(OUTPUT_DIR, "01e_bwgroup_dcc_s2s_uptake.csv"))
cat("\nSaved: 01e_bwgroup_dcc_s2s_uptake.csv\n")

cat("================================================================\n")


# -----------------------------------------------------------------------------
# 6. VARIABLE COMPLETENESS TABLE
# -----------------------------------------------------------------------------

# Define the variables of interest for the golden hour analysis
gh_vars <- tribble(
  ~variable,           ~label,
  # -- Timing / age
  "datetimeadmission", "Datetime of admission",
  "age",               "Age at admission (hours, numeric field)",
  "agecategory",       "Age category (text label)",
  # -- Demographics / birth
  "birthweight",       "Birthweight (g)",
  "bwgroup",           "Birthweight group (NBW/LBW/VLBW/ELBW/HBW)",
  "gestation",         "Gestational age (weeks)",
  "gestgroup",         "Gestational age group",
  "gender",            "Sex",
  "modedelivery",      "Mode of delivery",
  "typebirth",         "Type of birth (singleton/twin/triplet)",
  "inorout",           "Inborn / outborn",
  "birthplacesame",    "Born at SMCH",
  "readmission",       "Readmission flag",
  # -- Maternal
  "matageyrs",         "Maternal age (years)",
  "ansteroids",        "Antenatal steroids given",
  "mathivtest",        "Maternal HIV test",
  "prom",              "Prolonged rupture of membranes",
  # -- Golden hour interventions
  "delivinter",        "Delivery interventions (DCC/S2S/BF/None/U) -- int period only",
  "corddelay",         "Delayed cord clamping -- PGH only, absent at SMCH",
  "ikmc",              "Initial KMC / ESSC -- PGH only, absent at SMCH",
  # -- Clinical at admission
  "crybirth",          "Cried at birth",
  "babycrytriage",     "Baby crying at triage",
  "resus",             "Resuscitation at birth",
  "apgar1",            "Apgar score at 1 minute",
  "apgar5",            "Apgar score at 5 minutes",
  "apgar10",           "Apgar score at 10 minutes",
  "temperature",       "Temperature at admission (C)",
  "tempgroup",         "Temperature group",
  "hr",                "Heart rate at admission (bpm)",
  "rr",                "Respiratory rate at admission (bpm)",
  "satsair",           "SpO2 on air (%)",
  # -- Blood glucose
  "bloodsugarmmol",    "Blood glucose at admission (mmol/L)",
  "bsmonyn",           "Blood sugar monitoring done",
  # -- Feeding
  "feedsadm",          "Feeding type at admission",
  # -- Admission context
  "admreason",         "Primary reason for admission",
  "admittedfrom",      "Admitted from (ward/location)",
  # -- Outcome (discharge side)
  "neotreeoutcome",    "Neotree outcome (DC/NND/TRH/etc)",
  "causedeath",        "Cause of death (if NND)",
  "dischtemp",         "Temperature at discharge (C)",
  "datetimedischarge", "Datetime of discharge",
  "match_type",        "Match type (direct_match / unmatched)"
)

completeness_tbl <- gh_vars %>%
  rowwise() %>%
  mutate(
    n_total    = nrow(smch_inborn),
    n_complete = {
      # Coerce to character first so datetime/numeric columns don't error on
      # != "" comparison (read_csv may parse some columns as POSIXct etc.)
      col_chr <- as.character(smch_inborn[[variable]])
      sum(!is.na(col_chr) & col_chr != "" & col_chr != "NA", na.rm = TRUE)
    },
    pct_complete = round(100 * n_complete / n_total, 1),
    pct_missing  = round(100 - pct_complete, 1)
  ) %>%
  ungroup() %>%
  select(variable, label, n_complete, pct_complete, pct_missing)

cat("\n=== VARIABLE COMPLETENESS (SMCH inborn) ===\n")
print(completeness_tbl, n = Inf)

# Save to CSV
write_csv(completeness_tbl,
          file.path(OUTPUT_DIR, "01a_variable_completeness_smch_inborn.csv"))
cat("\nSaved: 01a_variable_completeness_smch_inborn.csv\n")


# -----------------------------------------------------------------------------
# 7. DESCRIPTIVE STATISTICS -- DEMOGRAPHICS AND BIRTH CHARACTERISTICS
# -----------------------------------------------------------------------------

cat("\n=== DESCRIPTIVE STATISTICS: DEMOGRAPHICS ===\n")

# Helper: summarise a numeric variable
num_summary <- function(x, label) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    cat(sprintf("  %-35s : no data\n", label))
    return(invisible(NULL))
  }
  cat(sprintf(
    "  %-35s : n=%d  mean=%.1f  SD=%.1f  median=%.1f  IQR=%.1f-%.1f  range=%.1f-%.1f\n",
    label, length(x), mean(x), sd(x), median(x),
    quantile(x, 0.25), quantile(x, 0.75), min(x), max(x)
  ))
}

# Helper: frequency table for a categorical variable
cat_summary <- function(df, col, label, top_n = 10) {
  cat(sprintf("\n  %s (%s):\n", label, col))
  tbl <- df %>%
    mutate(val = .data[[col]]) %>%
    mutate(val = ifelse(is.na(val) | val == "" | val == "NA", "_Missing_", val)) %>%
    count(val, sort = TRUE) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    slice_head(n = top_n)
  print(tbl, n = Inf)
}


# --- 7a. Intervention period -------------------------------------------------
cat("\n--- Intervention period (Nov 2024 - Oct 2025) ---\n")

num_summary(int_pop$birthweight_num,  "Birthweight (g)")
num_summary(int_pop$gestation_num,    "Gestational age (wks)")
num_summary(int_pop$matageyrs,        "Maternal age (yrs)")
num_summary(int_pop$age_num,          "Age at admission (hrs)")

cat_summary(int_pop, "bwgroup",     "Birthweight group")
cat_summary(int_pop, "gestgroup",   "Gestational age group")
cat_summary(int_pop, "gender",      "Sex")
cat_summary(int_pop, "modedelivery","Mode of delivery")
cat_summary(int_pop, "typebirth",   "Type of birth")
cat_summary(int_pop, "agecategory", "Age category at admission")
cat_summary(int_pop, "admreason",   "Primary admission reason")
cat_summary(int_pop, "admittedfrom","Admitted from")
cat_summary(int_pop, "ansteroids",  "Antenatal steroids")


# --- 7b. Historical period --------------------------------------------------
cat("\n--- Historical period (Jan 2022 - Oct 2024) ---\n")

num_summary(hist_pop$birthweight_num, "Birthweight (g)")
num_summary(hist_pop$gestation_num,   "Gestational age (wks)")
num_summary(hist_pop$age_num,         "Age at admission (hrs)")

cat_summary(hist_pop, "bwgroup",    "Birthweight group")
cat_summary(hist_pop, "gestgroup",  "Gestational age group")
cat_summary(hist_pop, "gender",     "Sex")


# -----------------------------------------------------------------------------
# 8. DESCRIPTIVE STATISTICS -- GOLDEN HOUR CLINICAL VARIABLES
# -----------------------------------------------------------------------------

cat("\n=== DESCRIPTIVE STATISTICS: GOLDEN HOUR CLINICAL VARIABLES ===\n")

for (pop_label in c("Intervention", "Historical")) {
  pop <- if (pop_label == "Intervention") int_pop else hist_pop
  cat(sprintf("\n--- %s period ---\n", pop_label))

  # Cry at birth and triage
  cat_summary(pop, "crybirth",      "Cried at birth")
  cat_summary(pop, "babycrytriage", "Baby crying at triage (TRUE/FALSE)")

  # Resuscitation
  cat_summary(pop, "resus", "Resuscitation at birth (multi-select)", top_n = 15)
  cat(sprintf(
    "  Any resuscitation (non-None/Norm): %d / %d (%.1f%%)\n",
    sum(pop$resus_any == TRUE, na.rm = TRUE),
    sum(!is.na(pop$resus_any)),
    100 * mean(pop$resus_any == TRUE, na.rm = TRUE)
  ))

  # DCC eligibility proxy
  cat(sprintf(
    "  DCC-eligible proxy (apgar1>6 OR no resus): %d / %d (%.1f%%)\n",
    sum(pop$dcc_eligible_proxy == TRUE, na.rm = TRUE),
    sum(!is.na(pop$dcc_eligible_proxy)),
    100 * mean(pop$dcc_eligible_proxy == TRUE, na.rm = TRUE)
  ))

  # Temperature
  num_summary(pop$temperature_num, "Temperature at admission (C)")
  cat_summary(pop, "tempgroup",    "Temperature group (1-degree bands)")
  cat_summary(pop, "temp_cat",     "Temperature category (WHO)")
  cat(sprintf(
    "  Normal temperature (36.5-37.5 C): %d / %d (%.1f%%)\n",
    sum(pop$temp_normal == TRUE, na.rm = TRUE),
    sum(!is.na(pop$temp_normal)),
    100 * mean(pop$temp_normal == TRUE, na.rm = TRUE)
  ))

  # Apgar
  num_summary(pop$apgar1_num,  "Apgar 1 min")
  num_summary(pop$apgar5_num,  "Apgar 5 min")
  num_summary(pop$apgar10_num, "Apgar 10 min")
  cat(sprintf(
    "  Apgar 1 min <= 6 (low): %d / %d (%.1f%%)\n",
    sum(pop$apgar1_num <= 6, na.rm = TRUE),
    sum(!is.na(pop$apgar1_num)),
    100 * mean(pop$apgar1_num <= 6, na.rm = TRUE)
  ))

  # Vitals
  num_summary(pop$hr_num,      "Heart rate (bpm)")
  num_summary(pop$rr_num,      "Respiratory rate (bpm)")
  num_summary(pop$satsair_num, "SpO2 on air (%)")

  # Blood glucose (very sparse)
  cat(sprintf(
    "  Blood glucose (mmol/L) recorded: %d / %d (%.1f%%)  [SPARSE -- see note]\n",
    sum(!is.na(pop$bloodsugarmmol_num)),
    nrow(pop),
    100 * mean(!is.na(pop$bloodsugarmmol_num))
  ))
  num_summary(pop$bloodsugarmmol_num, "Blood glucose (mmol/L)")
  cat_summary(pop, "bsmonyn", "Blood sugar monitoring done")

  # Feeding at admission
  cat_summary(pop, "feedsadm", "Feeding type at admission")
}


# -----------------------------------------------------------------------------
# 9. DESCRIPTIVE STATISTICS -- OUTCOMES
# -----------------------------------------------------------------------------

cat("\n=== DESCRIPTIVE STATISTICS: OUTCOMES ===\n")

for (pop_label in c("Intervention", "Historical")) {
  pop      <- if (pop_label == "Intervention") int_pop else hist_pop
  matched  <- pop %>% filter(matched == TRUE)

  cat(sprintf("\n--- %s period ---\n", pop_label))
  cat(sprintf("  Total inborn             : %d\n", nrow(pop)))
  cat(sprintf("  Matched (discharge data) : %d (%.1f%%)\n",
              nrow(matched), 100 * nrow(matched) / nrow(pop)))

  cat_summary(matched, "neotreeoutcome", "Neotree outcome")
  cat(sprintf(
    "  Neonatal mortality (NND) : %d / %d matched (%.1f%%)\n",
    sum(matched$neotreeoutcome == "NND", na.rm = TRUE),
    nrow(matched),
    100 * mean(matched$neotreeoutcome == "NND", na.rm = TRUE)
  ))

  # Cause of death
  if (sum(matched$neotreeoutcome == "NND", na.rm = TRUE) > 0) {
    cat_summary(
      matched %>% filter(neotreeoutcome == "NND"),
      "causedeath", "Cause of death (NND only)", top_n = 12
    )
  }

  # Discharge temperature
  dischtemp_num <- suppressWarnings(as.numeric(matched$dischtemp))
  num_summary(dischtemp_num, "Temperature at discharge (C)")
}


# -----------------------------------------------------------------------------
# 10. MONTHLY TRENDS -- TEMPERATURE AND ADMISSIONS (Q10 SETUP)
# -----------------------------------------------------------------------------

cat("\n=== MONTHLY TRENDS (Q10 SETUP -- SMCH inborn, 2022-2025) ===\n")

monthly_q10 <- q10_pop %>%
  group_by(adm_month, period) %>%
  summarise(
    n_admissions  = n(),
    n_matched     = sum(matched, na.rm = TRUE),
    n_temp        = sum(!is.na(temperature_num)),
    pct_normal_temp = round(100 * mean(temp_normal == TRUE, na.rm = TRUE), 1),
    mean_temp     = round(mean(temperature_num, na.rm = TRUE), 2),
    pct_nnd       = round(
      100 * mean(neotreeoutcome == "NND", na.rm = TRUE), 1),
    n_dcc_eligible = sum(dcc_eligible_proxy == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(adm_month)

cat("\nMonthly summary table:\n")
print(monthly_q10, n = Inf)

write_csv(monthly_q10,
          file.path(OUTPUT_DIR, "01b_monthly_trends_q10.csv"))
cat("\nSaved: 01b_monthly_trends_q10.csv\n")


# -----------------------------------------------------------------------------
# 11. SUBGROUP: LBW BABIES (RELEVANT TO LBW SUBGROUP ANALYSES)
# -----------------------------------------------------------------------------

cat("\n=== LBW SUBGROUP (SMCH inborn, intervention period) ===\n")

lbw_int <- int_pop %>%
  filter(bwgroup %in% c("LBW", "VLBW", "ELBW"))

cat(sprintf("  LBW babies (any LBW category): %d / %d (%.1f%%)\n",
            nrow(lbw_int), nrow(int_pop),
            100 * nrow(lbw_int) / nrow(int_pop)))

cat_summary(lbw_int, "bwgroup", "LBW subgroup breakdown")

cat(sprintf("  DCC-eligible proxy (LBW): %d / %d (%.1f%%)\n",
            sum(lbw_int$dcc_eligible_proxy == TRUE, na.rm = TRUE),
            sum(!is.na(lbw_int$dcc_eligible_proxy)),
            100 * mean(lbw_int$dcc_eligible_proxy == TRUE, na.rm = TRUE)))

num_summary(lbw_int$temperature_num, "Temperature at admission (LBW, C)")
cat(sprintf("  Normal temperature (LBW): %d / %d (%.1f%%)\n",
            sum(lbw_int$temp_normal == TRUE, na.rm = TRUE),
            sum(!is.na(lbw_int$temp_normal)),
            100 * mean(lbw_int$temp_normal == TRUE, na.rm = TRUE)))


# -----------------------------------------------------------------------------
# 12. SUMMARY TABLE: KEY COUNTS FOR THE MANUSCRIPT
# -----------------------------------------------------------------------------

cat("\n=== KEY NUMBERS SUMMARY ===\n")

summary_tbl <- tibble(
  metric = c(
    "Full master dataset (all facilities)",
    "SMCH total records",
    "SMCH inborn non-readmissions (all time)",
    "SMCH inborn -- intervention period (Nov2024-Oct2025)",
    "  of which matched to discharge",
    "  of which with outcome = DC",
    "  of which with outcome = NND",
    "SMCH inborn -- historical period (Jan2022-Oct2024)",
    "  of which matched to discharge",
    "Int period: delivinter recorded (any value incl. U)",
    "Int period: delivinter known (excl. U and missing)",
    "Int period: DCC received (of known)",
    "Int period: S2S received (of known)",
    "Int period: DCC + S2S both (of known)",
    "Int period: no intervention (of known)"
  ),
  n = c(
    nrow(df),
    nrow(all_smch),
    nrow(smch_inborn),
    nrow(int_pop),
    sum(int_pop$matched, na.rm = TRUE),
    sum(int_pop$neotreeoutcome == "DC", na.rm = TRUE),
    sum(int_pop$neotreeoutcome == "NND", na.rm = TRUE),
    nrow(hist_pop),
    sum(hist_pop$matched, na.rm = TRUE),
    sum(!is.na(int_pop$delivinter) & trimws(int_pop$delivinter) != ""),
    sum(int_pop$delivinter_known, na.rm = TRUE),
    sum(int_pop$dcc_received, na.rm = TRUE),
    sum(int_pop$s2s_received, na.rm = TRUE),
    sum(int_pop$dcc_received & int_pop$s2s_received, na.rm = TRUE),
    sum(int_pop$no_intervention, na.rm = TRUE)
  )
)

print(summary_tbl, n = Inf)

write_csv(summary_tbl,
          file.path(OUTPUT_DIR, "01c_key_counts_summary.csv"))
cat("\nSaved: 01c_key_counts_summary.csv\n")

cat("\n=== Script 01 complete ===\n")

# -- Close log -----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
