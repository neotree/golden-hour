# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 02: Outcome Analysis (Q8 and Q9)
# =============================================================================
# Purpose:
#   PRIMARY (Sections 5-6): Outcome analyses using actual DCC/ESSC exposure
#     from delivinter (introduced Oct 2024 at SMCH). Compares DCC-received vs
#     not-DCC, S2S-received vs not-S2S, and any-GH-intervention vs none, within
#     the intervention period. Includes LBW subgroup analyses.
#
#   SECONDARY (Sections 7-8): Period comparisons (intervention vs historical)
#     for NND (Q8) and normal admission temperature (Q9). These are secondary
#     because delivinter has no historical-period data, so the period comparison
#     cannot be exposure-based.
#
# Population notes:
#   NND is a discharge outcome -> analyses restricted to matched records
#   (match_type == "direct_match").
#   Admission temperature is recorded at admission -> temp analyses use all
#   inborn records (no discharge match required).
#   delivinter denominator: babies with delivinter_known == TRUE (i.e. delivinter
#   recorded and not "U"/Unknown). Historical period = 0 known records.
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(broom)     # for tidy() -- install.packages("broom") if needed

# -- Locate this script's directory -------------------------------------------
# Works when run via: Rscript 02_outcome_analysis.R
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
            "ZIM_db_master_joined_to_20260401.csv"),
  mustWork = FALSE
)

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -- Redirect all output to a log file ----------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "02_outcome_analysis_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2025-10-31")
HIST_START <- as.Date("2022-01-01")
HIST_END   <- as.Date("2024-10-31")


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
cat(sprintf("  Rows loaded: %d\n", nrow(raw)))


# -----------------------------------------------------------------------------
# 2. CLEANING AND DERIVED VARIABLES
# -----------------------------------------------------------------------------
# Reproduces the derivations from 01_descriptive_overview.R so this script
# is fully self-contained.

df <- raw %>%
  mutate(
    # -- Admission date/time --------------------------------------------------
    adm_dt    = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date  = as.Date(adm_dt),
    adm_month = floor_date(adm_dt, "month"),

    # -- Inborn flag (harmonised across script versions) ----------------------
    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                             ~ NA
    ),

    # -- Key clinical variables (numeric) -------------------------------------
    temperature_num = suppressWarnings(as.numeric(temperature)),
    apgar1_num      = suppressWarnings(as.numeric(apgar1)),
    birthweight_num = suppressWarnings(as.numeric(birthweight)),
    gestation_num   = suppressWarnings(as.numeric(gestation)),

    # -- Resuscitation: any non-None/Norm intervention ------------------------
    resus_any = case_when(
      resus %in% c("NONE", "Norm") ~ FALSE,
      is.na(resus) | resus == ""   ~ NA,
      TRUE                         ~ TRUE
    ),

    # -- DCC eligibility proxy ------------------------------------------------
    # Used for secondary period-comparison subgroup analyses.
    # Primary analyses use actual delivinter exposure (see below).
    dcc_eligible_proxy = case_when(
      is.na(apgar1_num) & is.na(resus_any) ~ NA,
      !is.na(apgar1_num) & apgar1_num > 6  ~ TRUE,
      resus_any == FALSE                    ~ TRUE,
      TRUE                                  ~ FALSE
    ),

    # -- Delivery interventions (delivinter): DCC, S2S, BF -------------------
    # Introduced Oct 2024; zero records in the historical period.
    # Multi-select in two encoding formats: {DCC,S2S} or ["DCC","S2S"].
    # grepl() handles both. "U" = Unknown (NA in all derived vars below).
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
    delivinter_known = (
      !is.na(delivinter) &
      trimws(delivinter) != "" &
      trimws(delivinter) != "U"
    ),

    # -- Normal admission temperature (WHO: 36.5-37.5 C) ---------------------
    temp_normal = case_when(
      is.na(temperature_num)                             ~ NA,
      temperature_num >= 36.5 & temperature_num <= 37.5 ~ TRUE,
      TRUE                                               ~ FALSE
    ),

    # -- Study period ---------------------------------------------------------
    period = case_when(
      adm_date >= INT_START  & adm_date <= INT_END   ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END  ~ "Historical",
      TRUE                                            ~ "Other"
    ),

    # -- Match and readmission flags ------------------------------------------
    matched        = match_type == "direct_match",
    is_readmission = readmission == "Y"
  )

cat("Derived variables created.\n")


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
cat(sprintf("  Intervention period (inborn)  : %d\n", nrow(int_pop)))
cat(sprintf("    of which matched            : %d\n",
            sum(int_pop$matched, na.rm = TRUE)))
cat(sprintf("  Historical period (inborn)    : %d\n", nrow(hist_pop)))
cat(sprintf("    of which matched            : %d\n",
            sum(hist_pop$matched, na.rm = TRUE)))
cat(sprintf("  Q10 combined (inborn)         : %d\n", nrow(q10_pop)))


# -----------------------------------------------------------------------------
# 4. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# -- chi_compare --------------------------------------------------------------
# Unadjusted 2x2 chi-square comparing a binary outcome (TRUE/FALSE column)
# by period (Historical vs Intervention).
# Returns a tidy data frame with rates, test statistic, and p-value.
chi_compare <- function(data, outcome_col, subgroup_label) {
  cat(sprintf("\n--- Chi-square: %s ---\n", subgroup_label))

  # Exclude records with NA outcome
  d <- data %>% filter(!is.na(.data[[outcome_col]]))

  # Contingency table: period (rows) x outcome (cols)
  tbl <- table(
    period  = factor(d$period, levels = c("Historical", "Intervention")),
    outcome = factor(d[[outcome_col]], levels = c(FALSE, TRUE),
                     labels = c("No", "Yes"))
  )
  cat("Contingency table (rows = period, cols = outcome):\n")
  print(tbl)
  cat(sprintf("  Row totals: Historical = %d, Intervention = %d\n",
              sum(tbl["Historical", ]), sum(tbl["Intervention", ])))

  test <- chisq.test(tbl, correct = FALSE)
  cat(sprintf("  Chi-square = %.4f  df = %d  p-value = %.4f\n",
              test$statistic, test$parameter, test$p.value))

  # Per-period rates
  rates <- d %>%
    group_by(period) %>%
    summarise(
      n_total   = n(),
      n_event   = sum(.data[[outcome_col]] == TRUE, na.rm = TRUE),
      rate_pct  = round(100 * n_event / n_total, 2),
      .groups   = "drop"
    ) %>%
    mutate(
      subgroup = subgroup_label,
      chi_sq   = round(as.numeric(test$statistic), 4),
      df       = as.integer(test$parameter),
      p_value  = round(test$p.value, 4)
    ) %>%
    arrange(factor(period, levels = c("Historical", "Intervention")))

  print(as.data.frame(rates))
  rates
}


# -- logistic_model -----------------------------------------------------------
# Adjusted logistic regression.
# Returns a list: $model (glm object) and $results (tidy OR table).
# period is always set to factor with Historical as reference.
logistic_model <- function(data, outcome_col, covariate_cols, analysis_label) {
  cat(sprintf("\n--- Logistic regression: %s ---\n", analysis_label))

  # Complete cases for all model variables
  model_vars <- c(outcome_col, covariate_cols)
  md <- data %>%
    select(all_of(model_vars)) %>%
    drop_na()

  cat(sprintf("  Observations with complete model data: %d / %d\n",
              nrow(md), nrow(data)))

  # Set period factor with Historical as reference
  if ("period" %in% covariate_cols) {
    md <- md %>%
      mutate(period = factor(period, levels = c("Historical", "Intervention")))
  }

  fmla <- as.formula(
    paste(outcome_col, "~", paste(covariate_cols, collapse = " + "))
  )
  cat(sprintf("  Formula : %s\n", deparse(fmla)))

  fit <- glm(fmla, data = md, family = binomial(link = "logit"))

  # Tidy table with ORs (exponentiated) and 95% profile-likelihood CIs
  results <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    rename(OR = estimate, CI_lower = conf.low, CI_upper = conf.high) %>%
    mutate(
      OR       = round(OR, 3),
      CI_lower = round(CI_lower, 3),
      CI_upper = round(CI_upper, 3),
      p.value  = round(p.value, 4),
      analysis = analysis_label
    ) %>%
    select(analysis, term, OR, CI_lower, CI_upper, p.value)

  cat("\n  Full model OR table:\n")
  print(as.data.frame(results), digits = 3)

  # Highlight the period coefficient specifically
  period_row <- results %>% filter(grepl("period", term))
  if (nrow(period_row) > 0) {
    cat(sprintf(
      "\n  *** Period effect (Intervention vs Historical):\n"
    ))
    cat(sprintf(
      "      OR = %.3f  (95%% CI: %.3f - %.3f)  p = %.4f\n",
      period_row$OR, period_row$CI_lower, period_row$CI_upper,
      period_row$p.value
    ))
  }

  list(model = fit, results = results)
}


# =============================================================================
# 5. PRIMARY OUTCOME ANALYSIS: DCC / S2S EXPOSURE vs NND
#    (intervention period only; denominator = delivinter_known)
# =============================================================================
# Compares outcomes between babies who received each Golden Hour intervention
# and those who received none, within the intervention period. This is the
# primary analysis for the DCC- and S2S-specific outcome questions.
#
# Population: int_pop, matched (NND requires discharge data), delivinter_known.
# Exclusion: delivinter == "U" or missing (exposure status unknown).
# Control group: no_intervention == TRUE (delivinter = NONE or Norm).

cat("\n")
cat("================================================================\n")
cat("  PRIMARY OUTCOME: DCC/S2S EXPOSURE vs NND (int period only)\n")
cat("================================================================\n")

# -- Subset: intervention period, matched, delivinter known -------------------
int_matched_known <- int_pop %>%
  filter(matched == TRUE, delivinter_known == TRUE) %>%
  mutate(nnd = neotreeoutcome == "NND")

cat(sprintf(
  "\nIntervention period, matched, delivinter known: %d records\n",
  nrow(int_matched_known)
))
cat(sprintf("  DCC received          : %d\n",
            sum(int_matched_known$dcc_received, na.rm = TRUE)))
cat(sprintf("  S2S received          : %d\n",
            sum(int_matched_known$s2s_received, na.rm = TRUE)))
cat(sprintf("  DCC + S2S             : %d\n",
            sum(int_matched_known$dcc_received &
                int_matched_known$s2s_received, na.rm = TRUE)))
cat(sprintf("  Any GH intervention   : %d\n",
            sum(int_matched_known$any_gh_intervention, na.rm = TRUE)))
cat(sprintf("  No intervention       : %d\n",
            sum(int_matched_known$no_intervention, na.rm = TRUE)))

# Helper: chi-square comparing a binary exposure (TRUE/FALSE) vs a binary outcome.
# outcome_col: name of the logical outcome column in data.
# outcome_labels: c("event_absent_label", "event_present_label")
exposure_chi <- function(data, exposure_col, outcome_col, label,
                         outcome_labels = c("No event", "Event")) {
  d <- data %>%
    filter(!is.na(.data[[exposure_col]]), !is.na(.data[[outcome_col]]))
  tbl <- table(
    exposure = factor(d[[exposure_col]], levels = c(FALSE, TRUE),
                      labels = c("No", "Yes")),
    outcome  = factor(d[[outcome_col]], levels = c(FALSE, TRUE),
                      labels = outcome_labels)
  )
  cat(sprintf("\n--- %s ---\n", label))
  print(tbl)
  rates <- d %>%
    group_by(exposed = .data[[exposure_col]]) %>%
    summarise(n = n(),
              n_event = sum(.data[[outcome_col]], na.rm = TRUE),
              pct_event = round(100 * n_event / n, 2),
              .groups = "drop")
  print(rates)
  test <- chisq.test(tbl, correct = FALSE)
  cat(sprintf("  Chi-sq = %.3f  p = %.4f\n", test$statistic, test$p.value))
  ev_col <- outcome_labels[2]
  no_col <- outcome_labels[1]
  tibble(exposure = label,
         outcome  = outcome_col,
         yes_n    = sum(tbl["Yes", ]),
         yes_nev  = tbl["Yes", ev_col],
         yes_pct  = round(100 * tbl["Yes", ev_col] / sum(tbl["Yes",]), 2),
         no_n     = sum(tbl["No", ]),
         no_nev   = tbl["No",  ev_col],
         no_pct   = round(100 * tbl["No",  ev_col] / sum(tbl["No",]),  2),
         chi_sq   = round(as.numeric(test$statistic), 3),
         p_value  = round(test$p.value, 4))
}

# Helper: logistic model for exposure vs a binary outcome.
# outcome_col must be a logical (TRUE/FALSE) column.
exposure_logistic <- function(data, exposure_col, outcome_col,
                              covariate_cols, label) {
  d <- data %>%
    select(all_of(c(outcome_col, exposure_col, covariate_cols))) %>%
    drop_na() %>%
    mutate(exposure = factor(.data[[exposure_col]],
                             levels = c(FALSE, TRUE),
                             labels = c("No", "Yes")))
  cat(sprintf("\n  Logistic regression: %s\n", label))
  cat(sprintf("  Complete cases: %d / %d\n", nrow(d), nrow(data)))
  covars_str <- paste(c("exposure", covariate_cols), collapse = " + ")
  fmla <- as.formula(paste(outcome_col, "~", covars_str))
  fit  <- glm(fmla, data = d, family = binomial(link = "logit"))
  res  <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    rename(OR = estimate, CI_lower = conf.low, CI_upper = conf.high) %>%
    mutate(across(c(OR, CI_lower, CI_upper), ~round(., 3)),
           p.value  = round(p.value, 4),
           analysis = label) %>%
    select(analysis, term, OR, CI_lower, CI_upper, p.value)
  exp_row <- res %>% filter(grepl("exposureYes", term))
  if (nrow(exp_row) > 0) {
    cat(sprintf("  OR (exposed vs unexposed): %.3f (95%% CI %.3f-%.3f)  p=%.4f\n",
                exp_row$OR, exp_row$CI_lower, exp_row$CI_upper, exp_row$p.value))
  }
  res
}

gh_covars <- c("birthweight_num", "gestation_num", "apgar1_num", "temperature_num")

# -- 5a. DCC vs no intervention -----------------------------------------------
cat("\n=== 5a. DCC received vs no intervention: NND ===\n")
dcc_vs_none <- int_matched_known %>%
  filter(dcc_received == TRUE | no_intervention == TRUE)
cat(sprintf("DCC vs None pool: %d records (DCC=%d, None=%d)\n",
            nrow(dcc_vs_none),
            sum(dcc_vs_none$dcc_received), sum(dcc_vs_none$no_intervention)))

q4_chi_dcc <- exposure_chi(dcc_vs_none, "dcc_received", "nnd", "DCC vs None: NND",
                           outcome_labels = c("Survived", "NND"))
q4_lr_dcc  <- exposure_logistic(dcc_vs_none, "dcc_received", "nnd", gh_covars,
                                "DCC vs None: NND adjusted")

# -- 5b. S2S vs no intervention -----------------------------------------------
cat("\n=== 5b. S2S received vs no intervention: NND ===\n")
s2s_vs_none <- int_matched_known %>%
  filter(s2s_received == TRUE | no_intervention == TRUE)
cat(sprintf("S2S vs None pool: %d records (S2S=%d, None=%d)\n",
            nrow(s2s_vs_none),
            sum(s2s_vs_none$s2s_received), sum(s2s_vs_none$no_intervention)))

q5_chi_s2s <- exposure_chi(s2s_vs_none, "s2s_received", "nnd", "S2S vs None: NND",
                           outcome_labels = c("Survived", "NND"))
q5_lr_s2s  <- exposure_logistic(s2s_vs_none, "s2s_received", "nnd", gh_covars,
                                "S2S vs None: NND adjusted")

# -- 5c. Any GH intervention vs none ------------------------------------------
cat("\n=== 5c. Any GH intervention (DCC or S2S) vs none: NND ===\n")
any_vs_none <- int_matched_known %>%
  filter(any_gh_intervention == TRUE | no_intervention == TRUE)
cat(sprintf("Any-GH vs None pool: %d records (any=%d, none=%d)\n",
            nrow(any_vs_none),
            sum(any_vs_none$any_gh_intervention), sum(any_vs_none$no_intervention)))

q6_chi_any <- exposure_chi(any_vs_none, "any_gh_intervention", "nnd",
                           "Any GH vs None: NND",
                           outcome_labels = c("Survived", "NND"))
q6_lr_any  <- exposure_logistic(any_vs_none, "any_gh_intervention", "nnd",
                                gh_covars, "Any GH vs None: NND adjusted")

# -- 5d. LBW subgroup: DCC vs no intervention ---------------------------------
cat("\n=== 5d. LBW subgroup: DCC vs no intervention: NND ===\n")
lbw_dcc_vs_none <- dcc_vs_none %>%
  filter(bwgroup %in% c("LBW", "VLBW", "ELBW"))
cat(sprintf("LBW DCC vs None: %d records\n", nrow(lbw_dcc_vs_none)))

if (nrow(lbw_dcc_vs_none) >= 10) {
  q7_chi_lbw <- exposure_chi(lbw_dcc_vs_none, "dcc_received", "nnd",
                             "LBW: DCC vs None: NND",
                             outcome_labels = c("Survived", "NND"))
  q7_lr_lbw  <- exposure_logistic(lbw_dcc_vs_none, "dcc_received", "nnd",
                                  gh_covars, "LBW: DCC vs None: NND adjusted")
} else {
  cat("  Insufficient records for LBW subgroup analysis.\n")
  q7_chi_lbw <- tibble()
  q7_lr_lbw  <- tibble()
}

# Save exposure chi-square results
q4567_chi <- bind_rows(q4_chi_dcc, q5_chi_s2s, q6_chi_any, q7_chi_lbw)
write_csv(q4567_chi,
          file.path(OUTPUT_DIR, "02a_exposure_nnd_unadjusted.csv"))
cat("\nSaved: 02a_exposure_nnd_unadjusted.csv\n")

# Save exposure logistic results
q4567_lr <- bind_rows(q4_lr_dcc, q5_lr_s2s, q6_lr_any, q7_lr_lbw)
write_csv(q4567_lr,
          file.path(OUTPUT_DIR, "02b_exposure_nnd_adjusted.csv"))
cat("\nSaved: 02b_exposure_nnd_adjusted.csv\n")


# =============================================================================
# 6. PRIMARY OUTCOME: DCC / S2S EXPOSURE vs NORMAL ADMISSION TEMPERATURE
# =============================================================================
# Temperature is an admission variable -- discharge match not required.
# Population: int_pop, delivinter_known.

cat("\n")
cat("================================================================\n")
cat("  PRIMARY OUTCOME: DCC/S2S EXPOSURE vs NORMAL ADM. TEMPERATURE\n")
cat("================================================================\n")

int_known_all <- int_pop %>%
  filter(delivinter_known == TRUE)

cat(sprintf("\nIntervention period, delivinter known (all inborn): %d\n",
            nrow(int_known_all)))

temp_covars <- c("birthweight_num", "gestation_num", "apgar1_num")

# -- 6a. DCC vs no intervention: normal temperature ---------------------------
cat("\n=== 6a. DCC vs no intervention: normal temperature ===\n")
dcc_vs_none_all <- int_known_all %>%
  filter(dcc_received == TRUE | no_intervention == TRUE)

q4t_chi <- exposure_chi(dcc_vs_none_all, "dcc_received", "temp_normal",
                        "DCC vs None: normal temp",
                        outcome_labels = c("Abnormal", "Normal"))
q4t_lr  <- exposure_logistic(dcc_vs_none_all, "dcc_received", "temp_normal",
                             temp_covars, "DCC vs None: normal temp adjusted")

# -- 6b. S2S vs no intervention: normal temperature ---------------------------
cat("\n=== 6b. S2S vs no intervention: normal temperature ===\n")
s2s_vs_none_all <- int_known_all %>%
  filter(s2s_received == TRUE | no_intervention == TRUE)

q5t_chi <- exposure_chi(s2s_vs_none_all, "s2s_received", "temp_normal",
                        "S2S vs None: normal temp",
                        outcome_labels = c("Abnormal", "Normal"))
q5t_lr  <- exposure_logistic(s2s_vs_none_all, "s2s_received", "temp_normal",
                             temp_covars, "S2S vs None: normal temp adjusted")

# -- 6c. Any GH vs none: normal temperature -----------------------------------
cat("\n=== 6c. Any GH vs none: normal temperature ===\n")
any_vs_none_all <- int_known_all %>%
  filter(any_gh_intervention == TRUE | no_intervention == TRUE)

q6t_chi <- exposure_chi(any_vs_none_all, "any_gh_intervention", "temp_normal",
                        "Any GH vs None: normal temp",
                        outcome_labels = c("Abnormal", "Normal"))
q6t_lr  <- exposure_logistic(any_vs_none_all, "any_gh_intervention", "temp_normal",
                             temp_covars, "Any GH vs None: normal temp adjusted")

# Save
bind_rows(q4t_chi, q5t_chi, q6t_chi) %>%
  write_csv(file.path(OUTPUT_DIR, "02c_exposure_temp_unadjusted.csv"))
cat("Saved: 02c_exposure_temp_unadjusted.csv\n")

bind_rows(q4t_lr, q5t_lr, q6t_lr) %>%
  write_csv(file.path(OUTPUT_DIR, "02d_exposure_temp_adjusted.csv"))
cat("Saved: 02d_exposure_temp_adjusted.csv\n")


# =============================================================================
# 7. SECONDARY: NND -- PERIOD COMPARISON (intervention vs historical)
# =============================================================================
# Retained as a secondary/sensitivity analysis. The period comparison captures
# all babies regardless of delivinter status and provides the historical context.

cat("\n")
cat("================================================================\n")
cat("  SECONDARY: NND -- PERIOD COMPARISON\n")
cat("================================================================\n")

# -- Matched combined dataset for period-comparison NND analysis --------------
# NND is a discharge-side outcome; restrict to direct_match records.
matched_combined <- q10_pop %>%
  filter(matched == TRUE) %>%
  mutate(nnd = neotreeoutcome == "NND")

cat(sprintf("\nMatched records for NND analysis: %d\n", nrow(matched_combined)))
cat(sprintf("  Historical  : %d matched\n",
            sum(matched_combined$period == "Historical")))
cat(sprintf("  Intervention: %d matched\n",
            sum(matched_combined$period == "Intervention")))

# Descriptive NND counts
cat("\nNND counts by period (matched records):\n")
matched_combined %>%
  group_by(period) %>%
  summarise(
    n_total = n(),
    n_nnd   = sum(nnd, na.rm = TRUE),
    nnd_pct = round(100 * n_nnd / n_total, 2),
    .groups = "drop"
  ) %>%
  arrange(factor(period, levels = c("Historical", "Intervention"))) %>%
  print()


# -- 5a. Unadjusted chi-square ------------------------------------------------

cat("\n--- Q8: Unadjusted analyses ---\n")

# Overall
q8_chi_overall <- chi_compare(
  matched_combined, "nnd",
  "Q8 NND - Overall (all matched, inborn)"
)

# DCC/ESSC-eligible subgroup
eligible_matched <- matched_combined %>% filter(dcc_eligible_proxy == TRUE)
cat(sprintf("\n  DCC-eligible subgroup (matched): %d records\n",
            nrow(eligible_matched)))

q8_chi_eligible <- chi_compare(
  eligible_matched, "nnd",
  "Q8 NND - DCC/ESSC-eligible subgroup (matched)"
)

# LBW subgroup (LBW + VLBW + ELBW)
lbw_matched <- matched_combined %>%
  filter(bwgroup %in% c("LBW", "VLBW", "ELBW"))
cat(sprintf("\n  LBW subgroup (matched): %d records\n", nrow(lbw_matched)))

q8_chi_lbw <- chi_compare(
  lbw_matched, "nnd",
  "Q8 NND - LBW subgroup (matched)"
)

# Combine chi-square results and save
q8_chi_all <- bind_rows(q8_chi_overall, q8_chi_eligible, q8_chi_lbw)
write_csv(q8_chi_all,
          file.path(OUTPUT_DIR, "02e_q8_nnd_period_unadjusted.csv"))
cat("\nSaved: 02e_q8_nnd_period_unadjusted.csv\n")


# -- 5b. Adjusted logistic regression -----------------------------------------

cat("\n--- Q8: Adjusted logistic regression ---\n")
cat("    Model: NND ~ period + birthweight_num + gestation_num +")
cat(" apgar1_num + temperature_num\n")

# Covariates for Q8 adjustment
q8_covars <- c("period", "birthweight_num", "gestation_num",
               "apgar1_num", "temperature_num")

q8_logistic <- logistic_model(
  matched_combined, "nnd", q8_covars,
  "Q8 adjusted: NND ~ period + BW + gestation + apgar1 + adm_temp"
)

write_csv(q8_logistic$results,
          file.path(OUTPUT_DIR, "02f_q8_nnd_period_adjusted.csv"))
cat("Saved: 02f_q8_nnd_period_adjusted.csv\n")


# =============================================================================
# 8. SECONDARY: NORMAL ADMISSION TEMPERATURE -- PERIOD COMPARISON
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  SECONDARY: NORMAL ADM. TEMPERATURE -- PERIOD COMPARISON\n")
cat("================================================================\n")

# Temperature is recorded at admission; discharge match is NOT required.
# Use all inborn records (q10_pop) for Q9.
# Note: birthweight group filter distinguishes LBW for the subgroup analyses.

cat(sprintf("\nAll inborn records for Q9 temperature analysis: %d\n",
            nrow(q10_pop)))
cat(sprintf("  Temperature recorded: %d / %d (%.1f%%)\n",
            sum(!is.na(q10_pop$temperature_num)),
            nrow(q10_pop),
            100 * mean(!is.na(q10_pop$temperature_num))))

# Descriptive counts
cat("\nNormal temperature counts by period (all inborn):\n")
q10_pop %>%
  filter(!is.na(temp_normal)) %>%
  group_by(period) %>%
  summarise(
    n_total     = n(),
    n_normal    = sum(temp_normal, na.rm = TRUE),
    normal_pct  = round(100 * n_normal / n_total, 2),
    .groups     = "drop"
  ) %>%
  arrange(factor(period, levels = c("Historical", "Intervention"))) %>%
  print()


# -- 6a. Unadjusted chi-square ------------------------------------------------

cat("\n--- Q9: Unadjusted analyses ---\n")

# Overall
q9_chi_overall <- chi_compare(
  q10_pop, "temp_normal",
  "Q9 Normal temp - Overall (all inborn)"
)

# DCC/ESSC-eligible subgroup
eligible_all <- q10_pop %>% filter(dcc_eligible_proxy == TRUE)
cat(sprintf("\n  DCC-eligible subgroup (all inborn): %d records\n",
            nrow(eligible_all)))

q9_chi_eligible <- chi_compare(
  eligible_all, "temp_normal",
  "Q9 Normal temp - DCC/ESSC-eligible subgroup"
)

# LBW subgroup
lbw_all <- q10_pop %>% filter(bwgroup %in% c("LBW", "VLBW", "ELBW"))
cat(sprintf("\n  LBW subgroup (all inborn): %d records\n", nrow(lbw_all)))

q9_chi_lbw <- chi_compare(
  lbw_all, "temp_normal",
  "Q9 Normal temp - LBW subgroup"
)

# Combine and save
q9_chi_all <- bind_rows(q9_chi_overall, q9_chi_eligible, q9_chi_lbw)
write_csv(q9_chi_all,
          file.path(OUTPUT_DIR, "02g_q9_temp_period_unadjusted.csv"))
cat("\nSaved: 02g_q9_temp_period_unadjusted.csv\n")


# -- 6b. Adjusted logistic regression -----------------------------------------

cat("\n--- Q9: Adjusted logistic regression ---\n")
cat("    Model: temp_normal ~ period + birthweight_num + gestation_num +")
cat(" apgar1_num\n")
cat("    (admission temperature is the outcome; excluded from predictors)\n")

# Covariates for Q9: note temperature_num is EXCLUDED (it is the outcome)
q9_covars <- c("period", "birthweight_num", "gestation_num", "apgar1_num")

q9_logistic <- logistic_model(
  q10_pop, "temp_normal", q9_covars,
  "Q9 adjusted: normal_temp ~ period + BW + gestation + apgar1"
)

write_csv(q9_logistic$results,
          file.path(OUTPUT_DIR, "02h_q9_temp_period_adjusted.csv"))
cat("Saved: 02h_q9_temp_period_adjusted.csv\n")


# =============================================================================
# 9. COMBINED RESULTS SUMMARY
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  COMBINED OR SUMMARY (for manuscript)\n")
cat("================================================================\n")

# Primary: exposure ORs
exposure_summary <- bind_rows(q4567_lr, q4t_lr, q5t_lr, q6t_lr) %>%
  filter(grepl("exposureYes", term)) %>%
  select(analysis, OR, CI_lower, CI_upper, p.value)

cat("\nPrimary exposure ORs:\n")
print(as.data.frame(exposure_summary))

write_csv(exposure_summary,
          file.path(OUTPUT_DIR, "02i_exposure_or_summary.csv"))
cat("Saved: 02i_exposure_or_summary.csv\n")

# Secondary: period ORs
period_summary <- bind_rows(
  q8_logistic$results %>% filter(grepl("Intervention", term)),
  q9_logistic$results %>% filter(grepl("Intervention", term))
) %>%
  select(analysis, term, OR, CI_lower, CI_upper, p.value)

cat("\nSecondary period ORs:\n")
print(as.data.frame(period_summary))

write_csv(period_summary,
          file.path(OUTPUT_DIR, "02j_period_or_summary.csv"))
cat("Saved: 02j_period_or_summary.csv\n")

cat("\n=== Script 02 complete ===\n")

# -- Close log -----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
