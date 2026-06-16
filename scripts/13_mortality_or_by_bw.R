# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 13: Mortality ODDS RATIO by intervention receipt, stratified by BW
# =============================================================================
# Purpose (Action 3, 2026-06-09 meeting):
#   Per-baby comparison of neonatal mortality (NND) between babies who DID and
#   DID NOT receive DCC / S2S(ESSC) during the intervention period, adjusted
#   for birthweight, gestation, Apgar-1 and admission temperature, and
#   STRATIFIED BY BIRTH WEIGHT (All / LBW / NBW).
#
#   This revives and updates the analysis that was shown at the 19-May meeting
#   (originally in Script 02, since retired as the primary method). Rachel asked
#   to revisit it specifically for the LBW subgroup, where uptake was much lower
#   (~20% rising to ~60%) so a per-baby contrast is more informative than the
#   population-level ITS.
#
#   IMPORTANT INTERPRETIVE CAVEAT (carry into any reporting):
#   This is observational and subject to INDICATION BIAS -- healthier babies are
#   both more likely to receive DCC/S2S (not needing immediate resuscitation)
#   and more likely to survive. Adjustment reduces but does not remove this.
#   The ITS (Scripts 04b/12) remains the primary QI evaluation method. Present
#   these ORs as supportive/secondary, not causal.
#
# MODEL (Action from 2026-06-16 meeting -- report UNADJUSTED and ADJUSTED side by side):
#   Unadjusted : NND ~ exposure
#   Adjusted   : NND ~ exposure + birthweight + gestation + apgar1 + temperature
#   Both binomial logistic regressions are fitted on the SAME complete-case rows
#   (i.e. babies with all adjustment covariates present), so the crude and
#   adjusted ORs are directly comparable on an identical sample. OR =
#   exp(beta_exposure) with profile-likelihood 95% CI. Exposure is compared
#   against the "no intervention" group only (delivinter recorded, not DCC/S2S/BF).
#
#   Rationale (Rachel, 2026-06-16): showing both demonstrates that the protective
#   association persists after adjusting for the factors a midwife uses to judge
#   how well a baby is (birthweight, gestation, Apgar) -- the expected pattern is
#   a larger crude effect that attenuates but stays in the same direction.
#
# EXPOSURES:   DCC vs none, S2S vs none, Any-GH (DCC or S2S) vs none
# STRATA:      All inborn / LBW (<2500 g) / NBW (>=2500 g)
#
# POPULATION:
#   Source   : master NeoTree file (ZIM_db_master_joined_to_20260401.csv).
#              NEONATAL UNIT ADMISSIONS only (not the maternal script).
#   Filter   : facility == "SMCH", inborn == TRUE, non-readmission,
#              match_type == "direct_match" (NND needs discharge data),
#              delivinter recorded and not "U".
#   Window   : intervention period Nov 2024 - Mar 2026 (extended; the 19-May
#              version used Nov 2024 - Oct 2025 -- change INT_END to reproduce).
#
# NOTE ON RESULTS vs MAY:
#   Extending the window weakens the DCC signal (was OR ~0.65, p~0.02 at 12
#   months; ~0.83, p~0.14 at 17 months) while S2S remains protective,
#   especially in LBW. State the window explicitly when reporting.
#
# OUTPUTS (all to 03-OUTPUTS/):
#   13a_mortality_or_full_models.csv   -- every coefficient; 'model' column =
#                                         unadjusted | adjusted
#   13b_mortality_or_summary.csv       -- per stratum/exposure: crude_OR + adj_OR
#                                         side by side (headline table)
#   13c_mortality_or_counts.csv        -- group sizes / NND counts / crude rates
#   13_mortality_or_by_bw_log.txt
#
# DSH note: ASCII-only.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(broom)

args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag),
                                      mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

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

LOG_FILE <- file.path(OUTPUT_DIR, "13_mortality_or_by_bw_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2026-03-31")   # extended; set 2025-10-31 to match 19-May

cat("=============================================================\n")
cat("  Script 13: Mortality OR by intervention receipt, by BW stratum\n")
cat("=============================================================\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat(sprintf("  Data   : %s\n", DATA_PATH))
cat("=============================================================\n\n")


# -----------------------------------------------------------------------------
# 1. LOAD AND DERIVE
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)
cat(sprintf("  Rows loaded: %d\n", nrow(raw)))

df <- raw %>%
  mutate(
    adm_dt   = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date = as.Date(adm_dt),

    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                             ~ NA
    ),
    is_readmission = readmission == "Y",
    matched        = match_type == "direct_match",

    birthweight_num = suppressWarnings(as.numeric(birthweight)),
    gestation_num   = suppressWarnings(as.numeric(gestation)),
    apgar1_num      = suppressWarnings(as.numeric(apgar1)),
    temperature_num = suppressWarnings(as.numeric(temperature)),

    nnd = neotreeoutcome == "NND",

    dcc_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"           ~ NA,
      grepl("DCC", delivinter, fixed = TRUE) ~ TRUE,
      TRUE                                   ~ FALSE
    ),
    s2s_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"           ~ NA,
      grepl("S2S", delivinter, fixed = TRUE) ~ TRUE,
      TRUE                                   ~ FALSE
    ),
    any_gh = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"           ~ NA,
      grepl("DCC|S2S", delivinter)           ~ TRUE,
      TRUE                                   ~ FALSE
    ),
    no_intervention = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"           ~ NA,
      grepl("DCC|S2S|BF", delivinter)        ~ FALSE,
      TRUE                                   ~ TRUE
    ),
    lbw = !is.na(birthweight_num) & birthweight_num < 2500,
    nbw = !is.na(birthweight_num) & birthweight_num >= 2500
  )

pop <- df %>%
  filter(facility == "SMCH",
         inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission),
         matched == TRUE,
         adm_date >= INT_START, adm_date <= INT_END)

cat(sprintf("\nIntervention NICU admissions (SMCH inborn, matched): %d\n", nrow(pop)))
cat(sprintf("  delivinter known : %d\n",
            sum(!is.na(pop$no_intervention))))
cat(sprintf("  LBW / NBW        : %d / %d\n",
            sum(pop$lbw, na.rm = TRUE), sum(pop$nbw, na.rm = TRUE)))

COVARS <- c("birthweight_num", "gestation_num", "apgar1_num", "temperature_num")


# -----------------------------------------------------------------------------
# 2. OR HELPER
# -----------------------------------------------------------------------------

run_or <- function(data, exposure_col, stratum_label, exposure_label) {
  d <- data %>%
    filter(.data[[exposure_col]] == TRUE | no_intervention == TRUE) %>%
    mutate(exposure = factor(if_else(.data[[exposure_col]] == TRUE, "Yes", "No"),
                             levels = c("No", "Yes")))

  # Complete cases for the ADJUSTED model. The unadjusted model is fitted on the
  # SAME rows so crude and adjusted ORs are directly comparable (identical n).
  md <- d %>% select(nnd, exposure, all_of(COVARS)) %>% drop_na()

  n_exp  <- sum(md$exposure == "Yes")
  n_unexp <- sum(md$exposure == "No")
  n_nnd  <- sum(md$nnd, na.rm = TRUE)

  label <- sprintf("%s | %s", stratum_label, exposure_label)
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Complete cases: %d (exposed %d, unexposed %d), NND events: %d\n",
              nrow(md), n_exp, n_unexp, n_nnd))

  counts <- tibble(
    stratum = stratum_label, comparison = exposure_label,
    n_total = nrow(md), n_exposed = n_exp, n_unexposed = n_unexp,
    n_nnd = n_nnd,
    crude_nnd_exposed   = round(100 * sum(md$nnd[md$exposure == "Yes"]) / max(n_exp, 1), 2),
    crude_nnd_unexposed = round(100 * sum(md$nnd[md$exposure == "No"])  / max(n_unexp, 1), 2)
  )

  na_summ <- tibble(
    stratum = stratum_label, comparison = exposure_label,
    n = nrow(md), n_nnd = n_nnd,
    crude_OR = NA_real_, crude_CI_lower = NA_real_,
    crude_CI_upper = NA_real_, crude_p = NA_real_,
    adj_OR = NA_real_, adj_CI_lower = NA_real_,
    adj_CI_upper = NA_real_, adj_p = NA_real_,
    note = "insufficient data"
  )

  if (n_nnd < 5 || nrow(md) < 30 || n_exp < 10 || n_unexp < 10) {
    cat("  >> Too few events/observations for a stable model. Skipped.\n")
    return(list(summary = na_summ, full = tibble(), counts = counts))
  }

  # Helper: extract the full OR table from a fitted model, tagging the model.
  or_tab <- function(fit, model_label) {
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(stratum = stratum_label, comparison = exposure_label,
             model = model_label,
             OR = round(estimate, 3),
             CI_lower = round(conf.low, 3), CI_upper = round(conf.high, 3),
             p.value = round(p.value, 4)) %>%
      select(stratum, comparison, model, term, OR, CI_lower, CI_upper, p.value)
  }

  # Unadjusted (crude) and adjusted models on the identical complete-case sample
  fit_crude <- glm(nnd ~ exposure,
                   data = md, family = binomial(link = "logit"))
  fit_adj   <- glm(nnd ~ exposure + birthweight_num + gestation_num +
                     apgar1_num + temperature_num,
                   data = md, family = binomial(link = "logit"))

  full  <- bind_rows(or_tab(fit_crude, "unadjusted"),
                     or_tab(fit_adj,   "adjusted"))
  c_exp <- full %>% filter(model == "unadjusted", term == "exposureYes")
  a_exp <- full %>% filter(model == "adjusted",   term == "exposureYes")

  cat(sprintf("  Unadjusted OR: %.3f (95%% CI %.3f-%.3f)  p = %.4f\n",
              c_exp$OR, c_exp$CI_lower, c_exp$CI_upper, c_exp$p.value))
  cat(sprintf("  Adjusted   OR: %.3f (95%% CI %.3f-%.3f)  p = %.4f\n",
              a_exp$OR, a_exp$CI_lower, a_exp$CI_upper, a_exp$p.value))

  summ <- tibble(
    stratum = stratum_label, comparison = exposure_label,
    n = nrow(md), n_nnd = n_nnd,
    crude_OR = c_exp$OR, crude_CI_lower = c_exp$CI_lower,
    crude_CI_upper = c_exp$CI_upper, crude_p = c_exp$p.value,
    adj_OR = a_exp$OR, adj_CI_lower = a_exp$CI_lower,
    adj_CI_upper = a_exp$CI_upper, adj_p = a_exp$p.value,
    note = "adjusted for BW + gestation + apgar1 + adm temp"
  )

  list(summary = summ, full = full, counts = counts)
}


# -----------------------------------------------------------------------------
# 3. RUN ALL STRATA x EXPOSURES
# -----------------------------------------------------------------------------

strata <- list(
  "All inborn"     = pop,
  "LBW (<2500g)"   = pop %>% filter(lbw == TRUE),
  "NBW (>=2500g)"  = pop %>% filter(nbw == TRUE)
)
exposures <- list("DCC vs none" = "dcc_received",
                  "S2S vs none" = "s2s_received",
                  "Any GH vs none" = "any_gh")

res <- list()
for (sname in names(strata)) {
  for (ename in names(exposures)) {
    res[[paste(sname, ename)]] <-
      run_or(strata[[sname]], exposures[[ename]], sname, ename)
  }
}

summary_tbl <- bind_rows(lapply(res, `[[`, "summary"))
full_tbl    <- bind_rows(lapply(res, `[[`, "full"))
counts_tbl  <- bind_rows(lapply(res, `[[`, "counts"))

write_csv(full_tbl,    file.path(OUTPUT_DIR, "13a_mortality_or_full_models.csv"))
write_csv(summary_tbl, file.path(OUTPUT_DIR, "13b_mortality_or_summary.csv"))
write_csv(counts_tbl,  file.path(OUTPUT_DIR, "13c_mortality_or_counts.csv"))

cat("\n=============================================================\n")
cat("  HEADLINE OR SUMMARY -- UNADJUSTED vs ADJUSTED (received vs none)\n")
cat("=============================================================\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# Compact crude -> adjusted view for quick reading / manuscript text
cat("\nCrude -> adjusted (exposure OR, 95% CI, p):\n")
summary_tbl %>%
  filter(!is.na(adj_OR)) %>%
  mutate(line = sprintf(
    "  %-14s %-15s : crude %.2f (%.2f-%.2f) p=%.3f  ->  adj %.2f (%.2f-%.2f) p=%.3f",
    stratum, comparison,
    crude_OR, crude_CI_lower, crude_CI_upper, crude_p,
    adj_OR, adj_CI_lower, adj_CI_upper, adj_p)) %>%
  pull(line) %>% cat(sep = "\n")
cat("\n")

cat("\nSaved: 13a_mortality_or_full_models.csv  (unadjusted + adjusted, all coefficients)\n")
cat("       13b_mortality_or_summary.csv       (crude vs adjusted exposure OR, side by side)\n")
cat("       13c_mortality_or_counts.csv        (group sizes / NND counts / crude rates)\n")

cat("\n=== Script 13 complete ===\n")
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
