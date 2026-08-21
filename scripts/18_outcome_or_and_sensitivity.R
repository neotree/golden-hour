# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# JULY 2026 REQUESTS -- Script 18: OUTCOME ODDS RATIOS (Table 4) + SENSITIVITY
# =============================================================================
# The single inferential analysis agreed for the outcomes paper (Rachel 14 July
# 15:20, incorporating Michelle's request to add sex and Hannah's over-adjustment
# concern). Per-baby logistic regression of death before discharge (NND) by the
# four-group exposure, Neither as the reference category.
#
# MODELS (fitted on the SAME complete-case sample so all columns are comparable):
#   Unadjusted        : nnd ~ exposure
#   Adjusted (BW)     : nnd ~ exposure + birthweight
#   Adjusted (BW+sex) : nnd ~ exposure + birthweight + sex     <- primary adjusted
#   (gestation, admission temperature and Apgar are DELIBERATELY excluded --
#    Apgar/resus are handled by the eligibility exclusion, resolving Hannah's
#    collinearity concern. Both adjusted columns are reported per Rachel's email,
#    which asked for BW-only and then BW+sex.)
#
# EXPOSURE (delivinter, four mutually exclusive groups):
#   Neither (reference) / DCC only / ESSC(S2S) only / Both DCC and ESSC
#   ESSC-only is reported as COUNTS ONLY (Tables 2/3) and OMITTED from the OR
#   models -- team decision 21 Jul (group too small; single-figure deaths).
#
# STRATA / POPULATIONS (Table 4 + sensitivity):
#   1. All eligible infants
#   2. LBW < 2500 g
#   3. SENSITIVITY: 1500-2499 g (babies most likely to benefit; Rachel's request)
#
# POPULATION: SMCH inborn, non-readmission, Nov 2024-Oct 2025, ELIGIBLE, with
#   delivinter known and direct_match (a discharge outcome is required for NND).
# ELIGIBILITY (team decision 21 Jul -- Marcia + Rachel): exclude 5-minute Apgar < 7
#   (apgar5). The resus/resustype limb is DROPPED (resus incompletely filled).
#   This supersedes the earlier "apgar1 < 6 OR resustype BVM/CPR" draft.
#
# INTERPRETIVE CAVEAT (carry into the manuscript): observational; residual
#   confounding by indication persists even after the eligibility exclusion and
#   adjustment. Report as supportive/secondary observational evidence, not causal.
#
# DATA-QUALITY HANDLING: see header of Script 17 / the variable-issues catalogue
#   (birthweight_g preference + range filter; sex normalisation; robust
#   delivinter/resustype matching).
#
# OUTPUTS (to ./outputs):
#   18a_or_full_models.csv   -- every coefficient, tagged model + stratum
#   18b_or_table4.csv        -- Table 4 layout: OR (95% CI) per group x model x stratum
#   18c_or_counts.csv        -- group sizes / death counts / crude rates per stratum
#   18_outcome_or_log.txt
#
# DSH note: ASCII-only.
# =============================================================================

library(tidyverse)
library(lubridate)
library(broom)

# -----------------------------------------------------------------------------
# 0. SETUP + TOGGLES
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
# Outputs are tagged by the data-file date so a cross-check run against the
# older master lands in its own subfolder (e.g. outputs/20260525, outputs/20260401).
OUTPUT_TAG <- sub("^.*_(\\d{8})\\.csv$", "\\1", DATA_FILE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", paste0(OUTPUT_TAG, "_apgar1"))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "18_outcome_or_log.txt")
log_con  <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")

# ELIGIBILITY (team decision, 21 Jul 2026 -- Marcia + Rachel):
#   Exclude babies with a LOW APGAR; DROP the resus/resustype limb, because the
#   resus fields were incompletely filled (assume undocumented = not needed).
#   Marcia specified the 5-MINUTE Apgar with a < 7 cut-off (this REPLACES the
#   1-minute Apgar < 6 in the first draft).
APGAR_FIELD    <- "apgar1"   # AUG 2026: 1-minute Apgar (Rachel 28 Jul)
APGAR_CUTOFF   <- 7          # exclude if Apgar < this value (apgar1 < 7 = 0-6)
EXCLUSION_MODE <- "apgar_only"   # "apgar_only" | "apgar_resustype" | "apgar_both"
# ESSC-only group: team agreed it is too small to model -- report it as counts in
# Tables 2/3 (Script 17) and OMIT its odds ratio here. FALSE = omit from models.
ANALYSE_ESSC_ONLY <- FALSE
# Minimum events/observations before a model is fitted (stability guard).
MIN_EVENTS <- 5

cat("=============================================================\n")
cat("  Script 18: Outcome odds ratios (Table 4) + sensitivity\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat(sprintf("  Exclusion mode : %s\n", EXCLUSION_MODE))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD + DERIVE  (shared derivation, mirrors Scripts 16/17)
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
    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    apgar1_num = if_else(!is.na(apgar1_num) & (apgar1_num < 0 | apgar1_num > 10), NA_real_, apgar1_num),
    apgar5_num = suppressWarnings(as.numeric(apgar5)),
    apgar5_num = if_else(!is.na(apgar5_num) & (apgar5_num < 0 | apgar5_num > 10), NA_real_, apgar5_num),
    apgar_excl = if (APGAR_FIELD == "apgar5") apgar5_num else apgar1_num,

    sex = factor(case_when(
      toupper(trimws(as.character(gender))) %in% c("F","FEMALE") ~ "F",
      toupper(trimws(as.character(gender))) %in% c("M","MALE")   ~ "M", TRUE ~ NA_character_),
      levels = c("F","M")),

    # BVM/CPR ascertainment from both fields (resustype is post-Mar-2025 only;
    # resus is 100% complete and carries the same multi-select tokens).
    rt_bvmcpr = !is.na(resustype) & grepl("BVM|CPR", resustype),
    rs_bvmcpr = !is.na(resus)     & grepl("BVM|CPR", resus),
    excl_apgar = !is.na(apgar_excl) & apgar_excl < APGAR_CUTOFF,
    excl_resus = if (EXCLUSION_MODE == "apgar_only") FALSE
                 else if (EXCLUSION_MODE == "apgar_both") (rt_bvmcpr | rs_bvmcpr)
                 else rt_bvmcpr,
    excluded   = excl_apgar | excl_resus,

    di_known = !(is.na(delivinter) | trimws(delivinter) == "" | trimws(delivinter) == "U"),
    has_dcc  = di_known & grepl("DCC", delivinter),
    has_s2s  = di_known & grepl("S2S", delivinter),
    exposure = factor(case_when(
      !di_known          ~ NA_character_,
      has_dcc & has_s2s  ~ "Both DCC and ESSC",
      has_dcc & !has_s2s ~ "DCC only",
      !has_dcc & has_s2s ~ "ESSC only",
      TRUE               ~ "Neither"),
      levels = c("Neither","DCC only","ESSC only","Both DCC and ESSC")),

    nnd = neotreeoutcome == "NND")

pop <- df %>%
  filter(facility == "SMCH", inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission),
         matched == TRUE,                       # discharge outcome required for NND
         adm_date >= INT_START, adm_date <= INT_END,
         excluded == FALSE, di_known == TRUE)

cat(sprintf("Analysis set (eligible, known exposure, matched): %d\n", nrow(pop)))
cat(sprintf("  deaths before discharge (NND): %d\n\n", sum(pop$nnd, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 2. MODEL HELPER -- crude + adjusted(BW) + adjusted(BW+sex) on one sample
# -----------------------------------------------------------------------------
or_ci <- function(fit) {
  # Wald 95% CI via confint.default -- no MASS dependency (avoids masking
  # dplyr::select) and runs anywhere. Difference from profile CIs is negligible
  # at these sample sizes. Swap to confint(fit) if profile CIs are preferred.
  est <- coef(fit)
  ci  <- confint.default(fit)
  tibble(term = names(est), OR = exp(est),
         CI_lower = exp(ci[, 1]), CI_upper = exp(ci[, 2]))
}

run_stratum <- function(data, stratum_label) {
  # complete-case on the fullest model's covariates so all 3 models share the sample
  md <- data %>% select(nnd, exposure, birthweight_num, sex) %>% drop_na()
  md$exposure <- droplevels(md$exposure)

  n_tot <- nrow(md); n_death <- sum(md$nnd)
  cat(sprintf("\n--- %s ---\n", stratum_label))
  cat(sprintf("  Complete cases: %d, deaths: %d\n", n_tot, n_death))
  grp_counts <- md %>% group_by(exposure) %>%
    summarise(n = n(), deaths = sum(nnd), rate = round(100*sum(nnd)/n(), 1), .groups = "drop") %>%
    mutate(stratum = stratum_label)
  print(as.data.frame(grp_counts), row.names = FALSE)

  # Team decision (21 Jul): ESSC-only is too small to model. Keep it in the counts
  # above (and in Tables 2/3), but drop it from the OR models. Removing a single
  # non-reference category does not change the other groups' ORs vs Neither.
  if (!ANALYSE_ESSC_ONLY) {
    md <- md %>% filter(exposure != "ESSC only")
    md$exposure <- droplevels(md$exposure)
    n_tot <- nrow(md); n_death <- sum(md$nnd)
    cat(sprintf("  (ESSC-only dropped from models; modelled n=%d, deaths=%d)\n",
                n_tot, n_death))
  }

  # guards: need Neither reference present, >=2 exposure levels, enough events
  if (nlevels(md$exposure) < 2 || n_death < MIN_EVENTS ||
      !("Neither" %in% levels(md$exposure))) {
    cat("  >> Too few events / groups for a stable model. Skipped.\n")
    return(list(full = tibble(), counts = grp_counts))
  }

  fits <- list(
    unadjusted        = glm(nnd ~ exposure, data = md, family = binomial()),
    `adjusted (BW)`   = glm(nnd ~ exposure + birthweight_num, data = md, family = binomial()),
    `adjusted (BW+sex)` = glm(nnd ~ exposure + birthweight_num + sex, data = md, family = binomial()))

  full <- imap_dfr(fits, function(fit, mlabel) {
    or_ci(fit) %>% filter(grepl("^exposure", term)) %>%
      mutate(model = mlabel, stratum = stratum_label,
             comparison = sub("^exposure", "", term),
             OR = round(OR, 3), CI_lower = round(CI_lower, 3), CI_upper = round(CI_upper, 3)) %>%
      select(stratum, model, comparison, OR, CI_lower, CI_upper)
  })

  cat("  OR vs Neither (exposure terms):\n")
  full %>% mutate(line = sprintf("    %-18s %-18s OR %.2f (%.2f-%.2f)",
                                 model, comparison, OR, CI_lower, CI_upper)) %>%
    pull(line) %>% cat(sep = "\n"); cat("\n")

  list(full = full, counts = grp_counts)
}

# -----------------------------------------------------------------------------
# 3. RUN THE THREE STRATA
# -----------------------------------------------------------------------------
strata <- list(
  "All eligible infants" = pop,
  "LBW <2500 g"          = pop %>% filter(!is.na(birthweight_num) & birthweight_num < 2500),
  "Sensitivity 1500-2499 g" = pop %>% filter(!is.na(birthweight_num) &
                                               birthweight_num >= 1500 & birthweight_num <= 2499))
res <- imap(strata, ~run_stratum(.x, .y))

full_tbl   <- bind_rows(map(res, "full"))
counts_tbl <- bind_rows(map(res, "counts"))
write_csv(full_tbl,   file.path(OUTPUT_DIR, "18a_or_full_models.csv"))
write_csv(counts_tbl, file.path(OUTPUT_DIR, "18c_or_counts.csv"))

# -----------------------------------------------------------------------------
# 4. TABLE 4 LAYOUT -- OR (95% CI) string, group x model, one block per stratum
# -----------------------------------------------------------------------------
table4 <- full_tbl %>%
  mutate(cell = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower, CI_upper)) %>%
  select(stratum, comparison, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell) %>%
  # ensure Neither reference row appears explicitly
  bind_rows(distinct(., stratum) %>%
              mutate(comparison = "Neither DCC nor ESSC",
                     unadjusted = "Reference",
                     `adjusted (BW)` = "Reference",
                     `adjusted (BW+sex)` = "Reference")) %>%
  mutate(comparison = factor(comparison,
           levels = c("Neither DCC nor ESSC","DCC only","ESSC only","Both DCC and ESSC"))) %>%
  arrange(stratum, comparison)

cat("\n================= TABLE 4: Odds ratios (vs Neither) =================\n")
print(as.data.frame(table4), row.names = FALSE)
write_csv(table4, file.path(OUTPUT_DIR, "18b_or_table4.csv"))

cat("\nSaved: 18a_or_full_models.csv, 18b_or_table4.csv, 18c_or_counts.csv\n")
cat("\n=== Script 18 complete ===\n")
cat("Reminder: report these as supportive observational evidence (indication\n")
cat("bias caveat), alongside the descriptive Tables 1-3.\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
