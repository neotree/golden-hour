# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 12: Interrupted Time Series (ITS) -- Neonatal Mortality
# =============================================================================
# Purpose:
#   Evaluate whether the Golden Hour intervention (November 2024) was
#   associated with a change in the monthly neonatal death rate (NND rate)
#   at SMCH, using the same segmented regression ITS framework as Script 04b.
#
# Outcome:
#   Monthly NND rate = NND / (NND + DC) * 100
#
#   Denominator: babies with a definitive discharge outcome (DC = discharged
#   alive, NND = neonatal death). Excluded from denominator:
#     - NaN / missing neotreeoutcome (not yet discharged or not recorded)
#     - TRO  (transferred out -- outcome unknown)
#     - TRH  (transferred to higher level -- outcome unknown)
#     - DAMA (discharged against medical advice -- outcome uncertain)
#     - ABS  (absconded -- outcome unknown)
#   These exclusions are applied consistently in both periods.
#
# Data:
#   Master NeoTree file (ZIM_db_master_joined_to_20260401.csv).
#   Filter to: facility == "SMCH", inborn == TRUE, match_type == "direct_match"
#   (direct_match ensures the discharge outcome neotreeoutcome is populated;
#   unmatched admissions do not have discharge data and cannot contribute to
#   the NND denominator).
#
# Analysis window (PRIMARY):
#   Pre-intervention : January 2023 - October 2024  (22 months)
#   Intervention     : November 2024 - March 2026   (17 months)
#   Total            : 39 months
#
#   Note: 2022 excluded from primary analysis consistent with Script 04b,
#   to avoid the anomalous temperature/data period. A sensitivity analysis
#   including 2022 is also run (see Section 7).
#
# ITS model (identical to Script 04b):
#   pct_nnd ~ time + intervention + time_after
#   Weighted by monthly denominator (DC + NND).
#   Ljung-Box autocorrelation check; Newey-West HAC SEs (lag=2) if p < 0.05.
#
# Subgroups:
#   (1) All SMCH inborn (primary)
#   (2) LBW: birthweight < 2500 g
#   (3) Normal BW: birthweight >= 2500 g
#   Birthweight from the NeoTree admission variable 'birthweight'.
#
# Outputs (all to 03-OUTPUTS/):
#   12_mortality_its_log.txt
#   12a_mortality_monthly.csv          -- monthly mortality rates (all + subgroups)
#   12b_its_coefficients.csv           -- ITS coefficients, all models
#   12c_mortality_its_plot.png         -- overall ITS plot
#   12d_lbw_mortality_its_plot.png     -- LBW subgroup ITS plot
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)
library(lmtest)    # coeftest()
library(sandwich)  # NeweyWest()

# -- Locate script directory ---------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag),
                                      mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

# -- Path configuration -------------------------------------------------------
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

# -- Log file ------------------------------------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "12_mortality_its_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour Mortality ITS Analysis -- SMCH Zimbabwe\n")
cat(sprintf("Run started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2026-03-31")
HIST_START <- as.Date("2023-01-01")   # primary: excludes 2022
HIST_END   <- as.Date("2024-10-31")
HIST_START_SENS <- as.Date("2022-01-01")  # sensitivity: includes 2022

# -- Minimum denominator per month to include in ITS --------------------------
MIN_N_DENOM <- 10   # months with fewer than 10 DC+NND excluded

# -- Colours -------------------------------------------------------------------
COL_HIST  <- "#4A90C4"
COL_INT   <- "#E07B32"
FILL_HIST <- "#DDEEF8"
FILL_INT  <- "#FEF0DC"
COL_CF    <- "#7B9E87"

cat("Configuration:\n")
cat(sprintf("  Data path  : %s\n", DATA_PATH))
cat(sprintf("  Output dir : %s\n", OUTPUT_DIR))
cat(sprintf("  ITS window : %s to %s (primary)\n", HIST_START, INT_END))
cat(sprintf("  ITS window : %s to %s (sensitivity)\n", HIST_START_SENS, INT_END))
cat(sprintf("  Int start  : %s\n", INT_START))


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("\nLoading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)

parse_probs <- problems(raw)
if (nrow(parse_probs) > 0) {
  cat(sprintf("  WARNING: %d parse problem(s) detected.\n", nrow(parse_probs)))
} else {
  cat("  No parse problems detected.\n")
}
cat(sprintf("  Rows loaded: %d\n", nrow(raw)))


# -----------------------------------------------------------------------------
# 2. CLEANING AND DERIVED VARIABLES
# -----------------------------------------------------------------------------

df <- raw %>%
  mutate(
    adm_dt    = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date  = as.Date(adm_dt),
    adm_month = floor_date(adm_dt, "month"),

    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                             ~ NA
    ),

    is_readmission = readmission == "Y",
    matched        = match_type == "direct_match",

    # -- Outcome derivation ---------------------------------------------------
    # neotreeoutcome: DC, NND, TRO, ABS, DAMA, TRH, NA
    # Include in denominator only DC and NND (definitive outcomes)
    nnd          = neotreeoutcome == "NND",
    discharged   = neotreeoutcome == "DC",
    has_outcome  = neotreeoutcome %in% c("NND", "DC"),  # denominator

    # -- Birthweight ----------------------------------------------------------
    bw_num = suppressWarnings(as.numeric(birthweight)),
    lbw    = !is.na(bw_num) & bw_num < 2500,
    nbw    = !is.na(bw_num) & bw_num >= 2500,

    # -- Period flag ----------------------------------------------------------
    period_primary = case_when(
      adm_date >= INT_START  & adm_date <= INT_END   ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END  ~ "Historical",
      TRUE                                            ~ "Other"
    ),
    period_sens = case_when(
      adm_date >= INT_START       & adm_date <= INT_END       ~ "Intervention",
      adm_date >= HIST_START_SENS & adm_date <= HIST_END      ~ "Historical",
      TRUE                                                     ~ "Other"
    )
  )


# -----------------------------------------------------------------------------
# 3. ANALYSIS POPULATION
# -----------------------------------------------------------------------------

# Primary population: SMCH, inborn, non-readmission, matched (has discharge data)
smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission),
    matched        == TRUE
  )

cat(sprintf("\nSMCH inborn matched records (all time): %d\n", nrow(smch_inborn)))

# Primary analysis window
q12_primary <- smch_inborn %>%
  filter(period_primary %in% c("Intervention", "Historical"))

cat(sprintf("Primary window (Jan 2023 - Mar 2026): %d records\n", nrow(q12_primary)))
cat(sprintf("  Pre-intervention : %d\n",
            sum(q12_primary$period_primary == "Historical")))
cat(sprintf("  Intervention     : %d\n",
            sum(q12_primary$period_primary == "Intervention")))

# Outcome summary
cat(sprintf("\nOutcome distribution (primary window):\n"))
print(table(outcome = q12_primary$neotreeoutcome, useNA = "ifany"))
cat(sprintf("\nHas definitive outcome (DC or NND): %d / %d (%.1f%%)\n",
            sum(q12_primary$has_outcome, na.rm = TRUE),
            nrow(q12_primary),
            100 * mean(q12_primary$has_outcome, na.rm = TRUE)))
cat(sprintf("NND count : %d\n", sum(q12_primary$nnd, na.rm = TRUE)))
cat(sprintf("DC  count : %d\n", sum(q12_primary$discharged, na.rm = TRUE)))
cat(sprintf("Overall NND rate (of DC+NND): %.1f%%\n",
            100 * sum(q12_primary$nnd, na.rm = TRUE) /
            sum(q12_primary$has_outcome, na.rm = TRUE)))

# Birthweight note
cat(sprintf("\nBirthweight recorded: %d / %d (%.1f%%)\n",
            sum(!is.na(q12_primary$bw_num)),
            nrow(q12_primary),
            100 * mean(!is.na(q12_primary$bw_num))))
cat(sprintf("LBW (<2500g): %d (%.1f%% of those with BW)\n",
            sum(q12_primary$lbw, na.rm = TRUE),
            100 * sum(q12_primary$lbw, na.rm = TRUE) /
            sum(!is.na(q12_primary$bw_num))))


# -----------------------------------------------------------------------------
# 4. BUILD MONTHLY TABLES
# -----------------------------------------------------------------------------

build_monthly <- function(data, period_var) {
  data %>%
    group_by(adm_month, period = .data[[period_var]]) %>%
    summarise(
      n_admissions = n(),
      n_dc         = sum(discharged, na.rm = TRUE),
      n_nnd        = sum(nnd,        na.rm = TRUE),
      n_outcome    = sum(has_outcome, na.rm = TRUE),  # DC + NND (denominator)
      n_tro_dama   = sum(neotreeoutcome %in% c("TRO","TRH","DAMA","ABS"),
                         na.rm = TRUE),
      n_missing    = sum(is.na(neotreeoutcome)),
      pct_nnd      = if_else(
        n_outcome > 0,
        round(100 * n_nnd / n_outcome, 2),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    filter(period %in% c("Intervention", "Historical")) %>%
    arrange(adm_month) %>%
    mutate(adm_month_date = as.Date(adm_month))
}

monthly_all <- build_monthly(q12_primary, "period_primary")
cat(sprintf("\nMonthly data (all inborn): %d months (%d pre, %d int)\n",
            nrow(monthly_all),
            sum(monthly_all$period == "Historical"),
            sum(monthly_all$period == "Intervention")))

# LBW subgroup
q12_lbw <- q12_primary %>% filter(lbw == TRUE)
monthly_lbw <- build_monthly(q12_lbw, "period_primary")
cat(sprintf("Monthly data (LBW): %d months\n", nrow(monthly_lbw)))

# Normal BW subgroup
q12_nbw <- q12_primary %>% filter(nbw == TRUE)
monthly_nbw <- build_monthly(q12_nbw, "period_primary")
cat(sprintf("Monthly data (NBW): %d months\n", nrow(monthly_nbw)))

# Monthly data printout
cat("\n--- Monthly NND rates (all inborn, primary window) ---\n")
print(monthly_all %>%
        select(adm_month_date, period, n_admissions, n_outcome,
               n_nnd, pct_nnd, n_missing), n = Inf)

# Combined monthly output
monthly_combined <- monthly_all %>%
  mutate(subgroup = "All inborn") %>%
  bind_rows(
    monthly_lbw %>% mutate(subgroup = "LBW (<2500g)"),
    monthly_nbw %>% mutate(subgroup = "NBW (>=2500g)")
  ) %>%
  select(subgroup, adm_month_date, period, n_admissions, n_outcome,
         n_nnd, n_dc, n_tro_dama, n_missing, pct_nnd)

write_csv(monthly_combined,
          file.path(OUTPUT_DIR, "12a_mortality_monthly.csv"))
cat("\nSaved: 12a_mortality_monthly.csv\n")


# -----------------------------------------------------------------------------
# 5. ITS MODEL FUNCTION
# -----------------------------------------------------------------------------

# Fits the segmented regression ITS, runs Ljung-Box, applies Newey-West if
# autocorrelation detected. Returns a list with model object, tidy coefficient
# tables (OLS and NW), Ljung-Box result, and annotated monthly data.

fit_its <- function(monthly, label, min_n = MIN_N_DENOM) {

  # Exclude months with too few denominator observations
  d <- monthly %>%
    filter(n_outcome >= min_n) %>%
    arrange(adm_month) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )

  n_pre  <- sum(d$intervention == 0)
  n_post <- sum(d$intervention == 1)
  cat(sprintf("\n%s\n", strrep("-", 60)))
  cat(sprintf("  ITS: %s\n", label))
  cat(sprintf("  Pre-intervention months : %d\n", n_pre))
  cat(sprintf("  Intervention months     : %d\n", n_post))
  cat(sprintf("  Months excluded (n<%d)  : %d\n",
              min_n, nrow(monthly) - nrow(d)))
  cat(sprintf("%s\n", strrep("-", 60)))

  if (n_pre < 3 || n_post < 2) {
    cat("  WARNING: Too few months for reliable ITS. Skipping.\n")
    return(NULL)
  }

  # Weighted OLS
  model <- lm(
    pct_nnd ~ time + intervention + time_after,
    data    = d,
    weights = n_outcome
  )

  cat("\n--- OLS summary ---\n")
  print(summary(model))

  # Ljung-Box autocorrelation test (lag 1)
  lb <- Box.test(residuals(model), lag = 1, type = "Ljung-Box")
  cat(sprintf("\n--- Ljung-Box test (lag 1) ---\n"))
  cat(sprintf("  Q = %.4f  df = %d  p = %.4f\n",
              lb$statistic, lb$parameter, lb$p.value))
  if (lb$p.value < 0.05) {
    cat("  >> Significant autocorrelation. Newey-West SEs applied.\n")
  } else {
    cat("  >> No significant autocorrelation. OLS SEs reliable.\n")
  }

  # Newey-West HAC SEs
  nw_test <- coeftest(model, vcov = NeweyWest(model, lag = 2, prewhite = FALSE))
  cat("\n--- Newey-West coefficients (lag=2) ---\n")
  print(nw_test)

  # Tidy coefficient tables
  extract_coefs <- function(mat, se_type) {
    data.frame(
      analysis  = label,
      term      = rownames(mat),
      estimate  = mat[, 1],
      std_error = mat[, 2],
      t_value   = mat[, 3],
      p_value   = mat[, 4],
      se_type   = se_type,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        ci_lower = estimate - 1.96 * std_error,
        ci_upper = estimate + 1.96 * std_error,
        across(where(is.numeric), ~ round(.x, 4))
      ) %>%
      select(analysis, term, estimate, std_error, ci_lower, ci_upper,
             t_value, p_value, se_type)
  }

  tidy_ols <- extract_coefs(summary(model)$coefficients, "OLS")
  tidy_nw  <- extract_coefs(as.matrix(nw_test), "Newey-West (lag=2)")

  # Interpretation
  cf_vec <- coef(model)
  b1 <- cf_vec["time"];         b2 <- cf_vec["intervention"]
  b3 <- cf_vec["time_after"];   b0 <- cf_vec["(Intercept)"]

  cat("\n--- Interpretation ---\n")
  cat(sprintf("  Pre-intervention trend (b1)    : %.3f %%-pts/month (%.2f /year)\n",
              b1, b1 * 12))
  cat(sprintf("  Level change at Nov 2024 (b2)  : %.3f %%-pts (immediate step)\n", b2))
  cat(sprintf("  Slope change post-int  (b3)    : %.3f %%-pts/month\n", b3))
  cat(sprintf("  Post-int slope (b1+b3)         : %.3f %%-pts/month\n", b1 + b3))

  # Fitted values + counterfactual
  ci_fit <- predict(model, interval = "confidence", level = 0.95)
  vcv    <- vcov(model)
  cf_se  <- sqrt(
    vcv["(Intercept)", "(Intercept)"] +
    d$time^2 * vcv["time", "time"] +
    2 * d$time * vcv["(Intercept)", "time"]
  )

  d <- d %>%
    mutate(
      fitted_its     = fitted(model),
      counterfactual = b0 + b1 * time,
      fit_lower      = ci_fit[, "lwr"],
      fit_upper      = ci_fit[, "upr"],
      cf_lower       = counterfactual - 1.96 * cf_se,
      cf_upper       = counterfactual + 1.96 * cf_se
    )

  list(
    model    = model,
    tidy_ols = tidy_ols,
    tidy_nw  = tidy_nw,
    lb       = lb,
    d        = d,
    b0 = b0, b1 = b1, b2 = b2, b3 = b3
  )
}


# -----------------------------------------------------------------------------
# 6. FIT ITS MODELS
# -----------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 64))
cat("\n  FITTING ITS MODELS\n")
cat(strrep("=", 64))
cat("\n")

its_all <- fit_its(monthly_all, "All inborn -- primary (Jan 2023 - Mar 2026)")
its_lbw <- fit_its(monthly_lbw, "LBW (<2500g) -- primary (Jan 2023 - Mar 2026)")
its_nbw <- fit_its(monthly_nbw, "NBW (>=2500g) -- primary (Jan 2023 - Mar 2026)")

# Sensitivity analysis: include 2022
q12_sens <- smch_inborn %>%
  filter(period_sens %in% c("Intervention", "Historical"))
monthly_all_sens <- build_monthly(q12_sens, "period_sens")
its_all_sens <- fit_its(monthly_all_sens,
                        "All inborn -- sensitivity (Jan 2022 - Mar 2026)")


# -----------------------------------------------------------------------------
# 7. SAVE COEFFICIENT TABLE
# -----------------------------------------------------------------------------

all_coefs <- bind_rows(
  if (!is.null(its_all))      bind_rows(its_all$tidy_ols, its_all$tidy_nw),
  if (!is.null(its_lbw))      bind_rows(its_lbw$tidy_ols, its_lbw$tidy_nw),
  if (!is.null(its_nbw))      bind_rows(its_nbw$tidy_ols, its_nbw$tidy_nw),
  if (!is.null(its_all_sens)) bind_rows(its_all_sens$tidy_ols, its_all_sens$tidy_nw)
)

write_csv(all_coefs, file.path(OUTPUT_DIR, "12b_its_coefficients.csv"))
cat("\nSaved: 12b_its_coefficients.csv\n")


# -----------------------------------------------------------------------------
# 8. PLOT FUNCTION
# -----------------------------------------------------------------------------

plot_mortality_its <- function(res, subtitle_text, cf_x_date) {

  if (is.null(res)) {
    cat("  Skipping plot: model was not fitted.\n")
    return(NULL)
  }

  d      <- res$d
  lb_p   <- res$lb$p.value
  d_pre  <- d %>% filter(intervention == 0)
  d_post <- d %>% filter(intervention == 1)

  inf_note <- if (lb_p < 0.05) {
    "Newey-West HAC SEs (autocorrelation detected, Ljung-Box p < 0.05)."
  } else {
    "OLS SEs (no significant autocorrelation, Ljung-Box p >= 0.05)."
  }

  # Counterfactual annotation y position
  cf_y <- res$b0 + res$b1 * d_post$time[
    which.min(abs(d_post$adm_month_date - cf_x_date))
  ]

  p <- ggplot(d, aes(x = adm_month_date)) +

    # Period shading
    annotate("rect",
             xmin = min(d$adm_month_date), xmax = INT_START,
             ymin = -Inf, ymax = Inf, fill = FILL_HIST, alpha = 1) +
    annotate("rect",
             xmin = INT_START, xmax = max(d$adm_month_date) + days(31),
             ymin = -Inf, ymax = Inf, fill = FILL_INT,  alpha = 1) +

    # Intervention line
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey45", linewidth = 0.55) +

    # Counterfactual ribbon + line
    geom_ribbon(data = d_post,
                aes(ymin = cf_lower, ymax = cf_upper),
                fill = COL_CF, alpha = 0.18) +
    geom_line(data = d_post,
              aes(y = counterfactual),
              colour = COL_CF, linetype = "dotted", linewidth = 1.0) +

    # Fitted ribbons
    geom_ribbon(data = d_pre,
                aes(ymin = fit_lower, ymax = fit_upper),
                fill = COL_HIST, alpha = 0.20) +
    geom_ribbon(data = d_post,
                aes(ymin = fit_lower, ymax = fit_upper),
                fill = COL_INT, alpha = 0.20) +

    # Fitted lines
    geom_line(data = d_pre,  aes(y = fitted_its),
              colour = COL_HIST, linewidth = 1.2) +
    geom_line(data = d_post, aes(y = fitted_its),
              colour = COL_INT, linewidth = 1.2) +

    # Observed points
    geom_point(aes(y = pct_nnd, colour = period, size = n_outcome),
               alpha = 0.80) +

    # Period labels
    annotate("text",
             x      = HIST_START + (INT_START - HIST_START) / 2,
             y      = max(d$pct_nnd, na.rm = TRUE) * 1.05,
             label  = "Historical period",
             colour = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
    annotate("text",
             x      = INT_START + (max(d$adm_month_date) - INT_START) / 2,
             y      = max(d$pct_nnd, na.rm = TRUE) * 1.05,
             label  = "Intervention period",
             colour = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +

    # Counterfactual label
    annotate("text",
             x      = cf_x_date,
             y      = cf_y - (max(d$pct_nnd, na.rm = TRUE) * 0.06),
             label  = "Counterfactual\n(pre-trend projected)",
             colour = COL_CF, hjust = 0.5, size = 2.7, fontface = "italic") +

    scale_colour_manual(
      name   = "Period",
      values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
    ) +
    scale_size_continuous(
      name   = "N\n(DC+NND)",
      range  = c(1.5, 5.5),
      breaks = c(50, 100, 200, 400)
    ) +
    scale_x_date(
      date_breaks = "3 months",
      date_labels = "%b\n%Y",
      expand      = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    labs(
      title    = "ITS -- Neonatal mortality rate (NND/discharged), SMCH",
      subtitle = subtitle_text,
      x        = NULL,
      y        = "% NND (of DC + NND)",
      caption  = paste0(
        "Segmented linear regression: pct_nnd ~ time + intervention + time_after,\n",
        "weighted by monthly denominator (DC + NND). ",
        "Denominator excludes TRO, TRH, DAMA, ABS, and missing outcomes.\n",
        "Dashed vertical line = intervention start (November 2024). ",
        "Dotted line = counterfactual (pre-intervention trend projected).\n",
        "Inference: ", inf_note
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11.5),
      plot.subtitle    = element_text(colour = "grey40", size = 8.8),
      plot.caption     = element_text(colour = "grey50", size = 7.2, hjust = 0),
      axis.text.x      = element_text(size = 8),
      axis.text.y      = element_text(size = 9),
      axis.title.y     = element_text(size = 9.5),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.key.size  = unit(0.5, "cm"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )

  return(p)
}


# -----------------------------------------------------------------------------
# 9. GENERATE AND SAVE PLOTS
# -----------------------------------------------------------------------------

cat("\nGenerating plots...\n")

# Overall ITS plot
p_all <- plot_mortality_its(
  its_all,
  subtitle_text = paste0(
    "All SMCH inborn (matched records), Jan 2023 - Mar 2026. ",
    "22 pre-intervention + 17 intervention months.\n",
    "Denominator: babies with definitive discharge outcome (DC or NND)."
  ),
  cf_x_date = as.Date("2025-09-01")
)

if (!is.null(p_all)) {
  ggsave(
    file.path(OUTPUT_DIR, "12c_mortality_its_plot.png"),
    plot = p_all, width = 13, height = 6, dpi = 300, units = "in"
  )
  cat("Saved: 12c_mortality_its_plot.png\n")
}

# LBW ITS plot
p_lbw <- plot_mortality_its(
  its_lbw,
  subtitle_text = paste0(
    "LBW subgroup (birthweight < 2500 g), SMCH inborn, Jan 2023 - Mar 2026. ",
    "22 pre-intervention + 17 intervention months.\n",
    "Denominator: LBW babies with definitive discharge outcome (DC or NND)."
  ),
  cf_x_date = as.Date("2025-09-01")
)

if (!is.null(p_lbw)) {
  ggsave(
    file.path(OUTPUT_DIR, "12d_lbw_mortality_its_plot.png"),
    plot = p_lbw, width = 13, height = 6, dpi = 300, units = "in"
  )
  cat("Saved: 12d_lbw_mortality_its_plot.png\n")
}


# -----------------------------------------------------------------------------
# 10. RESULTS SUMMARY
# -----------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 64))
cat("\n  MORTALITY ITS RESULTS SUMMARY\n")
cat(strrep("=", 64))
cat("\n")

print_summary <- function(res, label) {

  if (is.null(res)) {
    cat(sprintf("\n%s: model not fitted (insufficient data).\n", label))
    return(invisible(NULL))
  }

  cat(sprintf("\n--- %s ---\n", label))

  use_nw  <- res$lb$p.value < 0.05
  se_note <- ifelse(use_nw,
                    "Newey-West HAC SEs (autocorrelation detected)",
                    "OLS SEs (no significant autocorrelation)")
  coefs   <- if (use_nw) res$tidy_nw else res$tidy_ols

  cat(sprintf("  Standard errors: %s\n", se_note))
  cat(sprintf("  Ljung-Box p     : %.4f\n\n", res$lb$p.value))

  term_labels <- c(
    "time"         = "Pre-intervention trend (%-pts/month)",
    "intervention" = "Level change at Nov 2024 (%-pts, step)",
    "time_after"   = "Slope change post-intervention (%-pts/month)"
  )

  for (trm in names(term_labels)) {
    row <- coefs %>% filter(term == trm)
    if (nrow(row) == 1) {
      pval_str <- if (row$p_value < 0.001) "<0.001" else
                  sprintf("%.4f", row$p_value)
      sig_flag <- if (row$p_value < 0.05) " *" else
                  if (row$p_value < 0.10) " (marginal)" else ""
      cat(sprintf("  %s\n    Est: %.3f  95%% CI: [%.3f, %.3f]  p = %s%s\n\n",
                  term_labels[trm],
                  row$estimate, row$ci_lower, row$ci_upper,
                  pval_str, sig_flag))
    }
  }
}

print_summary(its_all,      "All inborn (primary)")
print_summary(its_lbw,      "LBW (<2500g)")
print_summary(its_nbw,      "NBW (>=2500g)")
print_summary(its_all_sens, "All inborn (sensitivity: includes 2022)")

cat("\n=== Script 12 complete ===\n")
cat(sprintf("Run completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Close log -----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
