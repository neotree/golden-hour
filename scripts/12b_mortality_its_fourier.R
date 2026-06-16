# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 12b: Fourier-Adjusted ITS -- Neonatal Mortality (Seasonality Correction)
# =============================================================================
# Purpose:
#   Formally test whether seasonal variation in neonatal mortality confounds
#   the ITS results from Script 12.  Adds a Fourier harmonic pair (sin_t,
#   cos_t) capturing the annual seasonal cycle to the segmented regression.
#   The intervention coefficients b2 and b3 then estimate the effect above
#   and beyond the seasonal pattern.
#
#   This script is the mortality analogue of Script 04c (Fourier-adjusted
#   temperature ITS) and follows exactly the same modelling strategy.
#
# Model (seasonally adjusted):
#   pct_nnd ~ time + intervention + time_after + sin_t + cos_t
#   where sin_t = sin(2 * pi * month_num / 12)
#         cos_t = cos(2 * pi * month_num / 12)
#   Weighted by monthly denominator (DC + NND).
#   Ljung-Box autocorrelation check; Newey-West HAC SEs (lag=2) if p < 0.05.
#
# Subgroups (same as Script 12):
#   (1) All SMCH inborn (primary)
#   (2) LBW: birthweight < 2500 g
#   (3) Normal BW: birthweight >= 2500 g
#
# Output filenames deliberately distinct from Script 12 outputs to avoid
# overwriting.  Outputs (all to 03-OUTPUTS/):
#   12b_mortality_fourier_log.txt
#   12b_mortality_fourier_monthly.csv        -- monthly data with Fourier terms
#   12b_mortality_fourier_coefficients.csv   -- unadjusted vs adjusted coefs
#   12b_mortality_fourier_plot_all.png       -- all-inborn two-panel plot
#   12b_mortality_fourier_plot_lbw.png       -- LBW two-panel plot
#
# DSH note: ASCII-only (no non-ASCII characters).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)
library(lmtest)    # coeftest()
library(sandwich)  # NeweyWest()
library(patchwork) # pA / pB stacking

# -- Locate script directory ---------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag),
                                      mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

# -- Paths --------------------------------------------------------------------
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
LOG_FILE <- file.path(OUTPUT_DIR, "12b_mortality_fourier_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour -- Fourier-Adjusted Mortality ITS (Script 12b)\n")
cat(strrep("=", 60))
cat(sprintf("\nRun started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Study boundaries ---------------------------------------------------------
INT_START       <- as.Date("2024-11-01")
INT_END         <- as.Date("2026-03-31")
HIST_START      <- as.Date("2023-01-01")   # primary window
HIST_END        <- as.Date("2024-10-31")
HIST_START_SENS <- as.Date("2022-01-01")   # sensitivity window
MIN_N_DENOM     <- 10

# -- Colours ------------------------------------------------------------------
COL_HIST  <- "#4A90C4"
COL_INT   <- "#E07B32"
COL_ADJ   <- "#7B4EA8"   # purple for Fourier-adjusted fit
COL_CF    <- "#7B9E87"
FILL_HIST <- "#DDEEF8"
FILL_INT  <- "#FEF0DC"

cat("Configuration:\n")
cat(sprintf("  Data        : %s\n", DATA_PATH))
cat(sprintf("  Output dir  : %s\n", OUTPUT_DIR))
cat(sprintf("  Primary ITS : %s to %s\n", HIST_START, INT_END))
cat(sprintf("  Sensitivity : %s to %s\n", HIST_START_SENS, INT_END))
cat(sprintf("  Int start   : %s\n\n", INT_START))


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)
parse_probs <- problems(raw)
if (nrow(parse_probs) > 0) {
  cat(sprintf("  WARNING: %d parse problem(s).\n", nrow(parse_probs)))
} else {
  cat("  No parse problems.\n")
}
cat(sprintf("  Rows: %d\n\n", nrow(raw)))


# -----------------------------------------------------------------------------
# 2. CLEAN AND FILTER
# -----------------------------------------------------------------------------

df <- raw %>%
  mutate(
    adm_dt   = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date = as.Date(adm_dt),
    adm_month = as.Date(floor_date(adm_dt, "month")),

    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                             ~ NA
    ),
    is_readmission = readmission == "Y",
    matched        = match_type == "direct_match",

    # Outcome
    nnd        = neotreeoutcome == "NND",
    discharged = neotreeoutcome == "DC",
    has_outcome = neotreeoutcome %in% c("NND", "DC"),

    # Birthweight
    bw_num = suppressWarnings(as.numeric(birthweight)),
    lbw    = !is.na(bw_num) & bw_num < 2500,
    nbw    = !is.na(bw_num) & bw_num >= 2500,

    # Period
    period_primary = case_when(
      adm_date >= INT_START  & adm_date <= INT_END  ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END ~ "Historical",
      TRUE                                           ~ "Other"
    ),
    period_sens = case_when(
      adm_date >= INT_START       & adm_date <= INT_END       ~ "Intervention",
      adm_date >= HIST_START_SENS & adm_date <= HIST_END      ~ "Historical",
      TRUE                                                     ~ "Other"
    )
  )

smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission),
    matched        == TRUE
  )

cat(sprintf("SMCH inborn matched: %d records\n\n", nrow(smch_inborn)))


# -----------------------------------------------------------------------------
# 3. BUILD MONTHLY TABLE (with Fourier terms)
# -----------------------------------------------------------------------------

build_monthly <- function(data, period_var) {
  data %>%
    group_by(adm_month, period = .data[[period_var]]) %>%
    summarise(
      n_admissions = n(),
      n_outcome    = sum(has_outcome, na.rm = TRUE),
      n_nnd        = sum(nnd, na.rm = TRUE),
      n_dc         = sum(discharged, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    filter(period %in% c("Intervention", "Historical")) %>%
    arrange(adm_month) %>%
    mutate(
      pct_nnd   = if_else(n_outcome > 0,
                          round(100 * n_nnd / n_outcome, 2), NA_real_),
      # Fourier terms: month_num = calendar month (1-12)
      month_num = as.integer(format(adm_month, "%m")),
      sin_t     = sin(2 * pi * month_num / 12),
      cos_t     = cos(2 * pi * month_num / 12)
    )
}

# Primary window
q_primary <- smch_inborn %>% filter(period_primary %in% c("Intervention", "Historical"))
monthly_all  <- build_monthly(q_primary, "period_primary")
monthly_lbw  <- build_monthly(q_primary %>% filter(lbw), "period_primary")
monthly_nbw  <- build_monthly(q_primary %>% filter(nbw), "period_primary")

# Sensitivity window
q_sens       <- smch_inborn %>% filter(period_sens %in% c("Intervention", "Historical"))
monthly_sens <- build_monthly(q_sens, "period_sens")

cat(sprintf("Monthly rows (all inborn, primary)  : %d (%d pre, %d int)\n",
            nrow(monthly_all),
            sum(monthly_all$period == "Historical"),
            sum(monthly_all$period == "Intervention")))
cat(sprintf("Monthly rows (LBW, primary)         : %d\n", nrow(monthly_lbw)))
cat(sprintf("Monthly rows (NBW, primary)         : %d\n", nrow(monthly_nbw)))
cat(sprintf("Monthly rows (all inborn, sensitiv.): %d\n\n", nrow(monthly_sens)))

# Save monthly data
monthly_all  %>% mutate(subgroup = "All inborn") %>%
  bind_rows(
    monthly_lbw  %>% mutate(subgroup = "LBW (<2500g)"),
    monthly_nbw  %>% mutate(subgroup = "NBW (>=2500g)"),
    monthly_sens %>% mutate(subgroup = "All inborn (sensitivity)")
  ) %>%
  select(subgroup, adm_month, period, n_outcome, n_nnd, pct_nnd,
         month_num, sin_t, cos_t) %>%
  write_csv(file.path(OUTPUT_DIR, "12b_mortality_fourier_monthly.csv"))
cat("Saved: 12b_mortality_fourier_monthly.csv\n\n")


# -----------------------------------------------------------------------------
# 4. ITS FIT FUNCTION (unadjusted) -- replicates Script 12 without Fourier
# -----------------------------------------------------------------------------

fit_unadj <- function(monthly, label) {
  d <- monthly %>%
    filter(n_outcome >= MIN_N_DENOM) %>%
    arrange(adm_month) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )
  if (sum(d$intervention == 0) < 3 || sum(d$intervention == 1) < 2) return(NULL)

  model <- lm(pct_nnd ~ time + intervention + time_after,
              data = d, weights = n_outcome)
  lb    <- Box.test(residuals(model), lag = 1, type = "Ljung-Box")

  # Always report NW as well for comparison
  nw_test <- coeftest(model, vcov = NeweyWest(model, lag = 2, prewhite = FALSE))

  extract <- function(mat, se_type) {
    data.frame(
      analysis  = label,
      model_type = "Unadjusted",
      term      = rownames(mat),
      estimate  = mat[, 1],
      std_error = mat[, 2],
      t_value   = mat[, 3],
      p_value   = mat[, 4],
      se_type   = se_type,
      stringsAsFactors = FALSE
    ) %>%
      mutate(ci_lower = estimate - 1.96 * std_error,
             ci_upper = estimate + 1.96 * std_error,
             across(where(is.numeric), ~ round(.x, 4)))
  }

  list(
    model    = model,
    tidy_ols = extract(summary(model)$coefficients, "OLS"),
    tidy_nw  = extract(as.matrix(nw_test), "Newey-West (lag=2)"),
    lb       = lb,
    d        = d %>% mutate(
      fitted_unadj   = fitted(model),
      counterfactual = coef(model)["(Intercept)"] + coef(model)["time"] * time
    )
  )
}


# -----------------------------------------------------------------------------
# 5. FOURIER ITS FIT FUNCTION
# -----------------------------------------------------------------------------

fit_fourier <- function(monthly, label) {

  d <- monthly %>%
    filter(n_outcome >= MIN_N_DENOM) %>%
    arrange(adm_month) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )

  n_pre  <- sum(d$intervention == 0)
  n_post <- sum(d$intervention == 1)

  cat(sprintf("\n%s\n", strrep("-", 64)))
  cat(sprintf("  Fourier ITS: %s\n", label))
  cat(sprintf("  Pre-intervention months : %d\n", n_pre))
  cat(sprintf("  Intervention months     : %d\n", n_post))
  cat(sprintf("%s\n", strrep("-", 64)))

  if (n_pre < 3 || n_post < 2) {
    cat("  WARNING: Too few months. Skipping.\n")
    return(NULL)
  }

  # Seasonally-adjusted model
  model_adj <- lm(
    pct_nnd ~ time + intervention + time_after + sin_t + cos_t,
    data    = d,
    weights = n_outcome
  )

  cat("\n--- Fourier-adjusted OLS summary ---\n")
  print(summary(model_adj))

  # Ljung-Box on adjusted residuals
  lb_adj <- Box.test(residuals(model_adj), lag = 1, type = "Ljung-Box")
  use_nw <- lb_adj$p.value < 0.05
  cat(sprintf("\n--- Ljung-Box (adjusted residuals, lag 1) ---\n"))
  cat(sprintf("  Q = %.4f  df = %d  p = %.4f\n",
              lb_adj$statistic, lb_adj$parameter, lb_adj$p.value))
  if (use_nw) {
    cat("  >> Significant autocorrelation. Newey-West SEs applied.\n")
  } else {
    cat("  >> No significant autocorrelation. OLS SEs reliable.\n")
  }

  nw_adj <- coeftest(model_adj,
                     vcov = NeweyWest(model_adj, lag = 2, prewhite = FALSE))
  cat("\n--- Newey-West coefficients (Fourier model, lag=2) ---\n")
  print(nw_adj)

  extract <- function(mat, se_type) {
    data.frame(
      analysis   = label,
      model_type = "Fourier-adjusted",
      term       = rownames(mat),
      estimate   = mat[, 1],
      std_error  = mat[, 2],
      t_value    = mat[, 3],
      p_value    = mat[, 4],
      se_type    = se_type,
      stringsAsFactors = FALSE
    ) %>%
      mutate(ci_lower = estimate - 1.96 * std_error,
             ci_upper = estimate + 1.96 * std_error,
             across(where(is.numeric), ~ round(.x, 4)))
  }

  tidy_ols <- extract(summary(model_adj)$coefficients, "OLS")
  tidy_nw  <- extract(as.matrix(nw_adj), "Newey-West (lag=2)")

  # Preferred coefs: NW if autocorrelation, OLS otherwise
  preferred <- if (use_nw) tidy_nw else tidy_ols

  # Key coefficients
  cf     <- coef(model_adj)
  b0_adj <- cf["(Intercept)"]
  b1_adj <- cf["time"]
  b2_adj <- cf["intervention"]
  b3_adj <- cf["time_after"]

  cat("\n--- Interpretation (Fourier-adjusted) ---\n")
  cat(sprintf("  b1 (pre-trend)    : %.3f pp/month\n", b1_adj))
  cat(sprintf("  b2 (level change) : %.3f pp  (p = %.4f)\n", b2_adj,
              preferred %>% filter(term == "intervention") %>% pull(p_value)))
  cat(sprintf("  b3 (slope change) : %.3f pp/month  (p = %.4f)\n", b3_adj,
              preferred %>% filter(term == "time_after") %>% pull(p_value)))
  cat(sprintf("  sin_t coef        : %.3f  (p = %.4f)\n", cf["sin_t"],
              preferred %>% filter(term == "sin_t") %>% pull(p_value)))
  cat(sprintf("  cos_t coef        : %.3f  (p = %.4f)\n", cf["cos_t"],
              preferred %>% filter(term == "cos_t") %>% pull(p_value)))

  # Fitted values from Fourier model on full range
  ci_fit <- predict(model_adj, interval = "confidence", level = 0.95)

  # Counterfactual: Fourier model without intervention terms
  d_cf <- d %>%
    mutate(intervention = 0L, time_after = 0L)
  cf_pred <- predict(model_adj, newdata = d_cf)

  # Seasonally-adjusted outcome: subtract Fourier contribution
  # seasonal_component = sin_coef * sin_t + cos_coef * cos_t
  sin_coef <- cf["sin_t"]
  cos_coef <- cf["cos_t"]
  d <- d %>%
    mutate(
      seasonal_component = sin_coef * sin_t + cos_coef * cos_t,
      pct_nnd_adj        = pct_nnd - seasonal_component,
      fitted_adj         = ci_fit[, "fit"],
      fit_lower          = ci_fit[, "lwr"],
      fit_upper          = ci_fit[, "upr"],
      counterfactual_adj = cf_pred
    )

  list(
    model    = model_adj,
    tidy_ols = tidy_ols,
    tidy_nw  = tidy_nw,
    preferred = preferred,
    lb_adj   = lb_adj,
    use_nw   = use_nw,
    d        = d,
    b0 = b0_adj, b1 = b1_adj, b2 = b2_adj, b3 = b3_adj,
    sin_coef = sin_coef, cos_coef = cos_coef
  )
}


# -----------------------------------------------------------------------------
# 6. FIT ALL MODELS
# -----------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 64))
cat("\n  FITTING FOURIER-ADJUSTED ITS MODELS\n")
cat(strrep("=", 64))

labels <- list(
  all  = "All inborn -- primary (Jan 2023 - Mar 2026)",
  lbw  = "LBW (<2500g) -- primary (Jan 2023 - Mar 2026)",
  nbw  = "NBW (>=2500g) -- primary (Jan 2023 - Mar 2026)",
  sens = "All inborn -- sensitivity (Jan 2022 - Mar 2026)"
)

unadj_all  <- fit_unadj(monthly_all,  labels$all)
unadj_lbw  <- fit_unadj(monthly_lbw,  labels$lbw)
unadj_nbw  <- fit_unadj(monthly_nbw,  labels$nbw)
unadj_sens <- fit_unadj(monthly_sens, labels$sens)

adj_all  <- fit_fourier(monthly_all,  labels$all)
adj_lbw  <- fit_fourier(monthly_lbw,  labels$lbw)
adj_nbw  <- fit_fourier(monthly_nbw,  labels$nbw)
adj_sens <- fit_fourier(monthly_sens, labels$sens)


# -----------------------------------------------------------------------------
# 7. SAVE COEFFICIENT COMPARISON TABLE
# -----------------------------------------------------------------------------

# Use preferred SEs for each model (NW if autocorrelation, else OLS)
preferred_se <- function(res_u, res_a) {
  u <- if (!is.null(res_u)) {
    if (res_u$lb$p.value < 0.05) res_u$tidy_nw else res_u$tidy_ols
  } else NULL

  a <- if (!is.null(res_a)) {
    if (res_a$use_nw) res_a$tidy_nw else res_a$tidy_ols
  } else NULL

  bind_rows(u, a)
}

all_coefs <- bind_rows(
  preferred_se(unadj_all,  adj_all),
  preferred_se(unadj_lbw,  adj_lbw),
  preferred_se(unadj_nbw,  adj_nbw),
  preferred_se(unadj_sens, adj_sens)
) %>%
  filter(term %in% c("time", "intervention", "time_after",
                     "sin_t", "cos_t")) %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE            ~ ""
    )
  )

write_csv(all_coefs, file.path(OUTPUT_DIR, "12b_mortality_fourier_coefficients.csv"))
cat("\nSaved: 12b_mortality_fourier_coefficients.csv\n")


# -----------------------------------------------------------------------------
# 8. RESULTS SUMMARY TABLE
# -----------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 64))
cat("\n  COMPARISON: UNADJUSTED vs FOURIER-ADJUSTED (preferred SEs)\n")
cat(strrep("=", 64))
cat("\n")

print_comparison <- function(res_u, res_a, label) {
  cat(sprintf("\n--- %s ---\n", label))

  if (is.null(res_u) || is.null(res_a)) {
    cat("  One or both models not fitted.\n")
    return(invisible(NULL))
  }

  pref_u <- if (res_u$lb$p.value < 0.05) res_u$tidy_nw else res_u$tidy_ols
  pref_a <- if (res_a$use_nw) res_a$tidy_nw else res_a$tidy_ols

  se_note_u <- if (res_u$lb$p.value < 0.05) "Newey-West" else "OLS"
  se_note_a <- if (res_a$use_nw) "Newey-West" else "OLS"

  cat(sprintf("  Unadjusted SEs  : %s  (Ljung-Box p = %.4f)\n",
              se_note_u, res_u$lb$p.value))
  cat(sprintf("  Fourier adj SEs : %s  (Ljung-Box p = %.4f)\n",
              se_note_a, res_a$lb_adj$p.value))
  cat("\n")

  for (trm in c("intervention", "time_after")) {
    row_u <- pref_u %>% filter(term == trm)
    row_a <- pref_a %>% filter(term == trm)
    lbl <- if (trm == "intervention") "b2 (level change)" else "b3 (slope change/mo)"
    fmt_row <- function(r) {
      if (nrow(r) == 0) return("  --")
      pstar <- if (r$p_value < 0.001) "<0.001 ***" else
               if (r$p_value < 0.01)  sprintf("%.4f **", r$p_value) else
               if (r$p_value < 0.05)  sprintf("%.4f *",  r$p_value) else
               sprintf("%.4f", r$p_value)
      sprintf("%.3f pp  [%.3f, %.3f]  p = %s",
              r$estimate, r$ci_lower, r$ci_upper, pstar)
    }
    cat(sprintf("  %-22s  Unadjusted : %s\n", lbl, fmt_row(row_u)))
    cat(sprintf("  %-22s  Adjusted   : %s\n", "",  fmt_row(row_a)))
    cat("\n")
  }

  # Fourier term significance (proxy for seasonal signal strength)
  sin_row <- pref_a %>% filter(term == "sin_t")
  cos_row <- pref_a %>% filter(term == "cos_t")
  if (nrow(sin_row) > 0 && nrow(cos_row) > 0) {
    cat(sprintf("  sin_t : %.3f (p = %.3f)   cos_t : %.3f (p = %.3f)\n",
                sin_row$estimate, sin_row$p_value,
                cos_row$estimate, cos_row$p_value))
  }
}

print_comparison(unadj_all,  adj_all,  labels$all)
print_comparison(unadj_lbw,  adj_lbw,  labels$lbw)
print_comparison(unadj_nbw,  adj_nbw,  labels$nbw)
print_comparison(unadj_sens, adj_sens, labels$sens)


# -----------------------------------------------------------------------------
# 9. PLOT FUNCTION
# -----------------------------------------------------------------------------

make_fourier_plot <- function(res_u, res_a, subtitle_a, subtitle_b) {

  if (is.null(res_u) || is.null(res_a)) {
    cat("  Skipping plot: model(s) not fitted.\n")
    return(NULL)
  }

  d_u <- res_u$d
  d_a <- res_a$d

  # ---- Panel A: raw data with Fourier-adjusted ITS fit --------------------
  d_pre_u  <- d_u %>% filter(intervention == 0)
  d_post_u <- d_u %>% filter(intervention == 1)
  d_pre_a  <- d_a %>% filter(intervention == 0)
  d_post_a <- d_a %>% filter(intervention == 1)

  # Annotation text: b2 and b3 unadjusted vs adjusted
  pref_u <- if (res_u$lb$p.value < 0.05) res_u$tidy_nw else res_u$tidy_ols
  pref_a <- if (res_a$use_nw) res_a$tidy_nw else res_a$tidy_ols

  b2_u <- pref_u %>% filter(term == "intervention")
  b3_u <- pref_u %>% filter(term == "time_after")
  b2_a <- pref_a %>% filter(term == "intervention")
  b3_a <- pref_a %>% filter(term == "time_after")

  annot_a <- sprintf(
    "Fourier-adjusted ITS:\nb2 = %+.3f pp  (p = %.3f)\nb3 = %+.3f pp/mo  (p = %.3f)",
    b2_a$estimate, b2_a$p_value, b3_a$estimate, b3_a$p_value
  )

  pA <- ggplot(d_a, aes(x = adm_month)) +
    annotate("rect",
             xmin = min(d_a$adm_month), xmax = INT_START,
             ymin = -Inf, ymax = Inf, fill = FILL_HIST, alpha = 1) +
    annotate("rect",
             xmin = INT_START, xmax = max(d_a$adm_month) + days(31),
             ymin = -Inf, ymax = Inf, fill = FILL_INT, alpha = 1) +
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey45", linewidth = 0.55) +

    # Unadjusted fitted lines (thin, grey)
    geom_line(data = d_pre_u,  aes(y = fitted_unadj),
              colour = "grey65", linewidth = 0.8, linetype = "solid") +
    geom_line(data = d_post_u, aes(y = fitted_unadj),
              colour = "grey65", linewidth = 0.8, linetype = "solid") +

    # Fourier-adjusted fit ribbons
    geom_ribbon(data = d_pre_a,
                aes(ymin = fit_lower, ymax = fit_upper),
                fill = COL_ADJ, alpha = 0.15) +
    geom_ribbon(data = d_post_a,
                aes(ymin = fit_lower, ymax = fit_upper),
                fill = COL_ADJ, alpha = 0.15) +

    # Counterfactual (Fourier-adjusted)
    geom_line(data = d_post_a,
              aes(y = counterfactual_adj),
              colour = COL_CF, linetype = "dotted", linewidth = 1.0) +

    # Fourier-adjusted fitted lines
    geom_line(data = d_pre_a,  aes(y = fitted_adj),
              colour = COL_ADJ, linewidth = 1.3) +
    geom_line(data = d_post_a, aes(y = fitted_adj),
              colour = COL_ADJ, linewidth = 1.3) +

    # Observed points
    geom_point(aes(y = pct_nnd, colour = period, size = n_outcome), alpha = 0.80) +

    annotate("label",
             x = INT_START + days(90),
             y = max(d_a$pct_nnd, na.rm = TRUE) * 0.82,
             label = annot_a, hjust = 0, size = 3.0,
             fill = "white", colour = "grey20") +

    scale_colour_manual(
      name   = "Period",
      values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
    ) +
    scale_size_continuous(name = "N (DC+NND)", range = c(1.5, 5.5),
                          breaks = c(50, 100, 200, 400)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y",
                 expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0.05, 0.12))) +
    labs(
      title    = "Panel A: Raw monthly NND rate with Fourier-adjusted ITS fit (purple)",
      subtitle = subtitle_a,
      x = NULL, y = "% NND (of DC + NND)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title     = element_text(face = "bold", size = 10.5),
      plot.subtitle  = element_text(colour = "grey40", size = 8.5),
      axis.text.x    = element_text(size = 8),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )

  # ---- Panel B: seasonally-adjusted NND rate with linear ITS lines --------
  pref_a_b <- pref_a  # same for the linear fits

  annot_b <- sprintf(
    "Adjusted for seasonal cycle (sin + cos).\nb2 = %+.3f pp  (p = %.3f)\nb3 = %+.3f pp/mo  (p = %.3f)\nCompare unadjusted:\nb2 = %+.3f pp  (p = %.3f)\nb3 = %+.3f pp/mo  (p = %.3f)",
    b2_a$estimate, b2_a$p_value, b3_a$estimate, b3_a$p_value,
    b2_u$estimate, b2_u$p_value, b3_u$estimate, b3_u$p_value
  )

  # Linear fits on adjusted outcome (remove Fourier component from fitted)
  cf_adj <- coef(res_a$model)
  d_adj_pre  <- d_a %>% filter(intervention == 0) %>%
    mutate(
      fit_linear_adj = cf_adj["(Intercept)"] + cf_adj["time"] * time
    )
  d_adj_post <- d_a %>% filter(intervention == 1) %>%
    mutate(
      fit_linear_adj = cf_adj["(Intercept)"] +
        cf_adj["time"] * time +
        cf_adj["intervention"] +
        cf_adj["time_after"] * time_after
    )

  pB <- ggplot(d_a, aes(x = adm_month)) +
    annotate("rect",
             xmin = min(d_a$adm_month), xmax = INT_START,
             ymin = -Inf, ymax = Inf, fill = FILL_HIST, alpha = 1) +
    annotate("rect",
             xmin = INT_START, xmax = max(d_a$adm_month) + days(31),
             ymin = -Inf, ymax = Inf, fill = FILL_INT, alpha = 1) +
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey45", linewidth = 0.55) +

    # Seasonally-adjusted points
    geom_point(aes(y = pct_nnd_adj, colour = period, size = n_outcome), alpha = 0.75) +

    # Linear ITS lines on adjusted data
    geom_line(data = d_adj_pre,  aes(y = fit_linear_adj),
              colour = COL_HIST, linewidth = 1.3) +
    geom_line(data = d_adj_post, aes(y = fit_linear_adj),
              colour = COL_INT,  linewidth = 1.3) +

    annotate("label",
             x = INT_START + days(90),
             y = max(d_a$pct_nnd_adj, na.rm = TRUE) * 0.80,
             label = annot_b, hjust = 0, size = 2.8,
             fill = "white", colour = "grey20") +

    scale_colour_manual(
      name   = "Period",
      values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
    ) +
    scale_size_continuous(name = "N (DC+NND)", range = c(1.5, 5.5),
                          breaks = c(50, 100, 200, 400)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y",
                 expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0.05, 0.12))) +
    labs(
      title    = "Panel B: Seasonally-adjusted NND rate with linear ITS lines",
      subtitle = subtitle_b,
      x = NULL, y = "% NND (seasonal component removed)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title     = element_text(face = "bold", size = 10.5),
      plot.subtitle  = element_text(colour = "grey40", size = 8.5),
      axis.text.x    = element_text(size = 8),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )

  pA / pB
}


# -----------------------------------------------------------------------------
# 10. GENERATE PLOTS
# -----------------------------------------------------------------------------

cat("\nGenerating plots...\n")

plot_all <- make_fourier_plot(
  unadj_all, adj_all,
  subtitle_a = paste0("All SMCH inborn (matched records), Jan 2023 - Mar 2026. ",
                      "Grey lines = unadjusted ITS; purple = Fourier-adjusted."),
  subtitle_b = paste0("Seasonal component (sin + cos) removed from observed rates. ",
                      "Linear ITS lines overlaid. Outcome: pct_nnd = NND / (DC+NND) * 100.")
)
if (!is.null(plot_all)) {
  ggsave(file.path(OUTPUT_DIR, "12b_mortality_fourier_plot_all.png"),
         plot = plot_all, width = 13, height = 9, dpi = 300, units = "in")
  cat("Saved: 12b_mortality_fourier_plot_all.png\n")
}

plot_lbw <- make_fourier_plot(
  unadj_lbw, adj_lbw,
  subtitle_a = paste0("LBW subgroup (birthweight < 2500 g), Jan 2023 - Mar 2026. ",
                      "Grey lines = unadjusted ITS; purple = Fourier-adjusted."),
  subtitle_b = paste0("Seasonal component removed. ",
                      "Key question: does seasonal variation explain the LBW b3 slope (p=0.043)?")
)
if (!is.null(plot_lbw)) {
  ggsave(file.path(OUTPUT_DIR, "12b_mortality_fourier_plot_lbw.png"),
         plot = plot_lbw, width = 13, height = 9, dpi = 300, units = "in")
  cat("Saved: 12b_mortality_fourier_plot_lbw.png\n")
}


# -----------------------------------------------------------------------------
# 11. FINAL SUMMARY
# -----------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 64))
cat("\n  SCRIPT 12b SUMMARY\n")
cat(strrep("=", 64))
cat("\n")
cat("Model: pct_nnd ~ time + intervention + time_after + sin_t + cos_t\n")
cat("       Weighted OLS. NW HAC SEs (lag=2) if Ljung-Box p < 0.05.\n\n")

coef_summary <- all_coefs %>%
  filter(term %in% c("intervention", "time_after")) %>%
  select(analysis, model_type, term, estimate, ci_lower, ci_upper, p_value, sig, se_type) %>%
  arrange(analysis, model_type, term)

print(coef_summary)

cat(sprintf("\n\n=== Script 12b complete ===\n"))
cat(sprintf("Run completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Outputs in: %s\n", OUTPUT_DIR))

# -- Close log ----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
