# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 04c: ITS with Seasonal Adjustment (Fourier Terms)
# =============================================================================
# Purpose:
#   Re-run the primary ITS analysis from Script 04b but with harmonic
#   (Fourier) seasonal terms added to the regression model. This controls
#   for the annual temperature cycle, allowing b2 and b3 to estimate the
#   intervention effect over and above expected seasonal variation.
#
# Motivation:
#   Scripts 04 and 04b produce a unadjusted ITS. The post-intervention period
#   spans multiple seasons (Nov 2024 - Mar 2026), so the slope estimate b3
#   reflects both the intervention trend and seasonal oscillation. Adding
#   sin/cos Fourier terms partitions variance correctly: the model captures
#   the seasonal rhythm explicitly, and b2/b3 become seasonally-adjusted
#   estimates. This is the methodologically appropriate sensitivity analysis
#   for any ITS with a seasonal outcome.
#
# Seasonal model:
#   pct_normal ~ time + intervention + time_after + sin_t + cos_t
#   where sin_t = sin(2 * pi * month_num / 12)
#         cos_t = cos(2 * pi * month_num / 12)
#   month_num = 1 (January) ... 12 (December)
#   This Fourier pair captures one full annual cycle using only 2 df.
#
# Unadjusted model (from 04b, reproduced here for direct comparison):
#   pct_normal ~ time + intervention + time_after
#
# Analysis windows (identical to 04b):
#   Primary     : January 2023 - March 2026 (39 months; 22 pre, 17 post)
#   Sensitivity : January 2022 - March 2026 (51 months; 34 pre, 17 post)
#
# Outputs (all to 03-OUTPUTS/, all prefixed 04c_):
#   04c_its_seasonality_correction_log.txt         -- full run log
#   04c_its_seasonality_correction_coefficients.csv -- adj + unadj coefficients
#   04c_its_seasonality_correction_monthly.csv      -- monthly data + fitted
#   04c_its_seasonality_correction_plot.png         -- two-panel figure
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)
library(lmtest)
library(sandwich)

args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

DATA_DIR    <- file.path(SCRIPT_DIR, "..", "00-DATA")
OUTPUTS_DIR <- file.path(SCRIPT_DIR, "..", "03-OUTPUTS")
dir.create(OUTPUTS_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUTPUTS_DIR, "04c_its_seasonality_correction_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
log <- function(...) {
  msg <- paste0(...)
  message(msg)
  writeLines(msg, log_con)
}

log("=============================================================")
log("Script 04c: ITS with Seasonal Adjustment (Fourier Terms)")
log("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("=============================================================")


# -----------------------------------------------------------------------------
# 1. PARAMETERS
# -----------------------------------------------------------------------------

INT_START       <- as.Date("2024-11-01")
INT_END         <- as.Date("2026-03-31")
HIST_START      <- as.Date("2023-01-01")   # primary window
HIST_START_SENS <- as.Date("2022-01-01")   # sensitivity window
NORMAL_LO       <- 36.5
NORMAL_HI       <- 37.5

COL_HIST <- "#4A90C4"
COL_INT  <- "#E07B32"
COL_CF   <- "#7B9E87"
COL_SADJ <- "#8B5CF6"    # purple for seasonally-adjusted points/line


# -----------------------------------------------------------------------------
# 2. LOAD DATA
# -----------------------------------------------------------------------------

data_files <- list.files(DATA_DIR, pattern = "ZIM_db_master.*\\.csv$", full.names = TRUE)
if (length(data_files) == 0) stop("No master data file found in 01-DATA/")
DATA_FILE <- data_files[length(data_files)]
log("Data file: ", basename(DATA_FILE))

df_raw <- read_csv(DATA_FILE, show_col_types = FALSE)
log("Rows loaded: ", nrow(df_raw))


# -----------------------------------------------------------------------------
# 3. DERIVE VARIABLES AND BUILD MONTHLY SERIES
# -----------------------------------------------------------------------------

df <- df_raw %>%
  mutate(
    adm_dt   = ymd_hms(datetimeadmission, quiet = TRUE),
    admdate  = as.Date(adm_dt),
    inborn   = inorout %in% c("Yes", "true", "True", "In"),
    is_readmission = readmission == "Y",
    temperature_num = suppressWarnings(as.numeric(temperature)),
    temp_normal = case_when(
      is.na(temperature_num)                                        ~ NA,
      temperature_num >= NORMAL_LO & temperature_num <= NORMAL_HI  ~ TRUE,
      TRUE                                                          ~ FALSE
    ),
    year_month = as.Date(floor_date(adm_dt, "month"))
  )

smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission)
  )

log("SMCH inborn non-readmission records: ", nrow(smch_inborn))

monthly <- smch_inborn %>%
  group_by(year_month) %>%
  summarise(
    n_total   = n(),
    n_temp    = sum(!is.na(temperature_num)),
    n_normal  = sum(temp_normal == TRUE, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    pct_normal = if_else(n_temp > 0, round(100 * n_normal / n_temp, 2), NA_real_)
  ) %>%
  arrange(year_month)

log("Monthly rows (full history): ", nrow(monthly))


# -----------------------------------------------------------------------------
# 4. ITS BUILDER FUNCTION WITH FOURIER TERMS
# -----------------------------------------------------------------------------

fit_its_fourier <- function(data, label) {
  log("")
  log("--- ", label, " ---")

  d <- data %>%
    filter(!is.na(pct_normal)) %>%
    mutate(
      month_num    = month(year_month),
      sin_t        = sin(2 * pi * month_num / 12),
      cos_t        = cos(2 * pi * month_num / 12),
      time         = as.numeric(year_month - min(year_month)) / 30.4375,
      intervention = as.integer(year_month >= INT_START),
      time_after   = pmax(0, as.numeric(year_month - INT_START) / 30.4375)
    )

  n_pre  <- sum(d$year_month < INT_START)
  n_post <- sum(d$year_month >= INT_START)
  log("Months: ", nrow(d), " (pre=", n_pre, ", post=", n_post, ")")

  # -- Unadjusted model (no seasonal terms) ----------------------------------
  fit_unadj <- lm(pct_normal ~ time + intervention + time_after,
                  data = d, weights = n_temp)

  # -- Seasonally adjusted model (Fourier terms) -----------------------------
  fit_adj <- lm(pct_normal ~ time + intervention + time_after + sin_t + cos_t,
                data = d, weights = n_temp)

  # -- Ljung-Box test on adjusted residuals -----------------------------------
  lb <- Box.test(residuals(fit_adj), lag = 1, type = "Ljung-Box")
  log("Ljung-Box p (adjusted model): ", round(lb$p.value, 4))
  use_nw <- lb$p.value < 0.05
  if (use_nw) {
    log("Autocorrelation detected -- using Newey-West HAC SEs (lag=2)")
    se_adj   <- coeftest(fit_adj,   vcov = NeweyWest(fit_adj,   lag = 2, prewhite = FALSE))
    se_unadj <- coeftest(fit_unadj, vcov = NeweyWest(fit_unadj, lag = 2, prewhite = FALSE))
  } else {
    log("No autocorrelation detected -- using OLS SEs")
    se_adj   <- coeftest(fit_adj)
    se_unadj <- coeftest(fit_unadj)
  }

  # -- Extract coefficient tables -------------------------------------------
  tidy_coef <- function(ct, model_label) {
    data.frame(
      term      = rownames(ct),
      estimate  = ct[, 1],
      std_error = ct[, 2],
      t_value   = ct[, 3],
      p_value   = ct[, 4],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        model    = model_label,
        window   = label,
        se_type  = if_else(use_nw, "Newey-West", "OLS"),
        sig      = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          p_value < 0.10  ~ ".",
          TRUE            ~ ""
        )
      )
  }

  coefs_adj   <- tidy_coef(se_adj,   "seasonally_adjusted")
  coefs_unadj <- tidy_coef(se_unadj, "unadjusted")

  for (nm in c("b1_time", "b2_intervention", "b3_time_after")) {
    row_adj   <- coefs_adj[coefs_adj$term %in% c("time","intervention","time_after"), ]
    row_unadj <- coefs_unadj[coefs_unadj$term %in% c("time","intervention","time_after"), ]
  }

  log("Adjusted model coefficients:")
  for (i in seq_len(nrow(coefs_adj))) {
    log(sprintf("  %-16s  est=%8.4f  se=%7.4f  p=%6.4f %s",
      coefs_adj$term[i], coefs_adj$estimate[i],
      coefs_adj$std_error[i], coefs_adj$p_value[i], coefs_adj$sig[i]))
  }
  log("Unadjusted model (b2, b3 for comparison):")
  b2u <- coefs_unadj[coefs_unadj$term == "intervention", ]
  b3u <- coefs_unadj[coefs_unadj$term == "time_after", ]
  log(sprintf("  b2 (unadj): est=%8.4f  p=%6.4f %s", b2u$estimate, b2u$p_value, b2u$sig))
  log(sprintf("  b3 (unadj): est=%8.4f  p=%6.4f %s", b3u$estimate, b3u$p_value, b3u$sig))

  # -- Fitted values and seasonally-adjusted series --------------------------
  d <- d %>%
    mutate(
      fitted_adj   = fitted(fit_adj),
      fitted_unadj = fitted(fit_unadj),
      # Counterfactual: remove b2 and b3 from adjusted model
      b2_hat = coef(fit_adj)["intervention"],
      b3_hat = coef(fit_adj)["time_after"],
      counterfactual = fitted_adj - b2_hat * intervention - b3_hat * time_after,
      # Seasonal component: sin_t and cos_t contribution
      seasonal_component = coef(fit_adj)["sin_t"] * sin_t +
                           coef(fit_adj)["cos_t"] * cos_t,
      # Seasonally-adjusted observed values
      # (remove seasonal component, re-center at series mean)
      pct_normal_adj = pct_normal - seasonal_component +
                       mean(seasonal_component, na.rm = TRUE)
    )

  list(
    data        = d,
    coefs_adj   = coefs_adj,
    coefs_unadj = coefs_unadj,
    fit_adj     = fit_adj,
    fit_unadj   = fit_unadj,
    use_nw      = use_nw
  )
}


# -----------------------------------------------------------------------------
# 5. RUN PRIMARY AND SENSITIVITY WINDOWS
# -----------------------------------------------------------------------------

monthly_full <- monthly %>%
  filter(year_month >= HIST_START_SENS & year_month <= INT_END)

monthly_primary <- monthly_full %>%
  filter(year_month >= HIST_START)

log("")
log("=== PRIMARY WINDOW (Jan 2023 - Mar 2026) ===")
res_primary <- fit_its_fourier(monthly_primary, "primary_Jan2023")

log("")
log("=== SENSITIVITY WINDOW (Jan 2022 - Mar 2026) ===")
res_sens <- fit_its_fourier(monthly_full, "sensitivity_Jan2022")


# -----------------------------------------------------------------------------
# 6. SAVE COEFFICIENTS CSV
# -----------------------------------------------------------------------------

coef_out <- bind_rows(
  res_primary$coefs_adj,
  res_primary$coefs_unadj,
  res_sens$coefs_adj,
  res_sens$coefs_unadj
) %>%
  select(window, model, se_type, term, estimate, std_error, t_value, p_value, sig)

write_csv(coef_out,
  file.path(OUTPUTS_DIR, "04c_its_seasonality_correction_coefficients.csv"))
log("")
log("Coefficients saved: 04c_its_seasonality_correction_coefficients.csv")


# -----------------------------------------------------------------------------
# 7. SAVE MONTHLY DATA CSV
# -----------------------------------------------------------------------------

monthly_out <- res_primary$data %>%
  select(year_month, n_total, n_temp, n_normal, pct_normal,
         month_num, sin_t, cos_t, time, intervention, time_after,
         fitted_adj, fitted_unadj, counterfactual, seasonal_component,
         pct_normal_adj)

write_csv(monthly_out,
  file.path(OUTPUTS_DIR, "04c_its_seasonality_correction_monthly.csv"))
log("Monthly data saved: 04c_its_seasonality_correction_monthly.csv")


# -----------------------------------------------------------------------------
# 8. PLOT: two-panel figure
# -----------------------------------------------------------------------------

plot_its_fourier <- function(res, title_suffix) {

  d <- res$data

  # --- Panel A: raw data + adjusted model fit ---
  pA <- ggplot(d, aes(x = year_month)) +
    # Background shading
    annotate("rect",
      xmin = min(d$year_month[d$year_month < INT_START]),
      xmax = INT_START,
      ymin = -Inf, ymax = Inf,
      fill = COL_HIST, alpha = 0.06) +
    annotate("rect",
      xmin = INT_START,
      xmax = max(d$year_month) + days(31),
      ymin = -Inf, ymax = Inf,
      fill = COL_INT, alpha = 0.06) +
    # Intervention line
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    # Raw data points
    geom_point(aes(y = pct_normal, size = n_temp,
                   colour = if_else(year_month >= INT_START, "Intervention", "Historical")),
               alpha = 0.7) +
    # Adjusted fitted line
    geom_line(aes(y = fitted_adj), colour = COL_SADJ, linewidth = 1.1, linetype = "solid") +
    # Counterfactual
    geom_line(aes(y = counterfactual), colour = COL_CF, linewidth = 0.9, linetype = "dashed") +
    scale_colour_manual(values = c("Historical" = COL_HIST, "Intervention" = COL_INT),
                        name = NULL) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    scale_y_continuous(limits = c(55, 100), labels = function(x) paste0(x, "%")) +
    labs(
      title = paste0("Panel A: Raw data with seasonally-adjusted model fit"),
      subtitle = paste0(title_suffix,
        " | Purple line = adjusted ITS fit (includes Fourier terms)",
        " | Dashed green = counterfactual"),
      x = NULL, y = "% Normal admission temperature"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )

  # --- Panel B: seasonally-adjusted data + linear ITS ---
  # Extract linear ITS coefficients from adjusted model
  cf <- coef(res$fit_adj)
  b2 <- cf["intervention"]
  b3 <- cf["time_after"]

  d2 <- d %>%
    mutate(
      linear_fit_pre  = cf["(Intercept)"] + cf["time"] * time,
      linear_fit_post = cf["(Intercept)"] + cf["time"] * time +
                        b2 * intervention + b3 * time_after,
      linear_cf       = cf["(Intercept)"] + cf["time"] * time
    )

  # Extract b2/b3 p-values for annotation
  coef_table <- if (res$use_nw) {
    coeftest(res$fit_adj, vcov = NeweyWest(res$fit_adj, lag = 2, prewhite = FALSE))
  } else {
    coeftest(res$fit_adj)
  }
  b2_p <- coef_table["intervention", "Pr(>|t|)"]
  b3_p <- coef_table["time_after",   "Pr(>|t|)"]

  fmt_p <- function(p) {
    if (p < 0.001) return("p<0.001***")
    if (p < 0.01)  return(paste0("p=", sprintf("%.3f", p), "**"))
    if (p < 0.05)  return(paste0("p=", sprintf("%.3f", p), "*"))
    paste0("p=", sprintf("%.3f", p))
  }

  annot <- paste0(
    "b2 (level change) = ", sprintf("%+.2f", b2), " pp  ", fmt_p(b2_p), "\n",
    "b3 (slope change) = ", sprintf("%+.3f", b3), " pp/mo  ", fmt_p(b3_p)
  )

  pB <- ggplot(d2, aes(x = year_month)) +
    annotate("rect",
      xmin = min(d2$year_month[d2$year_month < INT_START]),
      xmax = INT_START,
      ymin = -Inf, ymax = Inf,
      fill = COL_HIST, alpha = 0.06) +
    annotate("rect",
      xmin = INT_START,
      xmax = max(d2$year_month) + days(31),
      ymin = -Inf, ymax = Inf,
      fill = COL_INT, alpha = 0.06) +
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    # Seasonally-adjusted data points
    geom_point(aes(y = pct_normal_adj, size = n_temp,
                   colour = if_else(year_month >= INT_START, "Intervention", "Historical")),
               alpha = 0.7) +
    # Linear ITS trend lines (pre and post)
    geom_line(data = d2 %>% filter(year_month < INT_START),
              aes(y = linear_fit_pre),  colour = COL_HIST, linewidth = 1.1) +
    geom_line(data = d2 %>% filter(year_month >= INT_START),
              aes(y = linear_fit_post), colour = COL_INT,  linewidth = 1.1) +
    # Counterfactual (linear)
    geom_line(aes(y = linear_cf), colour = COL_CF, linewidth = 0.9, linetype = "dashed") +
    # Annotation box
    annotate("label",
      x = INT_START + days(90), y = 60,
      label = annot, hjust = 0, size = 3.2,
      fill = "white", colour = "grey20") +
    scale_colour_manual(values = c("Historical" = COL_HIST, "Intervention" = COL_INT),
                        name = NULL) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = paste0("Panel B: Seasonally-adjusted data and estimated intervention effect"),
      subtitle = paste0(title_suffix,
        " | Points = observed % minus seasonal component",
        " | Dashed green = counterfactual trend"),
      x = NULL, y = "Seasonally-adjusted % normal temperature"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "grey40")
    )

  list(pA = pA, pB = pB)
}

log("")
log("Generating plots ...")

panels_primary <- plot_its_fourier(res_primary,
  "Primary window: Jan 2023 - Mar 2026 (22 pre, 17 post months)")

panels_sens <- plot_its_fourier(res_sens,
  "Sensitivity window: Jan 2022 - Mar 2026 (34 pre, 17 post months)")

# Combine: 2 rows x 2 cols (primary top, sensitivity bottom)
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  combined <- gridExtra::arrangeGrob(
    panels_primary$pA, panels_primary$pB,
    panels_sens$pA,    panels_sens$pB,
    ncol = 2,
    top = "Script 04c: ITS with Fourier Seasonal Adjustment -- SMCH Golden Hour"
  )
  ggsave(
    file.path(OUTPUTS_DIR, "04c_its_seasonality_correction_plot.png"),
    combined, width = 16, height = 12, dpi = 150
  )
  log("4-panel plot saved (gridExtra layout)")
} else {
  # Fallback: save primary panels only, stacked
  library(patchwork)
  combined <- (panels_primary$pA / panels_primary$pB) +
    plot_annotation(
      title = "Script 04c: ITS with Fourier Seasonal Adjustment -- SMCH Golden Hour",
      subtitle = "Primary window: Jan 2023 - Mar 2026 | Sensitivity omitted (install gridExtra for 4-panel)"
    )
  ggsave(
    file.path(OUTPUTS_DIR, "04c_its_seasonality_correction_plot.png"),
    combined, width = 12, height = 12, dpi = 150
  )
  log("2-panel plot saved (patchwork layout; primary window only)")
}


# -----------------------------------------------------------------------------
# 9. SUMMARY COMPARISON TABLE
# -----------------------------------------------------------------------------

log("")
log("=============================================================")
log("COMPARISON: UNADJUSTED vs SEASONALLY-ADJUSTED (primary window)")
log("=============================================================")

get_b <- function(coef_df, term_name) {
  row <- coef_df[coef_df$term == term_name, ]
  if (nrow(row) == 0) return(list(est = NA, p = NA, sig = ""))
  list(est = row$estimate, p = row$p_value, sig = row$sig)
}

for (term in c("intervention", "time_after")) {
  b_adj   <- get_b(res_primary$coefs_adj,   term)
  b_unadj <- get_b(res_primary$coefs_unadj, term)
  label   <- if (term == "intervention") "b2 (level change)" else "b3 (slope change)"
  log(sprintf("%-22s  Unadjusted: %+7.3f (p=%5.3f%s)  |  Adjusted: %+7.3f (p=%5.3f%s)",
    label,
    b_unadj$est, b_unadj$p, b_unadj$sig,
    b_adj$est,   b_adj$p,   b_adj$sig))
}

log("")
log("Interpretation guide:")
log("  If adjusted b3 ~ unadjusted b3 (and both non-significant): seasonal")
log("    variation is not masking or creating a spurious intervention effect.")
log("  If adjusted b3 becomes significant and negative: the unadjusted model")
log("    underestimated the intervention benefit due to seasonal noise.")
log("  If adjusted b3 becomes non-significant when unadjusted was significant:")
log("    the previous result was driven by seasonal confounding.")
log("")
log("Script 04c complete.")
close(log_con)
