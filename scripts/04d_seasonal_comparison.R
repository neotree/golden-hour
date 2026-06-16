# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 04d: Same-Season Year-on-Year Temperature Comparison
# =============================================================================
# Purpose:
#   Compare normal admission temperature rates between the intervention period
#   and the same calendar months in pre-intervention years, removing seasonal
#   effects by design. For each post-intervention month (Nov 2024 - Mar 2026),
#   the comparison value is the mean of the same calendar month(s) observed in
#   the historical period (Jan 2023 - Oct 2024).
#
# Rationale:
#   The ITS model (Scripts 04, 04b) estimates a post-intervention slope (b3)
#   that mixes intervention effect with seasonal variation. Script 04c adds
#   Fourier terms to control for this statistically. Script 04d takes a
#   non-parametric, design-based approach: by matching on calendar month,
#   seasonal effects cancel out by construction. This is more transparent
#   and easier to communicate to clinical audiences.
#
# Method:
#   1. For each month in the intervention period (Nov 2024 - Mar 2026),
#      find all same-calendar-month observations in the primary historical
#      window (Jan 2023 - Oct 2024).
#   2. Compute: matched_diff = observed_int - historical_mean_same_month
#   3. Test: one-sample t-test on matched_diffs (H0: mean_diff = 0).
#   4. Season-level comparison (Zimbabwe seasons):
#        Warm/wet  : November - April
#        Cool/dry  : May - October
#   5. Repeat using full historical window (Jan 2022 - Oct 2024) as
#      sensitivity.
#
# Historical windows:
#   Primary     : January 2023 - October 2024
#   Sensitivity : January 2022 - October 2024
#
# Intervention period: November 2024 - March 2026 (17 months)
#
# Outputs (all to 03-OUTPUTS/, all prefixed 04d_):
#   04d_seasonal_comparison_log.txt       -- full run log
#   04d_seasonal_comparison_monthly.csv   -- month-matched comparison table
#   04d_seasonal_comparison_plot.png      -- three-panel figure
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)

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

LOG_FILE <- file.path(OUTPUTS_DIR, "04d_seasonal_comparison_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
log <- function(...) {
  msg <- paste0(...)
  message(msg)
  writeLines(msg, log_con)
}

log("=============================================================")
log("Script 04d: Same-Season Year-on-Year Temperature Comparison")
log("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("=============================================================")


# -----------------------------------------------------------------------------
# 1. PARAMETERS
# -----------------------------------------------------------------------------

INT_START       <- as.Date("2024-11-01")
INT_END         <- as.Date("2026-03-31")
HIST_START      <- as.Date("2023-01-01")   # primary historical window
HIST_START_SENS <- as.Date("2022-01-01")   # sensitivity historical window
NORMAL_LO       <- 36.5
NORMAL_HI       <- 37.5

# Zimbabwe seasons (Southern Hemisphere)
# Warm/wet: November - April (months 11, 12, 1, 2, 3, 4)
# Cool/dry: May - October    (months 5, 6, 7, 8, 9, 10)
WARM_MONTHS <- c(11, 12, 1, 2, 3, 4)
COOL_MONTHS <- c(5, 6, 7, 8, 9, 10)

COL_HIST <- "#4A90C4"
COL_INT  <- "#E07B32"
COL_DIFF <- "#C0392B"
COL_NULL <- "grey60"


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
# 3. BUILD MONTHLY SERIES
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
    n_total  = n(),
    n_temp   = sum(!is.na(temperature_num)),
    n_normal = sum(temp_normal == TRUE, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(
    pct_normal = if_else(n_temp > 0, round(100 * n_normal / n_temp, 2), NA_real_),
    cal_month  = month(year_month),
    cal_year   = year(year_month),
    season     = if_else(cal_month %in% WARM_MONTHS, "Warm (Nov-Apr)", "Cool (May-Oct)"),
    period     = case_when(
      year_month >= INT_START                                  ~ "Intervention",
      year_month >= HIST_START & year_month <  INT_START      ~ "Historical (primary)",
      year_month >= HIST_START_SENS & year_month < HIST_START ~ "Historical (sensitivity only)",
      TRUE                                                     ~ "Before analysis window"
    )
  ) %>%
  arrange(year_month) %>%
  filter(!is.na(pct_normal))

log("Monthly rows with temperature data: ", nrow(monthly))


# -----------------------------------------------------------------------------
# 4. MATCHED COMPARISON FUNCTION
# -----------------------------------------------------------------------------

run_comparison <- function(monthly_data, hist_window_start, window_label) {
  log("")
  log("--- ", window_label, " ---")

  historical <- monthly_data %>%
    filter(year_month >= hist_window_start & year_month < INT_START)

  intervention <- monthly_data %>%
    filter(year_month >= INT_START & year_month <= INT_END)

  log("Historical months available: ", nrow(historical))
  log("Intervention months: ", nrow(intervention))

  # Historical mean per calendar month (across all historical years)
  hist_by_month <- historical %>%
    group_by(cal_month) %>%
    summarise(
      hist_mean   = mean(pct_normal, na.rm = TRUE),
      hist_sd     = sd(pct_normal,   na.rm = TRUE),
      hist_n      = n(),
      hist_se     = if_else(n() > 1, sd(pct_normal, na.rm = TRUE) / sqrt(n()), NA_real_),
      hist_lo95   = hist_mean - qt(0.975, df = max(n() - 1, 1)) * hist_se,
      hist_hi95   = hist_mean + qt(0.975, df = max(n() - 1, 1)) * hist_se,
      .groups     = "drop"
    )

  # Match each intervention month to its historical counterpart
  matched <- intervention %>%
    left_join(hist_by_month, by = "cal_month") %>%
    mutate(
      diff      = pct_normal - hist_mean,
      diff_lo95 = pct_normal - hist_hi95,
      diff_hi95 = pct_normal - hist_lo95,
      window    = window_label
    )

  # -- One-sample t-test: H0: mean(diff) = 0 ---------------------------------
  diffs <- matched$diff[!is.na(matched$diff)]
  ttest <- t.test(diffs, mu = 0)

  log(sprintf("Mean difference (intervention - historical same month): %+.2f pp",
    mean(diffs)))
  log(sprintf("95%% CI: [%+.2f, %+.2f]",
    ttest$conf.int[1], ttest$conf.int[2]))
  log(sprintf("t(%d) = %.3f, p = %.4f%s",
    ttest$parameter, ttest$statistic, ttest$p.value,
    if (ttest$p.value < 0.05) " *" else ""))

  # -- Season-level comparison -----------------------------------------------
  log("")
  log("By season:")
  for (seas in c("Warm (Nov-Apr)", "Cool (May-Oct)")) {
    m_seas <- matched %>% filter(season == seas)
    h_seas <- historical %>% filter(season == seas)
    if (nrow(m_seas) == 0 || nrow(h_seas) == 0) {
      log(sprintf("  %-18s  no data", seas))
      next
    }
    d_seas <- m_seas$diff[!is.na(m_seas$diff)]
    if (length(d_seas) < 2) {
      log(sprintf("  %-18s  int_mean=%.1f%%  hist_mean=%.1f%%  diff=%+.1f pp (n<2, no test)",
        seas, mean(m_seas$pct_normal), mean(h_seas$pct_normal),
        mean(d_seas)))
      next
    }
    tt <- t.test(d_seas, mu = 0)
    log(sprintf("  %-18s  int_mean=%.1f%%  hist_mean=%.1f%%  diff=%+.1f pp  p=%.3f%s",
      seas,
      mean(m_seas$pct_normal),
      mean(h_seas$pct_normal),
      mean(d_seas),
      tt$p.value,
      if (tt$p.value < 0.05) " *" else ""))
  }

  # -- Month-by-month table --------------------------------------------------
  log("")
  log("Month-by-month:")
  for (i in seq_len(nrow(matched))) {
    r <- matched[i, ]
    log(sprintf("  %s  obs=%.1f%%  hist_mean=%.1f%%  diff=%+.1f pp",
      format(r$year_month, "%b %Y"),
      r$pct_normal,
      r$hist_mean,
      r$diff))
  }

  list(
    matched      = matched,
    historical   = historical,
    hist_by_month = hist_by_month,
    ttest        = ttest,
    window_label = window_label
  )
}


# -----------------------------------------------------------------------------
# 5. RUN PRIMARY AND SENSITIVITY
# -----------------------------------------------------------------------------

res_prim <- run_comparison(monthly, HIST_START,      "primary_Jan2023")
res_sens <- run_comparison(monthly, HIST_START_SENS,  "sensitivity_Jan2022")


# -----------------------------------------------------------------------------
# 6. SAVE CSV
# -----------------------------------------------------------------------------

monthly_out <- bind_rows(
  res_prim$matched %>%
    select(window, year_month, cal_month, cal_year, season, n_temp,
           pct_normal, hist_mean, hist_sd, hist_n, hist_se,
           hist_lo95, hist_hi95, diff, diff_lo95, diff_hi95),
  res_sens$matched %>%
    select(window, year_month, cal_month, cal_year, season, n_temp,
           pct_normal, hist_mean, hist_sd, hist_n, hist_se,
           hist_lo95, hist_hi95, diff, diff_lo95, diff_hi95)
)

write_csv(monthly_out,
  file.path(OUTPUTS_DIR, "04d_seasonal_comparison_monthly.csv"))
log("")
log("Monthly comparison saved: 04d_seasonal_comparison_monthly.csv")


# -----------------------------------------------------------------------------
# 7. PLOT
# -----------------------------------------------------------------------------

month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")

make_plots <- function(res) {

  d_int  <- res$matched
  d_hist <- res$hist_by_month %>%
    mutate(month_label = month_labels[cal_month])
  d_int2 <- d_int %>%
    mutate(
      month_label = format(year_month, "%b '%y"),
      month_short = format(year_month, "%b")
    )

  tt <- res$ttest
  test_annot <- sprintf(
    "Overall: mean diff = %+.1f pp\n95%% CI [%+.1f, %+.1f]\nt(%d) = %.2f, p = %.3f%s",
    mean(d_int$diff, na.rm = TRUE),
    tt$conf.int[1], tt$conf.int[2],
    tt$parameter, tt$statistic, tt$p.value,
    if (tt$p.value < 0.05) " *" else ""
  )

  # ---- Panel A: full time series with historical range by calendar month ----
  # Build reference band: for each intervention month, mark the historical range
  d_full <- bind_rows(
    res$historical %>%
      mutate(period2 = "Historical") %>%
      select(year_month, pct_normal, n_temp, season, period2),
    d_int %>%
      mutate(period2 = "Intervention") %>%
      select(year_month, pct_normal, n_temp, season, period2)
  )

  pA <- ggplot(d_full, aes(x = year_month)) +
    annotate("rect",
      xmin = INT_START, xmax = max(d_full$year_month) + days(31),
      ymin = -Inf, ymax = Inf, fill = COL_INT, alpha = 0.07) +
    geom_vline(xintercept = INT_START,
               linetype = "dashed", colour = "grey50", linewidth = 0.6) +
    geom_point(aes(y = pct_normal, size = n_temp,
                   colour = period2), alpha = 0.75) +
    scale_colour_manual(
      values = c("Historical" = COL_HIST, "Intervention" = COL_INT),
      name   = NULL
    ) +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    scale_y_continuous(limits = c(55, 100), labels = function(x) paste0(x, "%")) +
    labs(
      title = "Panel A: Monthly normal temperature -- historical vs intervention",
      subtitle = paste0(res$window_label,
        " | Vertical dashed line = November 2024 (intervention start)"),
      x = NULL, y = "% Normal admission temperature"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  # ---- Panel B: calendar-month dot plot (historical mean + intervention obs) -
  # d_int2 already has hist_mean/hist_lo95/hist_hi95 from run_comparison join
  d_cal <- d_int2

  pB <- ggplot() +
    # Historical mean and CI range
    geom_errorbar(data = d_cal,
      aes(x = factor(cal_month, labels = month_labels[sort(unique(cal_month))]),
          ymin = hist_lo95, ymax = hist_hi95),
      colour = COL_HIST, width = 0.3, linewidth = 0.8) +
    geom_point(data = d_cal,
      aes(x = factor(cal_month, labels = month_labels[sort(unique(cal_month))]),
          y = hist_mean),
      colour = COL_HIST, size = 3, shape = 18) +
    # Intervention observations
    geom_point(data = d_cal,
      aes(x = factor(cal_month, labels = month_labels[sort(unique(cal_month))]),
          y = pct_normal,
          colour = month_label),
      size = 4, alpha = 0.85) +
    # Connecting lines
    geom_segment(data = d_cal,
      aes(x    = factor(cal_month, labels = month_labels[sort(unique(cal_month))]),
          xend = factor(cal_month, labels = month_labels[sort(unique(cal_month))]),
          y    = hist_mean, yend = pct_normal,
          colour = month_label),
      linewidth = 0.7, alpha = 0.6) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "Panel B: Observed vs historical same-month mean",
      subtitle = "Blue diamond = historical mean (+/- 95% CI); coloured dots = intervention months",
      x = "Calendar month", y = "% Normal admission temperature",
      colour = "Intervention month"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  # ---- Panel C: difference (intervention minus historical mean) by month ---
  pC <- ggplot(d_cal,
    aes(x = year_month, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = COL_NULL, linewidth = 0.8) +
    geom_errorbar(aes(ymin = diff_lo95, ymax = diff_hi95),
                  colour = COL_DIFF, width = 12, linewidth = 0.7) +
    geom_point(aes(size = n_temp), colour = COL_DIFF, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x,
                colour = "grey40", fill = "grey85", alpha = 0.4,
                linewidth = 0.8, se = TRUE) +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    annotate("label",
      x = INT_START + days(10), y = min(d_cal$diff, na.rm = TRUE) * 0.85,
      label = test_annot, hjust = 0, size = 3,
      fill = "white", colour = "grey20") +
    labs(
      title = "Panel C: Difference (intervention month minus same-month historical mean)",
      subtitle = "Points above zero = intervention better than historical; error bars = propagated historical uncertainty",
      x = NULL, y = "Difference in % normal temperature (pp)"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  list(pA = pA, pB = pB, pC = pC)
}

log("")
log("Generating plots ...")

plots_prim <- make_plots(res_prim)
plots_sens <- make_plots(res_sens)

# Save primary window 3-panel plot
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined_prim <- (plots_prim$pA / plots_prim$pB / plots_prim$pC) +
    plot_annotation(
      title = "Script 04d: Same-Season Year-on-Year Comparison -- SMCH Golden Hour",
      subtitle = "Primary window (Jan 2023 historical baseline)"
    )
  ggsave(
    file.path(OUTPUTS_DIR, "04d_seasonal_comparison_plot.png"),
    combined_prim, width = 12, height = 15, dpi = 150
  )
  log("3-panel plot saved (patchwork): 04d_seasonal_comparison_plot.png")
} else {
  ggsave(
    file.path(OUTPUTS_DIR, "04d_seasonal_comparison_plot.png"),
    plots_prim$pC, width = 10, height = 6, dpi = 150
  )
  log("Difference panel saved (patchwork not available): 04d_seasonal_comparison_plot.png")
}


# -----------------------------------------------------------------------------
# 8. FINAL SUMMARY
# -----------------------------------------------------------------------------

log("")
log("=============================================================")
log("SUMMARY: SAME-SEASON YEAR-ON-YEAR COMPARISON")
log("=============================================================")
log("")
for (res in list(res_prim, res_sens)) {
  tt <- res$ttest
  log(paste0("Window: ", res$window_label))
  log(sprintf("  17 intervention months vs same calendar months in historical period"))
  log(sprintf("  Mean difference: %+.2f pp (95%% CI: [%+.2f, %+.2f])",
    mean(res$matched$diff, na.rm = TRUE),
    tt$conf.int[1], tt$conf.int[2]))
  log(sprintf("  t(%d) = %.3f, p = %.4f%s",
    tt$parameter, tt$statistic, tt$p.value,
    if (tt$p.value < 0.05) " *" else " (NS)"))
  log("")
}
log("Interpretation: A positive mean difference indicates intervention-period")
log("temperatures are higher than historical same-month values, suggesting")
log("improved thermal care. A non-significant difference indicates no detectable")
log("change after removing seasonal effects by month-matching.")
log("")
log("Script 04d complete.")
close(log_con)
