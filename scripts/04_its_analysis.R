# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 04: Interrupted Time Series (ITS) Analysis
# =============================================================================
# Purpose:
#   Evaluate whether the Golden Hour intervention (November 2024) was
#   associated with a change in the monthly proportion of babies with normal
#   admission temperature (36.5-37.5 C, WHO), using segmented regression ITS.
#
#   This addresses the limitation of the simple period comparison in Script 02
#   (Q9), which could not separate the intervention effect from the pre-existing
#   secular decline in temperature rates visible in the monthly trend plot.
#
# ITS model (segmented linear regression):
#   pct_normal ~ time + intervention + time_after
#
#   Where:
#     time         = sequential month number (1 = Jan 2022, 46 = Oct 2025)
#     intervention = 0 for all historical months, 1 for intervention months
#     time_after   = 0 for historical months, 1/2/3... for each intervention month
#
#   Coefficients:
#     (Intercept)  = estimated % normal temp at baseline (Jan 2022, t=1)
#     time         = pre-intervention monthly trend (%-pts per month)
#     intervention = immediate level change at Nov 2024 (%-pts step)
#     time_after   = change in monthly slope after intervention onset
#                    (post-intervention slope = time + time_after)
#
# Autocorrelation:
#   Monthly time series data are often autocorrelated. A Durbin-Watson test is
#   run on the OLS residuals. If significant autocorrelation is detected, the
#   primary inference uses Newey-West heteroscedasticity-and-autocorrelation-
#   consistent (HAC) standard errors (lag = 2). OLS and NW results are both
#   saved for transparency.
#
# Analyses run:
#   Primary   : full window, Jan 2022 - Oct 2025 (46 months)
#   Sensitivity: restricted window, Jan 2023 - Oct 2025 (34 months)
#                (excludes the anomalous early 2022 months with ~48-49%
#                 normal temperature, pending clinical input on whether
#                 these reflect a data quality issue or real ward conditions)
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)
library(lmtest)    # coeftest() for Newey-West inference
library(sandwich)  # NeweyWest() HAC standard errors
# Note: dwtest() does not support weighted regression; autocorrelation is
# instead assessed via Box.test() (Ljung-Box) on the weighted model residuals.

# -- Locate this script's directory -------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
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

# -- Redirect all output to a log file ----------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "04_its_analysis_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2026-03-31")   # extended from Oct 2025 to Mar 2026
HIST_START <- as.Date("2022-01-01")   # full primary window
HIST_START_SENS <- as.Date("2023-01-01")  # restricted sensitivity window
HIST_END   <- as.Date("2024-10-31")

# -- Minimum temperature recordings per month for inclusion -------------------
MIN_N_TEMP <- 5

# -- Plot colours (consistent with Script 03) ---------------------------------
COL_HIST  <- "#4A90C4"
COL_INT   <- "#E07B32"
FILL_HIST <- "#DDEEF8"
FILL_INT  <- "#FEF0DC"


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)

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

    temperature_num = suppressWarnings(as.numeric(temperature)),

    temp_normal = case_when(
      is.na(temperature_num)                             ~ NA,
      temperature_num >= 36.5 & temperature_num <= 37.5 ~ TRUE,
      TRUE                                               ~ FALSE
    ),

    period = case_when(
      adm_date >= INT_START  & adm_date <= INT_END   ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END  ~ "Historical",
      TRUE                                            ~ "Other"
    ),

    is_readmission = readmission == "Y"
  )

cat("Derived variables created.\n")


# -----------------------------------------------------------------------------
# 3. DEFINE ANALYSIS POPULATION
# -----------------------------------------------------------------------------

smch_inborn <- df %>%
  filter(
    facility      == "SMCH",
    inborn        == TRUE,
    is_readmission == FALSE | is.na(is_readmission)
  )

q10_pop <- smch_inborn %>%
  filter(period %in% c("Intervention", "Historical"))

cat(sprintf("\nQ10 population (SMCH inborn, Jan 2022 - Mar 2026): %d records\n",
            nrow(q10_pop)))


# -----------------------------------------------------------------------------
# 4. BUILD MONTHLY DATA
# -----------------------------------------------------------------------------

monthly_raw <- q10_pop %>%
  group_by(adm_month, period) %>%
  summarise(
    n_admissions    = n(),
    n_temp_recorded = sum(!is.na(temperature_num)),
    n_normal        = sum(temp_normal == TRUE, na.rm = TRUE),
    pct_normal      = if_else(
      n_temp_recorded > 0,
      round(100 * n_normal / n_temp_recorded, 2),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  filter(n_temp_recorded >= MIN_N_TEMP) %>%
  arrange(adm_month) %>%
  mutate(
    adm_month_date = as.Date(adm_month)
  )

cat(sprintf("\nMonthly data points available: %d\n", nrow(monthly_raw)))
cat(sprintf("  Historical months : %d\n", sum(monthly_raw$period == "Historical")))
cat(sprintf("  Intervention months: %d\n", sum(monthly_raw$period == "Intervention")))


# -----------------------------------------------------------------------------
# 5. ITS MODEL HELPER FUNCTION
# -----------------------------------------------------------------------------

# fit_its()
# Fits the segmented regression ITS model to a given monthly dataset.
# Assumes the dataset is sorted by adm_month (ascending) and covers a
# contiguous run of months with no gaps.
#
# Returns a list with:
#   $model      -- OLS lm() object
#   $its_data   -- input data augmented with ITS variables and fitted values
#   $dw         -- Durbin-Watson test result
#   $coef_ols   -- tidy OLS coefficient table
#   $coef_nw    -- tidy Newey-West corrected coefficient table
#   $label      -- descriptive label for the analysis
#
fit_its <- function(monthly_data, label) {

  cat(sprintf("\n%s\n%s\n", strrep("=", 64), label))
  cat(sprintf("%s\n\n", strrep("=", 64)))

  # -- Build ITS variables ---------------------------------------------------
  # time        : 1, 2, 3, ... N (sequential month index)
  # intervention: 0 pre-intervention, 1 post-intervention onset
  # time_after  : 0 pre-intervention, 1, 2, 3... months since intervention
  d <- monthly_data %>%
    arrange(adm_month) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )

  cat("ITS variable summary:\n")
  cat(sprintf("  Total months       : %d\n", nrow(d)))
  cat(sprintf("  Historical months  : %d (time 1 to %d)\n",
              sum(d$intervention == 0), max(d$time[d$intervention == 0])))
  cat(sprintf("  Intervention months: %d (time %d to %d)\n",
              sum(d$intervention == 1),
              min(d$time[d$intervention == 1]),
              max(d$time[d$intervention == 1])))
  cat(sprintf("  Intervention start : %s\n",
              format(min(as.Date(d$adm_month[d$intervention == 1])), "%B %Y")))

  # -- Fit weighted OLS segmented regression ---------------------------------
  # Weights = n_temp_recorded so that high-volume months contribute more.
  model <- lm(
    pct_normal ~ time + intervention + time_after,
    data    = d,
    weights = n_temp_recorded
  )

  cat("\n--- OLS model summary ---\n")
  print(summary(model))

  # -- Ljung-Box test for autocorrelation (lag 1) ----------------------------
  # Box.test() works on the raw residuals of any model including weighted lm.
  # Lag = 1 is appropriate for monthly series with potential AR(1) structure.
  lb <- Box.test(residuals(model), lag = 1, type = "Ljung-Box")
  cat(sprintf("\n--- Ljung-Box test for autocorrelation (lag 1) ---\n"))
  cat(sprintf("  Q statistic = %.4f  df = %d\n", lb$statistic, lb$parameter))
  cat(sprintf("  p-value     = %.4f\n", lb$p.value))
  if (lb$p.value < 0.05) {
    cat("  >> Significant autocorrelation detected. Newey-West SEs applied.\n")
  } else {
    cat("  >> No significant autocorrelation. OLS SEs are reliable.\n")
  }

  # -- Newey-West HAC standard errors (lag = 2 for monthly series) ----------
  # Reported regardless of DW result for transparency; use NW if DW p < 0.05.
  nw_test <- coeftest(model, vcov = NeweyWest(model, lag = 2, prewhite = FALSE))

  cat("\n--- Newey-West corrected coefficients (lag = 2) ---\n")
  print(nw_test)

  # -- Tidy coefficient tables -----------------------------------------------
  # Extract by column position rather than name to avoid make.names() mangling:
  # as.data.frame() converts "Std. Error" -> "Std..Error", "t value" -> "t.value",
  # etc., so rename() with backtick-quoted original names fails on some R versions.

  ols_mat  <- summary(model)$coefficients          # always 4 cols: Est, SE, t, p
  tidy_ols <- data.frame(
    analysis  = label,
    term      = rownames(ols_mat),
    estimate  = ols_mat[, 1],
    std_error = ols_mat[, 2],
    t_value   = ols_mat[, 3],
    p_value   = ols_mat[, 4],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      ci_lower = estimate - 1.96 * std_error,
      ci_upper = estimate + 1.96 * std_error,
      across(where(is.numeric), ~ round(.x, 4)),
      se_type  = "OLS"
    ) %>%
    select(analysis, term, estimate, std_error, ci_lower, ci_upper, t_value, p_value, se_type)

  nw_mat   <- as.matrix(nw_test)                   # always 4 cols: Est, SE, z, p
  tidy_nw  <- data.frame(
    analysis  = label,
    term      = rownames(nw_mat),
    estimate  = nw_mat[, 1],
    std_error = nw_mat[, 2],
    t_value   = nw_mat[, 3],
    p_value   = nw_mat[, 4],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      ci_lower = estimate - 1.96 * std_error,
      ci_upper = estimate + 1.96 * std_error,
      across(where(is.numeric), ~ round(.x, 4)),
      se_type  = "Newey-West (lag=2)"
    ) %>%
    select(analysis, term, estimate, std_error, ci_lower, ci_upper, t_value, p_value, se_type)

  # -- Plain-English interpretation ------------------------------------------
  cf  <- coef(model)
  b0  <- cf["(Intercept)"]
  b1  <- cf["time"]
  b2  <- cf["intervention"]
  b3  <- cf["time_after"]

  cat(sprintf("\n--- Interpretation ---\n"))
  cat(sprintf(
    "  Pre-intervention trend (b1)    : %.2f %%-pts/month (%.1f %%-pts/year)\n",
    b1, b1 * 12
  ))
  cat(sprintf(
    "  Level change at Nov 2024 (b2)  : %.2f %%-pts (immediate step)\n", b2
  ))
  cat(sprintf(
    "  Slope change post-int  (b3)    : %.2f %%-pts/month\n", b3
  ))
  cat(sprintf(
    "  Post-intervention slope (b1+b3): %.2f %%-pts/month (%.1f %%-pts/year)\n",
    b1 + b3, (b1 + b3) * 12
  ))

  # -- Fitted values and counterfactual --------------------------------------
  d <- d %>%
    mutate(
      fitted_its     = fitted(model),
      # Counterfactual: pre-intervention trend projected forward (no step, no slope change)
      counterfactual = b0 + b1 * time
    )

  list(
    model    = model,
    its_data = d,
    lb       = lb,
    coef_ols = tidy_ols,
    coef_nw  = tidy_nw,
    label    = label
  )
}


# -----------------------------------------------------------------------------
# 6. PRIMARY ITS ANALYSIS: FULL WINDOW (JAN 2022 - OCT 2025)
# -----------------------------------------------------------------------------

primary_data <- monthly_raw

its_primary <- fit_its(
  monthly_data = primary_data,
  label        = "ITS primary: Jan 2022 - Mar 2026 (51 months)"
)


# -----------------------------------------------------------------------------
# 7. SENSITIVITY ANALYSIS: RESTRICTED WINDOW (JAN 2023 - OCT 2025)
# -----------------------------------------------------------------------------
# Excludes the anomalous Jan-Feb 2022 months (~48-49% normal temperature)
# which may reflect data quality issues in the early months of Neotree
# recording at SMCH. Pending clinical confirmation from the SMCH team.

cat("\n")
cat("----------------------------------------------------------------\n")
cat("  NOTE: Sensitivity analysis excludes all of 2022 (starts Jan\n")
cat("  2023) pending clinical confirmation of data quality in that\n")
cat("  period (anomalously low normal-temp rates throughout 2022).\n")
cat("----------------------------------------------------------------\n")

sensitivity_data <- monthly_raw %>%
  filter(as.Date(adm_month) >= HIST_START_SENS)

its_sensitivity <- fit_its(
  monthly_data = sensitivity_data,
  label        = "ITS sensitivity: Jan 2023 - Mar 2026 (39 months)"
)


# -----------------------------------------------------------------------------
# 8. SAVE COEFFICIENT TABLES
# -----------------------------------------------------------------------------

# Combined table of all coefficients from both analyses
all_coefs <- bind_rows(
  its_primary$coef_ols,
  its_primary$coef_nw,
  its_sensitivity$coef_ols,
  its_sensitivity$coef_nw
)

write_csv(all_coefs,
          file.path(OUTPUT_DIR, "04a_its_coefficients.csv"))
cat("\nSaved: 04a_its_coefficients.csv\n")

# Monthly data with ITS variables and fitted values (primary)
write_csv(
  its_primary$its_data %>%
    select(adm_month_date, period, n_admissions, n_temp_recorded,
           n_normal, pct_normal, time, intervention, time_after,
           fitted_its, counterfactual),
  file.path(OUTPUT_DIR, "04b_its_monthly_fitted_primary.csv")
)
cat("Saved: 04b_its_monthly_fitted_primary.csv\n")

# Monthly data for sensitivity analysis
write_csv(
  its_sensitivity$its_data %>%
    select(adm_month_date, period, n_admissions, n_temp_recorded,
           n_normal, pct_normal, time, intervention, time_after,
           fitted_its, counterfactual),
  file.path(OUTPUT_DIR, "04c_its_monthly_fitted_sensitivity.csv")
)
cat("Saved: 04c_its_monthly_fitted_sensitivity.csv\n")


# -----------------------------------------------------------------------------
# 9. ITS PLOT (PRIMARY ANALYSIS)
# -----------------------------------------------------------------------------

cat("\nBuilding ITS plot (primary analysis)...\n")

pd <- its_primary$its_data

# Split fitted line into two segments for plotting (pre / post)
pd_pre  <- pd %>% filter(intervention == 0)
pd_post <- pd %>% filter(intervention == 1)

# Counterfactual: only plotted over the intervention period
pd_cf <- pd %>% filter(intervention == 1)

# Period x-range for shading
hist_xmin <- min(as.Date(pd$adm_month))
hist_xmax <- as.Date(INT_START)
int_xmin  <- as.Date(INT_START)
int_xmax  <- max(as.Date(pd$adm_month)) + days(31)

p_its <- ggplot(pd, aes(x = as.Date(adm_month))) +

  # -- Period shading ----------------------------------------------------------
  annotate("rect",
           xmin = hist_xmin, xmax = hist_xmax,
           ymin = -Inf,      ymax = Inf,
           fill = FILL_HIST, alpha = 1) +
  annotate("rect",
           xmin = int_xmin, xmax = int_xmax,
           ymin = -Inf,     ymax = Inf,
           fill = FILL_INT, alpha = 1) +

  # -- Intervention onset line ------------------------------------------------
  geom_vline(xintercept = INT_START,
             linetype = "dashed", colour = "grey45", linewidth = 0.55) +

  # -- Counterfactual (pre-intervention trend extended) -----------------------
  geom_line(data = pd_cf,
            aes(y = counterfactual),
            linetype = "dotted", colour = COL_HIST,
            linewidth = 0.9) +

  # -- Fitted ITS lines (pre and post) ----------------------------------------
  geom_line(data = pd_pre,
            aes(y = fitted_its),
            colour = COL_HIST, linewidth = 1.2) +
  geom_line(data = pd_post,
            aes(y = fitted_its),
            colour = COL_INT, linewidth = 1.2) +

  # -- Observed data points ---------------------------------------------------
  geom_point(aes(y = pct_normal,
                 colour = period,
                 size   = n_temp_recorded),
             alpha = 0.75) +

  # -- Period label annotations -----------------------------------------------
  annotate("text",
           x = as.Date("2023-01-15"), y = 96,
           label    = "Historical period",
           colour   = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x = as.Date("2025-04-01"), y = 96,
           label    = "Intervention\nperiod",
           colour   = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +

  # -- Counterfactual label ---------------------------------------------------
  annotate("text",
           x     = as.Date("2025-07-01"),
           y     = its_primary$its_data$counterfactual[
             which.min(abs(as.Date(its_primary$its_data$adm_month) -
                             as.Date("2025-07-01")))
           ] + 3,
           label  = "Counterfactual\n(pre-trend projected)",
           colour = COL_HIST, hjust = 0.5, size = 2.8, fontface = "italic") +

  # -- Scales -----------------------------------------------------------------
  scale_colour_manual(
    name   = "Period",
    values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
  ) +
  scale_size_continuous(
    name  = "Temps\nrecorded",
    range = c(1.5, 5),
    breaks = c(50, 100, 200, 400)
  ) +
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b\n%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(NA, 100),
    expand = expansion(mult = c(0.05, 0.08))
  ) +

  # -- Labels -----------------------------------------------------------------
  labs(
    title    = "Interrupted time series: normal admission temperature, SMCH Jan 2022 - Mar 2026",
    subtitle = paste0(
      "Solid lines = ITS fitted trend (pre- and post-intervention). ",
      "Dotted line = counterfactual (pre-intervention trend extended)."
    ),
    x       = NULL,
    y       = "% with normal admission temperature (36.5-37.5 C)",
    caption = paste0(
      "Segmented linear regression: pct_normal ~ time + intervention + time_after.",
      "\nWeighted by n temperatures recorded per month.",
      "\nDashed vertical line = intervention start (November 2024).",
      "\nPoint size proportional to number of temperature recordings."
    )
  ) +

  # -- Theme ------------------------------------------------------------------
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11.5),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    plot.caption     = element_text(colour = "grey50", size = 7.5, hjust = 0),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 9),
    axis.title.y     = element_text(size = 9.5),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.5, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

PNG_PATH <- file.path(OUTPUT_DIR, "04_its_plot_primary.png")
ggsave(
  filename = PNG_PATH,
  plot     = p_its,
  width    = 12,
  height   = 5.5,
  dpi      = 300,
  units    = "in"
)
cat(sprintf("Saved: 04_its_plot_primary.png  (%dx%d px @ 300 dpi)\n",
            12 * 300, round(5.5 * 300)))


# -----------------------------------------------------------------------------
# 9b. ITS PLOT (SENSITIVITY ANALYSIS: JAN 2023 - OCT 2025)
# -----------------------------------------------------------------------------

cat("\nBuilding ITS plot (sensitivity analysis)...\n")

pd_s <- its_sensitivity$its_data

pd_s_pre  <- pd_s %>% filter(intervention == 0)
pd_s_post <- pd_s %>% filter(intervention == 1)
pd_s_cf   <- pd_s %>% filter(intervention == 1)

hist_xmin_s <- min(as.Date(pd_s$adm_month))
hist_xmax_s <- as.Date(INT_START)
int_xmin_s  <- as.Date(INT_START)
int_xmax_s  <- max(as.Date(pd_s$adm_month)) + days(31)

# Coefficients for counterfactual label y-position
cf_s  <- coef(its_sensitivity$model)
b0_s  <- cf_s["(Intercept)"]
b1_s  <- cf_s["time"]

p_its_s <- ggplot(pd_s, aes(x = as.Date(adm_month))) +

  # -- Period shading ----------------------------------------------------------
  annotate("rect",
           xmin = hist_xmin_s, xmax = hist_xmax_s,
           ymin = -Inf,        ymax = Inf,
           fill = FILL_HIST, alpha = 1) +
  annotate("rect",
           xmin = int_xmin_s, xmax = int_xmax_s,
           ymin = -Inf,       ymax = Inf,
           fill = FILL_INT, alpha = 1) +

  # -- Intervention onset line ------------------------------------------------
  geom_vline(xintercept = INT_START,
             linetype = "dashed", colour = "grey45", linewidth = 0.55) +

  # -- Counterfactual (pre-intervention trend extended) -----------------------
  geom_line(data = pd_s_cf,
            aes(y = counterfactual),
            linetype = "dotted", colour = COL_HIST,
            linewidth = 0.9) +

  # -- Fitted ITS lines (pre and post) ----------------------------------------
  geom_line(data = pd_s_pre,
            aes(y = fitted_its),
            colour = COL_HIST, linewidth = 1.2) +
  geom_line(data = pd_s_post,
            aes(y = fitted_its),
            colour = COL_INT, linewidth = 1.2) +

  # -- Observed data points ---------------------------------------------------
  geom_point(aes(y = pct_normal,
                 colour = period,
                 size   = n_temp_recorded),
             alpha = 0.75) +

  # -- Period label annotations -----------------------------------------------
  annotate("text",
           x = as.Date("2023-12-15"), y = 96,
           label    = "Historical period\n(Jan 2023 - Oct 2024)",
           colour   = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x = as.Date("2025-04-01"), y = 96,
           label    = "Intervention\nperiod",
           colour   = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +

  # -- Counterfactual label ---------------------------------------------------
  annotate("text",
           x     = as.Date("2025-07-01"),
           y     = pd_s_cf$counterfactual[
             which.min(abs(as.Date(pd_s_cf$adm_month) -
                             as.Date("2025-07-01")))
           ] + 3,
           label  = "Counterfactual\n(pre-trend projected)",
           colour = COL_HIST, hjust = 0.5, size = 2.8, fontface = "italic") +

  # -- Scales -----------------------------------------------------------------
  scale_colour_manual(
    name   = "Period",
    values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
  ) +
  scale_size_continuous(
    name  = "Temps\nrecorded",
    range = c(1.5, 5),
    breaks = c(50, 100, 200, 400)
  ) +
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b\n%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(NA, 100),
    expand = expansion(mult = c(0.05, 0.08))
  ) +

  # -- Labels -----------------------------------------------------------------
  labs(
    title    = "ITS sensitivity analysis: normal admission temperature, SMCH Jan 2023 - Mar 2026",
    subtitle = paste0(
      "Restricted window: excludes 2022 (anomalous recording period; pending clinical review). ",
      "Dotted line = counterfactual."
    ),
    x       = NULL,
    y       = "% with normal admission temperature (36.5-37.5 C)",
    caption = paste0(
      "Segmented linear regression: pct_normal ~ time + intervention + time_after.",
      "\nWeighted by n temperatures recorded per month.",
      "\nDashed vertical line = intervention start (November 2024).",
      "\nPoint size proportional to number of temperature recordings."
    )
  ) +

  # -- Theme ------------------------------------------------------------------
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11.5),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    plot.caption     = element_text(colour = "grey50", size = 7.5, hjust = 0),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 9),
    axis.title.y     = element_text(size = 9.5),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.5, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

PNG_PATH_S <- file.path(OUTPUT_DIR, "04_its_plot_sensitivity.png")
ggsave(
  filename = PNG_PATH_S,
  plot     = p_its_s,
  width    = 12,
  height   = 5.5,
  dpi      = 300,
  units    = "in"
)
cat(sprintf("Saved: 04_its_plot_sensitivity.png  (%dx%d px @ 300 dpi)\n",
            12 * 300, round(5.5 * 300)))


# -----------------------------------------------------------------------------
# 10. SUMMARY FOR THE MANUSCRIPT
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  ITS RESULTS SUMMARY FOR MANUSCRIPT\n")
cat("================================================================\n")

for (its in list(its_primary, its_sensitivity)) {

  cat(sprintf("\n%s\n", its$label))

  # Use NW SEs if autocorrelation was significant
  use_nw <- its$lb$p.value < 0.05
  coefs  <- if (use_nw) its$coef_nw else its$coef_ols
  se_note <- if (use_nw) "Newey-West HAC SEs (autocorrelation detected)" else "OLS SEs (no autocorrelation)"
  cat(sprintf("  Standard errors used: %s\n\n", se_note))

  term_labels <- c(
    "time"         = "Pre-intervention monthly trend (%-pts/month)",
    "intervention" = "Level change at Nov 2024 (%-pts, immediate step)",
    "time_after"   = "Slope change post-intervention (%-pts/month)"
  )

  for (trm in names(term_labels)) {
    row <- coefs %>% filter(term == trm)
    if (nrow(row) == 1) {
      cat(sprintf(
        "  %s\n    Estimate: %.2f (95%% CI: %.2f to %.2f)  p = %.4f\n\n",
        term_labels[trm],
        row$estimate, row$ci_lower, row$ci_upper, row$p_value
      ))
    }
  }
}

cat("\n=== Script 04 complete ===\n")

# -- Close log -----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
