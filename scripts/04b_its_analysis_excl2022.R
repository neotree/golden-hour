# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 04b: ITS Analysis -- Pre-intervention window excludes 2022
# =============================================================================
# Purpose:
#   Repeat the segmented regression ITS from Script 04 using a pre-intervention
#   window that starts January 2023, excluding all 2022 data.
#
# Rationale for excluding 2022:
#   (a) Clinical/data-quality: temperature recording rates in 2022 were
#       anomalously low (~48-55% normal vs ~75-85% in 2023-2024). The clinical
#       explanation has not yet been confirmed by the SMCH team.
#   (b) Statistical: the power analysis (Script 05) showed that including 2022
#       inflates the pre-intervention residual sigma from ~122 to ~196 %-pts,
#       roughly halving power for the level-change estimand (b2).
#   (c) ITS design balance: the ratio of pre- to post-intervention months is
#       22:12 when 2022 is excluded, compared to 34:12 when included. A more
#       balanced pre/post split improves precision for b2 and b3.
#
# This script therefore treats the Jan 2023 window as the PRIMARY analysis.
# The full-window results (Jan 2022) remain available in Script 04 for
# reference and are the appropriate sensitivity check once the SMCH team
# confirms or refutes the 2022 data quality concern.
#
# Analysis window:
#   Pre-intervention : January 2023 - October 2024  (22 months)
#   Intervention     : November 2024 - March 2026   (17 months)
#   Total            : 39 months
#
# ITS model (identical to Script 04):
#   pct_normal ~ time + intervention + time_after
#   Weighted by monthly n temperature recordings.
#   Ljung-Box test for autocorrelation; Newey-West HAC SEs (lag = 2) applied
#   if autocorrelation is detected.
#
# Outputs (04b_ prefix; all to 03-OUTPUTS/):
#   04b_its_coefficients.csv          -- OLS and NW coefficient tables
#   04b_its_monthly_fitted.csv        -- monthly data with ITS vars + fitted values
#   04b_its_plot.png                  -- ITS time series plot
#   04b_its_analysis_log.txt          -- full log
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

args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
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

LOG_FILE <- file.path(OUTPUT_DIR, "04b_its_analysis_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour ITS Analysis (excluding 2022) -- SMCH Zimbabwe\n")
cat(sprintf("Run started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2026-03-31")   # extended from Oct 2025 to Mar 2026
HIST_START <- as.Date("2023-01-01")   # excludes all 2022 data
HIST_END   <- as.Date("2024-10-31")

MIN_N_TEMP <- 5

COL_HIST  <- "#4A90C4"
COL_INT   <- "#E07B32"
FILL_HIST <- "#DDEEF8"
FILL_INT  <- "#FEF0DC"
COL_CF    <- "#7B9E87"


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
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


# -----------------------------------------------------------------------------
# 3. ANALYSIS POPULATION
# -----------------------------------------------------------------------------

smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission)
  )

q10_pop <- smch_inborn %>%
  filter(period %in% c("Intervention", "Historical"))

cat(sprintf("\nAnalysis population (SMCH inborn, Jan 2023 - Mar 2026): %d records\n",
            nrow(q10_pop)))


# -----------------------------------------------------------------------------
# 4. BUILD MONTHLY DATA
# -----------------------------------------------------------------------------

monthly <- q10_pop %>%
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
  mutate(adm_month_date = as.Date(adm_month))

cat(sprintf("\nMonthly data points: %d total (%d pre-intervention, %d intervention)\n",
            nrow(monthly),
            sum(monthly$period == "Historical"),
            sum(monthly$period == "Intervention")))


# -----------------------------------------------------------------------------
# 5. FIT ITS MODEL
# -----------------------------------------------------------------------------

cat(sprintf("\n%s\n", strrep("=", 64)))
cat("ITS PRIMARY ANALYSIS: Jan 2023 - Mar 2026 (22 pre + 17 post)\n")
cat(sprintf("%s\n\n", strrep("=", 64)))

# Build ITS variables
d <- monthly %>%
  arrange(adm_month) %>%
  mutate(
    time         = row_number(),
    intervention = as.integer(period == "Intervention"),
    time_after   = cumsum(intervention)
  )

cat("ITS variable summary:\n")
cat(sprintf("  Total months        : %d\n", nrow(d)))
cat(sprintf("  Pre-intervention    : %d (time 1 to %d)\n",
            sum(d$intervention == 0), max(d$time[d$intervention == 0])))
cat(sprintf("  Intervention        : %d (time %d to %d)\n",
            sum(d$intervention == 1),
            min(d$time[d$intervention == 1]),
            max(d$time[d$intervention == 1])))
cat(sprintf("  Pre-intervention start: %s\n",
            format(min(as.Date(d$adm_month[d$intervention == 0])), "%B %Y")))
cat(sprintf("  Intervention start  : %s\n",
            format(min(as.Date(d$adm_month[d$intervention == 1])), "%B %Y")))

# Fit weighted OLS segmented regression
model <- lm(
  pct_normal ~ time + intervention + time_after,
  data    = d,
  weights = n_temp_recorded
)

cat("\n--- OLS model summary ---\n")
print(summary(model))

# Ljung-Box autocorrelation test (lag 1)
lb <- Box.test(residuals(model), lag = 1, type = "Ljung-Box")
cat(sprintf("\n--- Ljung-Box test for autocorrelation (lag 1) ---\n"))
cat(sprintf("  Q statistic = %.4f  df = %d\n", lb$statistic, lb$parameter))
cat(sprintf("  p-value     = %.4f\n", lb$p.value))
if (lb$p.value < 0.05) {
  cat("  >> Significant autocorrelation detected. Newey-West SEs applied.\n")
} else {
  cat("  >> No significant autocorrelation. OLS SEs are reliable.\n")
}

# Newey-West HAC standard errors (lag = 2)
nw_test <- coeftest(model, vcov = NeweyWest(model, lag = 2, prewhite = FALSE))
cat("\n--- Newey-West corrected coefficients (lag = 2) ---\n")
print(nw_test)

# Coefficient extraction by column position (robust to column-name mangling)
ols_mat <- summary(model)$coefficients
tidy_ols <- data.frame(
  analysis  = "ITS primary (excl 2022): Jan 2023 - Oct 2025",
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
  select(analysis, term, estimate, std_error, ci_lower, ci_upper,
         t_value, p_value, se_type)

nw_mat  <- as.matrix(nw_test)
tidy_nw <- data.frame(
  analysis  = "ITS primary (excl 2022): Jan 2023 - Oct 2025",
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
  select(analysis, term, estimate, std_error, ci_lower, ci_upper,
         t_value, p_value, se_type)

# Plain-English interpretation
cf <- coef(model)
b1 <- cf["time"]
b2 <- cf["intervention"]
b3 <- cf["time_after"]

cat("\n--- Interpretation ---\n")
cat(sprintf("  Pre-intervention trend (b1)    : %.2f %%-pts/month (%.1f %%-pts/year)\n",
            b1, b1 * 12))
cat(sprintf("  Level change at Nov 2024 (b2)  : %.2f %%-pts (immediate step)\n", b2))
cat(sprintf("  Slope change post-int  (b3)    : %.2f %%-pts/month\n", b3))
cat(sprintf("  Post-intervention slope (b1+b3): %.2f %%-pts/month (%.1f %%-pts/year)\n",
            b1 + b3, (b1 + b3) * 12))

# Fitted values and counterfactual
b0 <- cf["(Intercept)"]
d <- d %>%
  mutate(
    fitted_its     = fitted(model),
    counterfactual = b0 + b1 * time
  )

# CI bands (OLS, for plotting)
ci_fit <- predict(model, interval = "confidence", level = 0.95)
vcv    <- vcov(model)
cf_se  <- sqrt(
  vcv["(Intercept)", "(Intercept)"] +
  d$time^2 * vcv["time", "time"] +
  2 * d$time * vcv["(Intercept)", "time"]
)
d <- d %>%
  mutate(
    fit_lower = ci_fit[, "lwr"],
    fit_upper = ci_fit[, "upr"],
    cf_lower  = counterfactual - 1.96 * cf_se,
    cf_upper  = counterfactual + 1.96 * cf_se
  )


# -----------------------------------------------------------------------------
# 6. SAVE OUTPUTS
# -----------------------------------------------------------------------------

all_coefs <- bind_rows(tidy_ols, tidy_nw)
write_csv(all_coefs, file.path(OUTPUT_DIR, "04b_its_coefficients.csv"))
cat("\nSaved: 04b_its_coefficients.csv\n")

write_csv(
  d %>% select(adm_month_date, period, n_admissions, n_temp_recorded,
               n_normal, pct_normal, time, intervention, time_after,
               fitted_its, counterfactual, fit_lower, fit_upper,
               cf_lower, cf_upper),
  file.path(OUTPUT_DIR, "04b_its_monthly_fitted.csv")
)
cat("Saved: 04b_its_monthly_fitted.csv\n")


# -----------------------------------------------------------------------------
# 7. ITS PLOT
# -----------------------------------------------------------------------------

cat("\nBuilding ITS plot...\n")

d_pre  <- d %>% filter(intervention == 0)
d_post <- d %>% filter(intervention == 1)
d_cf   <- d %>% filter(intervention == 1)

hist_xmin <- min(d$adm_month_date)
hist_xmax <- INT_START
int_xmin  <- INT_START
int_xmax  <- max(d$adm_month_date) + days(31)

p_its <- ggplot(d, aes(x = adm_month_date)) +

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

  # -- Counterfactual CI ribbon and line --------------------------------------
  geom_ribbon(data = d_cf,
              aes(ymin = cf_lower, ymax = cf_upper),
              fill = COL_CF, alpha = 0.18) +
  geom_line(data = d_cf,
            aes(y = counterfactual),
            colour = COL_CF, linetype = "dotted", linewidth = 1.0) +

  # -- Fitted ITS CI ribbons --------------------------------------------------
  geom_ribbon(data = d_pre,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_HIST, alpha = 0.20) +
  geom_ribbon(data = d_post,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_INT, alpha = 0.20) +

  # -- Fitted ITS lines -------------------------------------------------------
  geom_line(data = d_pre,
            aes(y = fitted_its),
            colour = COL_HIST, linewidth = 1.2) +
  geom_line(data = d_post,
            aes(y = fitted_its),
            colour = COL_INT, linewidth = 1.2) +

  # -- Observed data points ---------------------------------------------------
  geom_point(aes(y = pct_normal, colour = period, size = n_temp_recorded),
             alpha = 0.80) +

  # -- Labels -----------------------------------------------------------------
  annotate("text",
           x = as.Date("2023-12-01"), y = 97,
           label  = "Historical period\n(Jan 2023 - Oct 2024)",
           colour = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x = as.Date("2025-07-01"), y = 97,
           label  = "Intervention\nperiod",
           colour = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x      = as.Date("2025-09-15"),
           y      = d_cf$counterfactual[
             which.min(abs(d_cf$adm_month_date - as.Date("2025-09-15")))
           ] - 4.5,
           label  = "Counterfactual\n(pre-trend projected)",
           colour = COL_CF, hjust = 0.5, size = 2.7, fontface = "italic") +

  # -- Scales -----------------------------------------------------------------
  scale_colour_manual(
    name   = "Period",
    values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
  ) +
  scale_size_continuous(
    name   = "Temps\nrecorded",
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
    limits = c(NA, 100),
    expand = expansion(mult = c(0.05, 0.10))
  ) +

  labs(
    title    = "ITS (excl. 2022): normal admission temperature, SMCH Jan 2023 - Mar 2026",
    subtitle = paste0(
      "22 pre-intervention months (Jan 2023 - Oct 2024) + 17 intervention months. ",
      "Shaded ribbons = 95% OLS CI on fitted trends."
    ),
    x       = NULL,
    y       = "% with normal admission temperature (36.5-37.5 C)",
    caption = paste0(
      "Segmented linear regression: pct_normal ~ time + intervention + time_after, ",
      "weighted by n temperatures recorded per month.\n",
      "Dashed vertical line = intervention start (November 2024). ",
      "Dotted line = counterfactual (pre-intervention trend projected forward).\n",
      "Point size proportional to n temperature recordings. ",
      if (lb$p.value < 0.05)
        "Inference: Newey-West HAC SEs (autocorrelation detected, Ljung-Box p < 0.05)."
      else
        "Inference: OLS SEs (no significant autocorrelation, Ljung-Box p >= 0.05)."
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

ggsave(
  file.path(OUTPUT_DIR, "04b_its_plot.png"),
  plot = p_its, width = 13, height = 6, dpi = 300, units = "in"
)
cat(sprintf("Saved: 04b_its_plot.png\n"))


# -----------------------------------------------------------------------------
# 8. MANUSCRIPT SUMMARY
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  ITS RESULTS SUMMARY (excl. 2022)\n")
cat("================================================================\n\n")

use_nw  <- lb$p.value < 0.05
coefs   <- if (use_nw) tidy_nw else tidy_ols
se_note <- ifelse(use_nw,
                  "Newey-West HAC SEs (autocorrelation detected)",
                  "OLS SEs (no significant autocorrelation)")
cat(sprintf("  Standard errors: %s\n\n", se_note))

term_labels <- c(
  "time"         = "Pre-intervention monthly trend (%-pts/month)",
  "intervention" = "Level change at Nov 2024 (%-pts, immediate step)",
  "time_after"   = "Slope change post-intervention (%-pts/month)"
)

for (trm in names(term_labels)) {
  row <- coefs %>% filter(term == trm)
  if (nrow(row) == 1) {
    pval_str <- if (row$p_value < 0.001) "<0.001" else
                sprintf("%.4f", row$p_value)
    cat(sprintf("  %s\n    Estimate: %.2f (95%% CI: %.2f to %.2f)  p = %s\n\n",
                term_labels[trm],
                row$estimate, row$ci_lower, row$ci_upper, pval_str))
  }
}

cat("  Note: compare with Script 04 (full window, Jan 2022) to assess\n")
cat("  sensitivity of results to inclusion of the 2022 anomalous period.\n")

cat(sprintf("\nRun completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=== Script 04b complete ===\n")

sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
