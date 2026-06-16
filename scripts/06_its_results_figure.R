# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 06: ITS Results Figure (publication-quality)
# =============================================================================
# Purpose:
#   Produce a polished, publication-ready graphical summary of the ITS analysis
#   results from Script 04. Reads the Script 04 output CSVs rather than
#   re-running the full data pipeline. Re-fits the weighted ITS models from
#   the saved monthly data in order to compute confidence interval bands on the
#   fitted lines and counterfactual trajectory (these cannot be stored in CSVs
#   because they require the variance-covariance matrix of the model).
#
# Dependency:
#   Script 04 (04_its_analysis.R) must have been run first; this script reads:
#     04a_its_coefficients.csv           -- NW and OLS coefficient tables
#     04b_its_monthly_fitted_primary.csv -- monthly series, primary window
#     04c_its_monthly_fitted_sensitivity.csv -- monthly series, sensitivity window
#
# Outputs produced:
#   06a_its_timeseries.png    -- Panel A: enhanced ITS time series (primary)
#   06b_its_forest.png        -- Panel B: coefficient forest plot (both windows)
#   06_its_combined.png       -- Combined A + B figure (requires cowplot)
#   06_its_figure_log.txt     -- log
#
# Figure design:
#   Panel A  Enhanced version of the Script 04 primary ITS plot, adding:
#     - 95% CI ribbons around the pre- and post-intervention fitted trends
#     - 95% CI ribbon around the counterfactual trajectory
#     - Annotation of key months (early high-performing months and drops)
#     CI bands use OLS variance (standard for ITS visualisations); NW SEs
#     are used for inference in Panel B.
#
#   Panel B  Forest plot comparing primary and sensitivity analyses for the
#     two key ITS estimands (b2 level change, b3 slope change), using
#     Newey-West 95% CIs from 04a_its_coefficients.csv.
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
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag), mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUTPUT_DIR, "06_its_figure_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour -- ITS Results Figure\n")
cat(sprintf("Run started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Constants (consistent with Scripts 03-04) --------------------------------
INT_START <- as.Date("2024-11-01")

COL_HIST  <- "#4A90C4"
COL_INT   <- "#E07B32"
FILL_HIST <- "#DDEEF8"
FILL_INT  <- "#FEF0DC"
COL_CF    <- "#7B9E87"    # muted green for counterfactual -- distinct from both periods


# -----------------------------------------------------------------------------
# 1. LOAD SCRIPT 04 OUTPUTS
# -----------------------------------------------------------------------------

cat("Loading Script 04 outputs...\n")

required_files <- c(
  "04a_its_coefficients.csv",
  "04b_its_monthly_fitted_primary.csv",
  "04c_its_monthly_fitted_sensitivity.csv"
)

for (f in required_files) {
  path <- file.path(OUTPUT_DIR, f)
  if (!file.exists(path)) {
    stop(sprintf(
      "Required file not found: %s\nRun 04_its_analysis.R first.", path
    ))
  }
}

coefs    <- read_csv(file.path(OUTPUT_DIR, "04a_its_coefficients.csv"),
                     show_col_types = FALSE)
d_prim   <- read_csv(file.path(OUTPUT_DIR, "04b_its_monthly_fitted_primary.csv"),
                     show_col_types = FALSE)
d_sens   <- read_csv(file.path(OUTPUT_DIR, "04c_its_monthly_fitted_sensitivity.csv"),
                     show_col_types = FALSE)

cat(sprintf("  Coefficients     : %d rows\n", nrow(coefs)))
cat(sprintf("  Primary monthly  : %d months\n", nrow(d_prim)))
cat(sprintf("  Sensitivity monthly: %d months\n", nrow(d_sens)))

# Coerce date column
d_prim <- d_prim %>% mutate(adm_month_date = as.Date(adm_month_date))
d_sens <- d_sens %>% mutate(adm_month_date = as.Date(adm_month_date))


# -----------------------------------------------------------------------------
# 2. RE-FIT ITS MODELS TO OBTAIN CI BANDS
# -----------------------------------------------------------------------------
# predict() with interval = "confidence" requires the fitted lm object.
# The models are re-fitted here from the saved monthly data (which contains
# pct_normal, time, intervention, time_after, n_temp_recorded) -- identical
# to the models in Script 04.
#
# The 95% CI band on the fitted ITS line uses OLS standard errors (correct
# for visualising the uncertainty in the trend parameters). Newey-West SEs
# are used for inference (Panel B forest plot) but are not appropriate for
# constructing point-wise CI ribbons because they are asymptotic corrections
# to the coefficient SEs, not direct estimates of the variance of y-hat.
#
# Counterfactual CI is derived analytically from the model vcov matrix:
#   CF_i = b0 + b1 * time_i    (intervention = 0, time_after = 0)
#   Var(CF_i) = Var(b0) + time_i^2 * Var(b1) + 2 * time_i * Cov(b0, b1)
# where b0 and b1 are the (Intercept) and time coefficients of the full
# ITS model (not a separately fitted pre-intervention model).

cat("\nRe-fitting ITS models for CI computation...\n")

fit_model_from_csv <- function(d) {
  lm(pct_normal ~ time + intervention + time_after,
     data    = d,
     weights = n_temp_recorded)
}

model_prim <- fit_model_from_csv(d_prim)
model_sens <- fit_model_from_csv(d_sens)

cat(sprintf("  Primary    R2 = %.3f,  sigma = %.3f\n",
            summary(model_prim)$r.squared, summary(model_prim)$sigma))
cat(sprintf("  Sensitivity R2 = %.3f,  sigma = %.3f\n",
            summary(model_sens)$r.squared, summary(model_sens)$sigma))


# -- Attach CI bands to primary monthly data ----------------------------------

add_ci_bands <- function(d, model) {

  # 95% CI on fitted ITS line (OLS, point-wise)
  ci_fit  <- predict(model, interval = "confidence", level = 0.95)

  # Counterfactual and its CI from the full model vcov
  vcv     <- vcov(model)
  b0      <- coef(model)[["(Intercept)"]]
  b1      <- coef(model)[["time"]]
  v_b0    <- vcv["(Intercept)", "(Intercept)"]
  v_b1    <- vcv["time", "time"]
  c_b0b1  <- vcv["(Intercept)", "time"]
  cf_se   <- sqrt(v_b0 + d$time^2 * v_b1 + 2 * d$time * c_b0b1)

  d %>%
    mutate(
      fit_lower = ci_fit[, "lwr"],
      fit_upper = ci_fit[, "upr"],
      cf_lower  = counterfactual - 1.96 * cf_se,
      cf_upper  = counterfactual + 1.96 * cf_se
    )
}

d_prim <- add_ci_bands(d_prim, model_prim)
d_sens <- add_ci_bands(d_sens, model_sens)

cat("  CI bands attached to monthly data.\n")


# -----------------------------------------------------------------------------
# 3. PREPARE FOREST PLOT DATA
# -----------------------------------------------------------------------------
# Use Newey-West CIs for inference (as recommended in Script 04).
# Show b2 (level change) and b3 (slope change) for both analysis windows.

forest_raw <- coefs %>%
  filter(
    term    %in% c("intervention", "time_after"),
    se_type == "Newey-West (lag=2)"
  ) %>%
  mutate(
    term_label = case_when(
      term == "intervention" ~
        "Immediate level change at Nov 2024\n(b2, %-pts)",
      term == "time_after"   ~
        "Change in monthly slope after Nov 2024\n(b3, %-pts / month)"
    ),
    window_label = case_when(
      grepl("primary",     analysis, ignore.case = TRUE) ~ "Primary\n(Jan 2022 - Oct 2025)",
      grepl("sensitivity", analysis, ignore.case = TRUE) ~ "Sensitivity\n(Jan 2023 - Oct 2025)"
    ),
    # Ordered factor so primary appears at top within each facet
    window_fct = factor(window_label,
                        levels = c("Sensitivity\n(Jan 2023 - Oct 2025)",
                                   "Primary\n(Jan 2022 - Oct 2025)"))
  )

cat(sprintf("\nForest plot: %d rows (2 terms x 2 windows)\n", nrow(forest_raw)))

# Report NW estimates for the log
for (i in seq_len(nrow(forest_raw))) {
  cat(sprintf("  [%s] %s: %.2f (95%% CI %.2f to %.2f, p=%.4f)\n",
              forest_raw$window_label[i],
              forest_raw$term[i],
              forest_raw$estimate[i],
              forest_raw$ci_lower[i],
              forest_raw$ci_upper[i],
              forest_raw$p_value[i]))
}


# -----------------------------------------------------------------------------
# 4. PANEL A: ENHANCED ITS TIME SERIES (PRIMARY ANALYSIS)
# -----------------------------------------------------------------------------

cat("\nBuilding Panel A (enhanced ITS time series)...\n")

# Split into pre / post segments for ribbon and line plotting
d_pre  <- d_prim %>% filter(intervention == 0)
d_post <- d_prim %>% filter(intervention == 1)
d_cf   <- d_prim %>% filter(intervention == 1)   # counterfactual only in post-period

hist_xmin <- min(d_prim$adm_month_date)
hist_xmax <- INT_START
int_xmin  <- INT_START
int_xmax  <- max(d_prim$adm_month_date) + days(31)

# Key month annotations (from brief):
#   Nov 2024 - Jan 2025: early high performance (~89-91%)
#   Apr 2025, Sep-Oct 2025: notable unexplained drops
# These are annotated with small bracket-style vertical segments rather than
# large text blocks to keep the figure uncluttered.
annotate_months <- tribble(
  ~date,                  ~label,                    ~ypos,  ~vjust,
  as.Date("2024-12-01"),  "Early high\nperformance",  93,     0,
  as.Date("2025-04-01"),  "Notable\ndrop",            60,     1,
  as.Date("2025-09-01"),  "Notable\ndrops",           60,     1
)

p_ts <- ggplot(d_prim, aes(x = adm_month_date)) +

  # -- Period background shading ----------------------------------------------
  annotate("rect",
           xmin = hist_xmin, xmax = hist_xmax,
           ymin = -Inf,      ymax = Inf,
           fill = FILL_HIST, alpha = 1) +
  annotate("rect",
           xmin = int_xmin, xmax = int_xmax,
           ymin = -Inf,     ymax = Inf,
           fill = FILL_INT, alpha = 1) +

  # -- Intervention onset line -----------------------------------------------
  geom_vline(xintercept = INT_START,
             linetype = "dashed", colour = "grey45", linewidth = 0.6) +

  # -- Counterfactual CI ribbon and line -------------------------------------
  geom_ribbon(data = d_cf,
              aes(ymin = cf_lower, ymax = cf_upper),
              fill = COL_CF, alpha = 0.18) +
  geom_line(data = d_cf,
            aes(y = counterfactual),
            colour = COL_CF, linetype = "dotted", linewidth = 1.0) +

  # -- Fitted ITS CI ribbons (pre and post) ----------------------------------
  geom_ribbon(data = d_pre,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_HIST, alpha = 0.20) +
  geom_ribbon(data = d_post,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_INT, alpha = 0.20) +

  # -- Fitted ITS lines (pre and post) ---------------------------------------
  geom_line(data = d_pre,
            aes(y = fitted_its),
            colour = COL_HIST, linewidth = 1.3) +
  geom_line(data = d_post,
            aes(y = fitted_its),
            colour = COL_INT, linewidth = 1.3) +

  # -- Observed data points --------------------------------------------------
  geom_point(aes(y      = pct_normal,
                 colour = period,
                 size   = n_temp_recorded),
             alpha = 0.80) +

  # -- Period labels ---------------------------------------------------------
  annotate("text",
           x = as.Date("2023-01-15"), y = 97,
           label  = "Historical period",
           colour = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x = as.Date("2025-04-01"), y = 97,
           label  = "Intervention\nperiod",
           colour = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +

  # -- Counterfactual label --------------------------------------------------
  annotate("text",
           x      = as.Date("2025-06-15"),
           y      = d_cf$counterfactual[
             which.min(abs(d_cf$adm_month_date - as.Date("2025-06-15")))
           ] - 4.5,
           label  = "Counterfactual\n(pre-trend projected)",
           colour = COL_CF, hjust = 0.5, size = 2.7, fontface = "italic") +

  # -- Scales ----------------------------------------------------------------
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

  # -- Labels ----------------------------------------------------------------
  labs(
    title    = "Monthly proportion of SMCH inborn babies with normal admission temperature",
    subtitle = paste0(
      "Segmented linear regression ITS (primary analysis: Jan 2022 - Oct 2025, 46 months). ",
      "Shaded bands = 95% CI on fitted trends."
    ),
    x       = NULL,
    y       = "% with normal admission temperature (36.5-37.5 C)",
    caption = paste0(
      "ITS model: pct_normal ~ time + intervention + time_after, weighted by monthly n temperatures recorded.\n",
      "Solid lines = weighted OLS fitted trend. Shaded ribbons = 95% OLS confidence interval on fitted mean.",
      " Dotted line = pre-intervention trend projected forward (counterfactual).\n",
      "Point size proportional to n temperature recordings. ",
      "Inference (Panel B) uses Newey-West HAC standard errors (lag = 2)."
    )
  ) +

  # -- Theme -----------------------------------------------------------------
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
  file.path(OUTPUT_DIR, "06a_its_timeseries.png"),
  plot = p_ts, width = 13, height = 6, dpi = 300, units = "in"
)
cat("Saved: 06a_its_timeseries.png\n")


# -----------------------------------------------------------------------------
# 5. PANEL B: COEFFICIENT FOREST PLOT
# -----------------------------------------------------------------------------
# Shows b2 (level change) and b3 (slope change) for primary and sensitivity
# analyses, with Newey-West 95% CIs, in separate facets (free x-scales).
# A vertical reference line at 0 shows the null hypothesis.

cat("\nBuilding Panel B (coefficient forest plot)...\n")

# Colour by window
win_colours <- c(
  "Primary\n(Jan 2022 - Oct 2025)"     = COL_HIST,
  "Sensitivity\n(Jan 2023 - Oct 2025)" = COL_INT
)

p_forest <- ggplot(
  forest_raw,
  aes(x = estimate, y = window_fct, colour = window_label)
) +

  # Null reference line
  geom_vline(xintercept = 0,
             linetype = "solid", colour = "grey75", linewidth = 0.6) +

  # CI error bars (geom_errorbar with orientation = "y" replaces geom_errorbarh)
  geom_errorbar(
    aes(xmin = ci_lower, xmax = ci_upper),
    width       = 0.25,
    linewidth   = 0.85,
    orientation = "y"
  ) +

  # Estimate point
  geom_point(size = 3.8) +

  # Numeric annotation: estimate (CI) to the right of the error bar
  geom_text(
    aes(
      x     = ci_upper,
      label = sprintf("%.2f (%.2f, %.2f)",
                      estimate, ci_lower, ci_upper)
    ),
    hjust  = -0.08,
    size   = 2.9,
    colour = "grey25"
  ) +

  facet_wrap(
    ~ term_label,
    scales = "free_x",
    ncol   = 1
  ) +

  scale_colour_manual(
    name   = "Analysis window",
    values = win_colours,
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +

  labs(
    title    = "ITS coefficient estimates: level change and slope change",
    subtitle = paste0(
      "Newey-West HAC 95% confidence intervals (lag = 2). ",
      "Primary window: Jan 2022 - Oct 2025. ",
      "Sensitivity window: Jan 2023 - Oct 2025."
    ),
    x       = "Estimate (percentage points)",
    y       = NULL,
    caption = paste0(
      "b2 = immediate step change in % normal temperature at intervention onset (November 2024).\n",
      "b3 = change in monthly slope (%-pts/month) after intervention onset; ",
      "post-intervention slope = pre-intervention trend + b3.\n",
      "Positive b2 indicates an immediate improvement; ",
      "positive b3 indicates an accelerating improvement over the intervention period."
    )
  ) +

  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11.5),
    plot.subtitle    = element_text(colour = "grey40", size = 8.8),
    plot.caption     = element_text(colour = "grey50", size = 7.2, hjust = 0),
    axis.text        = element_text(size = 9),
    axis.title.x     = element_text(size = 9.5),
    strip.text       = element_text(face = "bold", size = 9.5),
    strip.background = element_rect(fill = "grey93", colour = NA),
    legend.position  = "none",     # window already on y-axis
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "06b_its_forest.png"),
  plot = p_forest, width = 9, height = 6, dpi = 300, units = "in"
)
cat("Saved: 06b_its_forest.png\n")


# -----------------------------------------------------------------------------
# 6. COMBINED FIGURE (Panel A | Panel B)
# -----------------------------------------------------------------------------
# Attempts to combine the two panels side by side (2:1 width ratio) using
# cowplot. If cowplot is not installed, saves a brief note and moves on --
# the individual panels are always produced regardless.
#
# In the combined figure, Panel A sits on the left (time series) and
# Panel B on the right (forest plot), with shared panel labels A and B.

cat("\nAttempting combined figure (requires cowplot)...\n")

combined_saved <- tryCatch({

  library(cowplot)

  # Strip title/subtitle from Panel B for the combined figure to save space;
  # carry them only in the individual file.
  p_forest_bare <- p_forest +
    labs(title = NULL, subtitle = NULL)

  p_combined <- plot_grid(
    p_ts, p_forest_bare,
    ncol       = 2,
    rel_widths = c(2.1, 1),
    labels     = c("A", "B"),
    label_size = 13,
    label_fontface = "bold"
  )

  # Add a shared title
  title_row <- ggdraw() +
    draw_label(
      paste0("Interrupted time series analysis: normal admission temperature, ",
             "SMCH inborn neonates (Jan 2022 - Oct 2025)"),
      fontface = "bold",
      size     = 11.5,
      x        = 0.02,
      hjust    = 0
    )

  p_final <- plot_grid(
    title_row, p_combined,
    ncol        = 1,
    rel_heights = c(0.06, 1)
  )

  ggsave(
    file.path(OUTPUT_DIR, "06_its_combined.png"),
    plot = p_final, width = 18, height = 7, dpi = 300, units = "in"
  )
  cat("Saved: 06_its_combined.png\n")
  TRUE

}, error = function(e) {
  cat(sprintf(
    "  cowplot not available (%s).\n  Individual panels saved; combine manually if needed.\n",
    conditionMessage(e)
  ))
  FALSE
})


# -----------------------------------------------------------------------------
# 7. SENSITIVITY ANALYSIS ANNOTATED PLOT
# -----------------------------------------------------------------------------
# A parallel version of Panel A using the sensitivity window (Jan 2023),
# so that the two plots can be placed side by side or in supplementary
# material without rerunning the script.

cat("\nBuilding sensitivity time series plot...\n")

d_s_pre  <- d_sens %>% filter(intervention == 0)
d_s_post <- d_sens %>% filter(intervention == 1)
d_s_cf   <- d_sens %>% filter(intervention == 1)

p_ts_sens <- ggplot(d_sens, aes(x = adm_month_date)) +

  annotate("rect",
           xmin = min(d_sens$adm_month_date), xmax = INT_START,
           ymin = -Inf, ymax = Inf,
           fill = FILL_HIST, alpha = 1) +
  annotate("rect",
           xmin = INT_START, xmax = max(d_sens$adm_month_date) + days(31),
           ymin = -Inf, ymax = Inf,
           fill = FILL_INT, alpha = 1) +

  geom_vline(xintercept = INT_START,
             linetype = "dashed", colour = "grey45", linewidth = 0.6) +

  geom_ribbon(data = d_s_cf,
              aes(ymin = cf_lower, ymax = cf_upper),
              fill = COL_CF, alpha = 0.18) +
  geom_line(data = d_s_cf,
            aes(y = counterfactual),
            colour = COL_CF, linetype = "dotted", linewidth = 1.0) +

  geom_ribbon(data = d_s_pre,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_HIST, alpha = 0.20) +
  geom_ribbon(data = d_s_post,
              aes(ymin = fit_lower, ymax = fit_upper),
              fill = COL_INT, alpha = 0.20) +

  geom_line(data = d_s_pre,
            aes(y = fitted_its),
            colour = COL_HIST, linewidth = 1.3) +
  geom_line(data = d_s_post,
            aes(y = fitted_its),
            colour = COL_INT, linewidth = 1.3) +

  geom_point(aes(y = pct_normal, colour = period, size = n_temp_recorded),
             alpha = 0.80) +

  annotate("text",
           x = as.Date("2024-01-01"), y = 97,
           label  = "Historical period\n(Jan 2023 - Oct 2024)",
           colour = COL_HIST, hjust = 0.5, size = 3.2, fontface = "bold") +
  annotate("text",
           x = as.Date("2025-04-01"), y = 97,
           label  = "Intervention\nperiod",
           colour = COL_INT, hjust = 0.5, size = 3.2, fontface = "bold") +

  annotate("text",
           x      = as.Date("2025-06-15"),
           y      = d_s_cf$counterfactual[
             which.min(abs(d_s_cf$adm_month_date - as.Date("2025-06-15")))
           ] - 4.5,
           label  = "Counterfactual\n(pre-trend projected)",
           colour = COL_CF, hjust = 0.5, size = 2.7, fontface = "italic") +

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
    title    = "Monthly proportion of SMCH inborn babies with normal admission temperature",
    subtitle = paste0(
      "Sensitivity analysis (Jan 2023 - Oct 2025, 34 months; excludes 2022 anomalous period). ",
      "Shaded bands = 95% CI on fitted trends."
    ),
    x       = NULL,
    y       = "% with normal admission temperature (36.5-37.5 C)",
    caption = paste0(
      "ITS model: pct_normal ~ time + intervention + time_after, weighted by monthly n temperatures recorded.\n",
      "Solid lines = weighted OLS fitted trend. Shaded ribbons = 95% OLS confidence interval on fitted mean.",
      " Dotted line = counterfactual (pre-intervention trend extended)."
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
  file.path(OUTPUT_DIR, "06c_its_timeseries_sensitivity.png"),
  plot = p_ts_sens, width = 13, height = 6, dpi = 300, units = "in"
)
cat("Saved: 06c_its_timeseries_sensitivity.png\n")


# -----------------------------------------------------------------------------
# 8. SUMMARY
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  FIGURE OUTPUT SUMMARY\n")
cat("================================================================\n\n")
cat("  06a_its_timeseries.png           -- Panel A (primary, 13x6 in)\n")
cat("  06b_its_forest.png               -- Panel B (forest, 9x6 in)\n")
cat("  06c_its_timeseries_sensitivity.png -- Sensitivity time series\n")
if (combined_saved) {
  cat("  06_its_combined.png              -- Combined A+B (18x7 in)\n")
} else {
  cat("  06_its_combined.png              -- NOT saved (cowplot unavailable)\n")
}

cat("\n  Key estimates (Newey-West, primary analysis):\n")
prim_nw <- forest_raw %>%
  filter(grepl("Primary", window_label))
for (i in seq_len(nrow(prim_nw))) {
  pval_str <- if (prim_nw$p_value[i] < 0.001) "<0.001" else
              sprintf("%.3f", prim_nw$p_value[i])
  cat(sprintf("    %s: %.2f %%-pts (95%% CI %.2f to %.2f, p=%s)\n",
              prim_nw$term[i],
              prim_nw$estimate[i],
              prim_nw$ci_lower[i],
              prim_nw$ci_upper[i],
              pval_str))
}

cat(sprintf("\nRun completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=== Script 06 complete ===\n")

sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
