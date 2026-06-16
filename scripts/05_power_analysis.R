# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 05: ITS Power Analysis (simulation-based)
# =============================================================================
# Purpose:
#   Estimate the statistical power of the segmented regression ITS (Script 04)
#   to detect clinically meaningful changes in the monthly proportion of SMCH
#   inborn babies with normal admission temperature (36.5-37.5 C, WHO).
#
# Two estimands:
#   b2  -- immediate level change at November 2024 (%-pts, step function)
#   b3  -- change in monthly slope after intervention onset (%-pts/month)
#
# Design:
#   Primary window   : 46 months (Jan 2022 - Oct 2025; 34 pre, 12 post)
#   Sensitivity window: 34 months (Jan 2023 - Oct 2025; 22 pre, 12 post)
#
# Method:
#   1. Load the actual monthly series (same pipeline as Script 04).
#   2. Fit a pre-intervention linear trend to estimate residual sigma --
#      the month-to-month variability in normal-temp rates unexplained by
#      the secular trend. This parameterises the simulation noise.
#   3. For each candidate effect size on a grid, simulate N_SIM monthly
#      datasets using:
#        y_sim_i = (b0 + b1*t_i + b2*intv_i + b3*ta_i) + epsilon_i
#        epsilon_i ~ N(0, sigma^2 / w_i)
#      where w_i is the observed monthly temperature recording count (the
#      ITS weights). Outcomes are clamped to [0, 100] to stay on the % scale.
#   4. Fit the weighted ITS model to each simulated dataset; record whether
#      b2 (and/or b3) is significant at alpha = 0.05 (two-sided).
#   5. Power = proportion of significant results across N_SIM simulations.
#
# Autocorrelation sensitivity (b2 curve, primary window only):
#   The baseline simulation assumes independent residuals (AR(0)), giving an
#   upper bound on power. An additional sensitivity run simulates AR(1)
#   residuals (rho = 0.2, 0.4) to illustrate the impact of serial correlation.
#
# Clinical reference points:
#   A level change (b2) of >=5 %-pts or a slope change (b3) of >=0.5
#   %-pts/month would be clinically meaningful in this setting.
#   If Script 04 has already been run, the observed b2 and b3 from that
#   analysis are read from 04a_its_coefficients.csv and annotated on the plots.
#
# Outputs (all to 03-OUTPUTS/):
#   05a_power_b2.csv              -- b2 power table (all windows/rho)
#   05b_power_b3.csv              -- b3 power table (both windows)
#   05c_power_heatmap_primary.csv -- joint power grid (primary window)
#   05d_power_level_change.png    -- power curves for b2
#   05e_power_slope_change.png    -- power curves for b3
#   05f_power_heatmap.png         -- joint power heatmap (primary window)
#   05_power_log.txt              -- analysis log
#
# Run time: approximately 3-6 minutes depending on hardware.
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

set.seed(20241101)   # reproducible; seed = intervention start date

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

# -- Log file -----------------------------------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "05_power_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour ITS Power Analysis -- SMCH Zimbabwe\n")
cat(sprintf("Run started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Study period boundaries (consistent with Script 04) ----------------------
INT_START       <- as.Date("2024-11-01")
INT_END         <- as.Date("2025-10-31")
HIST_START      <- as.Date("2022-01-01")
HIST_START_SENS <- as.Date("2023-01-01")
HIST_END        <- as.Date("2024-10-31")

MIN_N_TEMP <- 5      # minimum monthly temperature recordings for inclusion

# -- Simulation parameters ----------------------------------------------------
N_SIM  <- 1000       # simulations per effect size cell
ALPHA  <- 0.05       # two-sided significance level

# Grid for level change (b2): 0 to 20 %-pts in steps of 0.5
B2_GRID <- seq(0, 20, by = 0.5)

# Grid for slope change (b3): 0 to 2 %-pts/month in steps of 0.1
B3_GRID <- seq(0, 2, by = 0.1)

# Coarser joint grid for the heatmap (smaller N_SIM to keep runtime reasonable)
B2_HEAT   <- seq(0, 15, by = 1)
B3_HEAT   <- seq(0, 1.5, by = 0.1)
N_SIM_HEAT <- 500

# AR(1) rho values for the autocorrelation sensitivity section
RHO_VALS <- c(0, 0.2, 0.4)

# -- Colours (consistent with Scripts 03-04) ----------------------------------
COL_PRIM <- "#4A90C4"
COL_SENS <- "#E07B32"


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
# 2. CLEANING AND DERIVED VARIABLES (identical to Script 04)
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
# 3. ANALYSIS POPULATION AND MONTHLY DATA (identical to Script 04)
# -----------------------------------------------------------------------------

smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission)
  )

q10_pop <- smch_inborn %>%
  filter(period %in% c("Intervention", "Historical"))

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
  mutate(adm_month_date = as.Date(adm_month))

cat(sprintf("\nMonthly data: %d months total (%d historical, %d intervention)\n",
            nrow(monthly_raw),
            sum(monthly_raw$period == "Historical"),
            sum(monthly_raw$period == "Intervention")))


# -----------------------------------------------------------------------------
# 4. BUILD ITS DESIGN MATRICES FOR BOTH WINDOWS
# -----------------------------------------------------------------------------
# The simulation uses the actual observed monthly structure (same time points,
# same weights) so that power estimates reflect the real study design.

build_its_vars <- function(monthly_data) {
  monthly_data %>%
    arrange(adm_month) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )
}

its_prim <- build_its_vars(monthly_raw)
its_sens <- build_its_vars(
  monthly_raw %>% filter(as.Date(adm_month) >= HIST_START_SENS)
)

cat(sprintf("\nPrimary    window: %d months (%d pre, %d post)\n",
            nrow(its_prim),
            sum(its_prim$intervention == 0),
            sum(its_prim$intervention == 1)))
cat(sprintf("Sensitivity window: %d months (%d pre, %d post)\n",
            nrow(its_sens),
            sum(its_sens$intervention == 0),
            sum(its_sens$intervention == 1)))


# -----------------------------------------------------------------------------
# 5. ESTIMATE RESIDUAL SIGMA FROM PRE-INTERVENTION TREND
# -----------------------------------------------------------------------------
# Fit a simple linear trend to the pre-intervention months only (no intervention
# dummy variables). The residual SD from this model -- sigma -- captures the
# month-to-month variability in temperature rates above the secular trend, and
# is used to parameterise the simulation noise in all subsequent sections.
#
# In the weighted lm, summary()$sigma is:
#   sqrt( sum(w_i * resid_i^2) / df )
# This is the residual SD on the outcome scale (%-pts), accounting for the
# heteroscedasticity induced by the weights.

cat("\n--- Estimating residual sigma from pre-intervention data ---\n")

pre_prim <- its_prim %>% filter(intervention == 0)
pre_sens <- its_sens %>% filter(intervention == 0)

pre_model_prim <- lm(pct_normal ~ time,
                     data = pre_prim, weights = n_temp_recorded)
pre_model_sens <- lm(pct_normal ~ time,
                     data = pre_sens, weights = n_temp_recorded)

sigma_prim <- summary(pre_model_prim)$sigma
sigma_sens <- summary(pre_model_sens)$sigma

b0_prim <- coef(pre_model_prim)[["(Intercept)"]]
b1_prim <- coef(pre_model_prim)[["time"]]
b0_sens <- coef(pre_model_sens)[["(Intercept)"]]
b1_sens <- coef(pre_model_sens)[["time"]]

cat(sprintf("  Primary    sigma: %.3f %%-pts\n", sigma_prim))
cat(sprintf("  Sensitivity sigma: %.3f %%-pts\n", sigma_sens))
cat(sprintf("  Primary    pre-trend: %.2f + %.4f * t  (%.2f %%-pts/yr)\n",
            b0_prim, b1_prim, b1_prim * 12))
cat(sprintf("  Sensitivity pre-trend: %.2f + %.4f * t  (%.2f %%-pts/yr)\n",
            b0_sens, b1_sens, b1_sens * 12))


# -----------------------------------------------------------------------------
# 6. SIMULATION HELPER FUNCTION
# -----------------------------------------------------------------------------
#
# run_power_sim()
#
# For a given ITS design and noise level, simulates N_SIM monthly datasets
# under an assumed level change (b2) and slope change (b3), fits the weighted
# ITS model each time, and returns the power for b2, b3, and either.
#
# Arguments:
#   its_data -- data frame with columns: time, intervention, time_after,
#               n_temp_recorded (used as weights). pct_normal is NOT used;
#               the expected mu is constructed from b0, b1, b2, b3.
#   b0, b1   -- pre-intervention baseline from the pre-intervention model
#   b2       -- assumed level change at intervention onset (%-pts)
#   b3       -- assumed slope change after intervention onset (%-pts/month)
#   sigma    -- residual SD from the pre-intervention model
#   n_sim    -- number of simulated datasets
#   rho      -- AR(1) coefficient for serial correlation (0 = independent)
#
# Returns a named numeric vector with three elements:
#   power_b2  -- P(b2 significant)
#   power_b3  -- P(b3 significant)
#   power_any -- P(b2 OR b3 significant)

run_power_sim <- function(its_data, b0, b1, b2, b3, sigma, n_sim, rho = 0) {

  n   <- nrow(its_data)
  tv  <- its_data$time
  iv  <- its_data$intervention
  tav <- its_data$time_after
  wv  <- its_data$n_temp_recorded

  # Expected outcome under the assumed effect sizes
  mu <- b0 + b1 * tv + b2 * iv + b3 * tav

  # Noise SD per month: sigma / sqrt(weight) under the weighted-regression model
  noise_sd <- sigma / sqrt(wv)

  sig_b2  <- logical(n_sim)
  sig_b3  <- logical(n_sim)

  for (i in seq_len(n_sim)) {

    # Simulate residuals (independent or AR(1))
    if (rho == 0) {
      eps <- rnorm(n, 0, noise_sd)
    } else {
      # AR(1) with marginal SD = noise_sd:
      # e_t = rho * e_{t-1} + sqrt(1 - rho^2) * z_t  (unit-variance AR(1))
      # eps_t = noise_sd_t * e_t
      z    <- rnorm(n)
      e    <- numeric(n)
      e[1] <- z[1]
      for (k in 2:n) {
        e[k] <- rho * e[k - 1] + sqrt(1 - rho^2) * z[k]
      }
      eps <- noise_sd * e
    }

    y_sim <- pmin(pmax(mu + eps, 0), 100)   # clamp to [0, 100]

    fit <- lm(y_sim ~ tv + iv + tav, weights = wv)
    cf  <- summary(fit)$coefficients

    sig_b2[i] <- cf["iv",  "Pr(>|t|)"] < ALPHA
    sig_b3[i] <- cf["tav", "Pr(>|t|)"] < ALPHA
  }

  c(
    power_b2  = mean(sig_b2),
    power_b3  = mean(sig_b3),
    power_any = mean(sig_b2 | sig_b3)
  )
}


# -----------------------------------------------------------------------------
# 7. POWER CURVES FOR LEVEL CHANGE (b2)
# -----------------------------------------------------------------------------
# Slope change fixed at b3 = 0. Runs for:
#   (a) Primary window, independent residuals (AR rho = 0)
#   (b) Sensitivity window, independent residuals
#   (c) Primary window, AR(1) at rho = 0.2 and 0.4 (autocorrelation sensitivity)

cat("\n=== Power curves: level change (b2), b3 = 0 ===\n")

power_b2_rows <- list()

for (rho in RHO_VALS) {
  cat(sprintf("  rho = %.1f  [primary window]...\n", rho))
  for (eff in B2_GRID) {
    res <- run_power_sim(
      its_data = its_prim,
      b0 = b0_prim, b1 = b1_prim,
      b2 = eff, b3 = 0,
      sigma = sigma_prim, n_sim = N_SIM, rho = rho
    )
    power_b2_rows[[length(power_b2_rows) + 1]] <- data.frame(
      window    = "Primary (Jan 2022)",
      rho       = rho,
      b2        = eff,
      power_b2  = res[["power_b2"]],
      power_b3  = res[["power_b3"]],
      power_any = res[["power_any"]]
    )
  }
}

cat("  rho = 0.0  [sensitivity window]...\n")
for (eff in B2_GRID) {
  res <- run_power_sim(
    its_data = its_sens,
    b0 = b0_sens, b1 = b1_sens,
    b2 = eff, b3 = 0,
    sigma = sigma_sens, n_sim = N_SIM, rho = 0
  )
  power_b2_rows[[length(power_b2_rows) + 1]] <- data.frame(
    window    = "Sensitivity (Jan 2023)",
    rho       = 0,
    b2        = eff,
    power_b2  = res[["power_b2"]],
    power_b3  = res[["power_b3"]],
    power_any = res[["power_any"]]
  )
}

power_b2 <- bind_rows(power_b2_rows)

write_csv(power_b2, file.path(OUTPUT_DIR, "05a_power_b2.csv"))
cat("  Saved: 05a_power_b2.csv\n")


# -----------------------------------------------------------------------------
# 8. POWER CURVES FOR SLOPE CHANGE (b3)
# -----------------------------------------------------------------------------
# Level change fixed at b2 = 0.  Both windows, independent residuals only.

cat("\n=== Power curves: slope change (b3), b2 = 0 ===\n")

power_b3_rows <- list()

for (win_label in c("Primary (Jan 2022)", "Sensitivity (Jan 2023)")) {
  use_prim <- win_label == "Primary (Jan 2022)"
  cat(sprintf("  %s...\n", win_label))
  for (eff in B3_GRID) {
    res <- run_power_sim(
      its_data = if (use_prim) its_prim else its_sens,
      b0       = if (use_prim) b0_prim  else b0_sens,
      b1       = if (use_prim) b1_prim  else b1_sens,
      b2 = 0, b3 = eff,
      sigma    = if (use_prim) sigma_prim else sigma_sens,
      n_sim = N_SIM, rho = 0
    )
    power_b3_rows[[length(power_b3_rows) + 1]] <- data.frame(
      window   = win_label,
      b3       = eff,
      power_b3 = res[["power_b3"]]
    )
  }
}

power_b3 <- bind_rows(power_b3_rows)

write_csv(power_b3, file.path(OUTPUT_DIR, "05b_power_b3.csv"))
cat("  Saved: 05b_power_b3.csv\n")


# -----------------------------------------------------------------------------
# 9. JOINT POWER HEATMAP (primary window)
# -----------------------------------------------------------------------------
# Power to detect b2 OR b3 significant, over a grid of both effect sizes.
# Uses a reduced N_SIM_HEAT for tractable run time.

cat(sprintf("\n=== Joint power heatmap (primary window): %d x %d grid, %d sims/cell ===\n",
            length(B2_HEAT), length(B3_HEAT), N_SIM_HEAT))

heat_grid   <- expand.grid(b2 = B2_HEAT, b3 = B3_HEAT)
n_cells     <- nrow(heat_grid)
power_any_v <- numeric(n_cells)

for (ci in seq_len(n_cells)) {
  if (ci %% 20 == 0) cat(sprintf("  Cell %d / %d\n", ci, n_cells))
  res <- run_power_sim(
    its_data = its_prim,
    b0 = b0_prim, b1 = b1_prim,
    b2 = heat_grid$b2[ci],
    b3 = heat_grid$b3[ci],
    sigma = sigma_prim, n_sim = N_SIM_HEAT, rho = 0
  )
  power_any_v[ci] <- res[["power_any"]]
}

power_heat <- heat_grid %>% mutate(power_any = power_any_v)

write_csv(power_heat, file.path(OUTPUT_DIR, "05c_power_heatmap_primary.csv"))
cat("  Saved: 05c_power_heatmap_primary.csv\n")


# -----------------------------------------------------------------------------
# 10. CHECK FOR OBSERVED EFFECT SIZES FROM SCRIPT 04
# -----------------------------------------------------------------------------
# If 04a_its_coefficients.csv exists, extract the observed b2 and b3 from
# the primary Newey-West model for annotation on the power plots.

b2_obs <- NA_real_
b3_obs <- NA_real_

coef_csv <- file.path(OUTPUT_DIR, "04a_its_coefficients.csv")
if (file.exists(coef_csv)) {
  cat("\nReading Script 04 coefficient results for plot annotation...\n")
  coefs04 <- read_csv(coef_csv, show_col_types = FALSE)
  prim_nw <- coefs04 %>%
    filter(grepl("primary", analysis, ignore.case = TRUE),
           se_type == "Newey-West (lag=2)")
  b2_row <- prim_nw %>% filter(term == "intervention")
  b3_row <- prim_nw %>% filter(term == "time_after")
  if (nrow(b2_row) == 1) {
    b2_obs <- abs(b2_row$estimate)   # absolute value for power curve x-axis
    cat(sprintf("  Observed b2 (level change): %.2f %%-pts\n", b2_row$estimate))
  }
  if (nrow(b3_row) == 1) {
    b3_obs <- abs(b3_row$estimate)
    cat(sprintf("  Observed b3 (slope change): %.2f %%-pts/month\n", b3_row$estimate))
  }
} else {
  cat("\nScript 04 output not found -- run 04_its_analysis.R to annotate plots.\n")
}


# -----------------------------------------------------------------------------
# 11. PLOTS
# -----------------------------------------------------------------------------

cat("\nBuilding plots...\n")

# Helper: 80% power reference line and label
ref_line <- function(xmax) {
  list(
    geom_hline(yintercept = 0.80, linetype = "dashed",
               colour = "grey60", linewidth = 0.5),
    annotate("text", x = xmax * 0.97, y = 0.825,
             label = "80% power", colour = "grey50",
             hjust = 1, size = 3)
  )
}

# -- 11a. Power curves for level change (b2) ----------------------------------

# Subset the main curves (independent residuals only) for the primary plot layer
b2_main <- power_b2 %>% filter(rho == 0)

# AR(1) sensitivity curves (primary window, rho > 0)
b2_ar <- power_b2 %>%
  filter(window == "Primary (Jan 2022)", rho > 0) %>%
  mutate(ar_label = sprintf("Primary, AR(1) rho=%.1f", rho))

p_b2 <- ggplot() +

  ref_line(max(B2_GRID)) +

  # AR(1) sensitivity (dashed, behind main curves)
  geom_line(
    data = b2_ar,
    aes(x = b2, y = power_b2, group = ar_label, linetype = ar_label),
    colour = COL_PRIM, linewidth = 0.75, alpha = 0.65
  ) +

  # Main independent-residual curves (both windows)
  geom_line(
    data = b2_main,
    aes(x = b2, y = power_b2, colour = window),
    linewidth = 1.2
  ) +
  geom_point(
    data = b2_main,
    aes(x = b2, y = power_b2, colour = window),
    size = 1.5, alpha = 0.8
  ) +

  # Observed b2 from Script 04 (if available)
  {
    if (!is.na(b2_obs)) {
      list(
        geom_vline(xintercept = b2_obs, linetype = "dotdash",
                   colour = "grey30", linewidth = 0.6),
        annotate("text", x = b2_obs + 0.3, y = 0.12,
                 label = sprintf("Observed\nb2 = %.1f", b2_obs),
                 colour = "grey30", hjust = 0, size = 2.8)
      )
    }
  } +

  scale_colour_manual(
    name   = "Window",
    values = c("Primary (Jan 2022)" = COL_PRIM,
               "Sensitivity (Jan 2023)" = COL_SENS)
  ) +
  scale_linetype_manual(
    name   = "AR(1) sensitivity\n(primary window)",
    values = c("Primary, AR(1) rho=0.2" = "dashed",
               "Primary, AR(1) rho=0.4" = "dotted")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 20, by = 2),
    labels = function(x) paste0(x, " %pts")
  ) +
  labs(
    title    = "Power to detect a level change (b2) at November 2024",
    subtitle = paste0(
      "Solid lines: independent residuals. ",
      "Dashed/dotted: primary window with AR(1) (rho = 0.2, 0.4). ",
      "Slope change b3 = 0.\n",
      "Alpha = 0.05, N_SIM = ", N_SIM, " per cell."
    ),
    x = "Assumed level change (%-pts, absolute)",
    y = "Power"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11.5),
    plot.subtitle   = element_text(colour = "grey40", size = 8.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "05d_power_level_change.png"),
  plot = p_b2, width = 11, height = 5.5, dpi = 300, units = "in"
)
cat("Saved: 05d_power_level_change.png\n")


# -- 11b. Power curves for slope change (b3) ----------------------------------

p_b3 <- ggplot(power_b3, aes(x = b3, y = power_b3, colour = window)) +

  ref_line(max(B3_GRID)) +

  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.8) +

  {
    if (!is.na(b3_obs)) {
      list(
        geom_vline(xintercept = b3_obs, linetype = "dotdash",
                   colour = "grey30", linewidth = 0.6),
        annotate("text", x = b3_obs + 0.03, y = 0.12,
                 label = sprintf("Observed\nb3 = %.2f", b3_obs),
                 colour = "grey30", hjust = 0, size = 2.8)
      )
    }
  } +

  scale_colour_manual(
    name   = "Window",
    values = c("Primary (Jan 2022)" = COL_PRIM,
               "Sensitivity (Jan 2023)" = COL_SENS)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 2, by = 0.2),
    labels = function(x) paste0(x, " %pts/mo")
  ) +
  labs(
    title    = "Power to detect a slope change (b3) after November 2024",
    subtitle = paste0(
      "Level change b2 = 0. Independent residuals assumed. ",
      "Alpha = 0.05, N_SIM = ", N_SIM, " per cell."
    ),
    x = "Assumed slope change (%-pts/month, absolute)",
    y = "Power"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11.5),
    plot.subtitle   = element_text(colour = "grey40", size = 8.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "05e_power_slope_change.png"),
  plot = p_b3, width = 11, height = 5.5, dpi = 300, units = "in"
)
cat("Saved: 05e_power_slope_change.png\n")


# -- 11c. Joint power heatmap (primary window) --------------------------------

p_heat <- ggplot(power_heat, aes(x = b2, y = b3, fill = power_any)) +

  geom_tile(colour = "white", linewidth = 0.15) +

  # 80% power contour
  geom_contour(aes(z = power_any), breaks = 0.80,
               colour = "white", linewidth = 1.0, linetype = "dashed") +

  # Annotate observed effect sizes if available
  {
    if (!is.na(b2_obs) && !is.na(b3_obs)) {
      list(
        geom_point(
          data = data.frame(b2 = b2_obs, b3 = b3_obs),
          aes(x = b2, y = b3), inherit.aes = FALSE,
          shape = 4, size = 5, colour = "white", stroke = 1.5
        ),
        annotate("text",
                 x = b2_obs + 0.5, y = b3_obs + 0.05,
                 label = "Observed\neffect",
                 colour = "white", hjust = 0, size = 2.8)
      )
    }
  } +

  scale_fill_gradient2(
    name     = "Power\n(b2 OR b3)",
    low      = "#F7FBFF",
    mid      = "#6BAED6",
    high     = "#08306B",
    midpoint = 0.5,
    limits   = c(0, 1),
    labels   = percent_format(accuracy = 1)
  ) +
  scale_x_continuous(
    name   = "Level change b2 (%-pts)",
    breaks = B2_HEAT[B2_HEAT %% 3 == 0]
  ) +
  scale_y_continuous(
    name   = "Slope change b3 (%-pts/month)",
    breaks = round(B3_HEAT[round(B3_HEAT * 10) %% 3 == 0], 1)
  ) +
  labs(
    title    = "Joint power: primary ITS window (Jan 2022 - Oct 2025)",
    subtitle = paste0(
      "Power to detect b2 OR b3 significant at alpha = 0.05. ",
      "Dashed white contour = 80% power. N_SIM = ", N_SIM_HEAT, " per cell."
    ),
    caption  = paste0(
      "46 months (34 pre-intervention, 12 post-intervention). ",
      "Weights = observed monthly temperature recording counts. ",
      "Independent residuals assumed.\n",
      if (!is.na(b2_obs) && !is.na(b3_obs))
        sprintf("Cross (x) = observed effect from Script 04 (b2=%.1f, b3=%.2f).", b2_obs, b3_obs)
      else
        "Run 04_its_analysis.R to annotate observed effect size."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(colour = "grey40", size = 8.5),
    plot.caption  = element_text(colour = "grey50", size = 7.5, hjust = 0),
    panel.grid    = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "05f_power_heatmap.png"),
  plot = p_heat, width = 10, height = 6.5, dpi = 300, units = "in"
)
cat("Saved: 05f_power_heatmap.png\n")


# -----------------------------------------------------------------------------
# 12. SUMMARY TABLE
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  POWER ANALYSIS SUMMARY\n")
cat("================================================================\n\n")

cat(sprintf("  Pre-intervention residual sigma:\n"))
cat(sprintf("    Primary    window : %.2f %%-pts\n", sigma_prim))
cat(sprintf("    Sensitivity window: %.2f %%-pts\n\n", sigma_sens))

cat(sprintf("  N_SIM = %d  |  alpha = %.2f (two-sided)\n\n", N_SIM, ALPHA))

# -- b2 power table -----------------------------------------------------------
cat("  Power to detect level change (b2), b3 = 0:\n")
cat(sprintf("    %-8s  %-24s  %-24s  %-18s\n",
            "b2 (%-pts)", "Primary rho=0", "Primary rho=0.2",
            "Sensitivity rho=0"))
cat(sprintf("    %s\n", strrep("-", 78)))
for (eff in c(2, 3, 5, 7, 10, 15)) {
  r_p0 <- power_b2 %>%
    filter(window == "Primary (Jan 2022)",    rho == 0,   abs(b2 - eff) < 0.01) %>%
    pull(power_b2)
  r_p2 <- power_b2 %>%
    filter(window == "Primary (Jan 2022)",    rho == 0.2, abs(b2 - eff) < 0.01) %>%
    pull(power_b2)
  r_s0 <- power_b2 %>%
    filter(window == "Sensitivity (Jan 2023)", rho == 0,  abs(b2 - eff) < 0.01) %>%
    pull(power_b2)
  fmt <- function(x) if (length(x) == 1) sprintf("%.0f%%", x * 100) else "  --"
  cat(sprintf("    %-10.0f  %-24s  %-24s  %s\n",
              eff, fmt(r_p0), fmt(r_p2), fmt(r_s0)))
}

# -- b3 power table -----------------------------------------------------------
cat(sprintf("\n  Power to detect slope change (b3), b2 = 0:\n"))
cat(sprintf("    %-15s  %-20s  %s\n",
            "b3 (%-pts/mo)", "Primary (rho=0)", "Sensitivity (rho=0)"))
cat(sprintf("    %s\n", strrep("-", 56)))
for (eff in c(0.3, 0.5, 1.0, 1.5, 2.0)) {
  r_p <- power_b3 %>%
    filter(window == "Primary (Jan 2022)",    abs(b3 - eff) < 0.01) %>%
    pull(power_b3)
  r_s <- power_b3 %>%
    filter(window == "Sensitivity (Jan 2023)", abs(b3 - eff) < 0.01) %>%
    pull(power_b3)
  fmt <- function(x) if (length(x) == 1) sprintf("%.0f%%", x * 100) else "  --"
  cat(sprintf("    %-17.1f  %-20s  %s\n", eff, fmt(r_p), fmt(r_s)))
}

# -- 80% power thresholds -----------------------------------------------------
cat("\n  Minimum detectable effect at 80% power (approximate):\n")

mde_b2_prim <- power_b2 %>%
  filter(window == "Primary (Jan 2022)", rho == 0) %>%
  filter(power_b2 >= 0.80) %>%
  slice_min(b2, n = 1) %>%
  pull(b2)
mde_b2_sens <- power_b2 %>%
  filter(window == "Sensitivity (Jan 2023)", rho == 0) %>%
  filter(power_b2 >= 0.80) %>%
  slice_min(b2, n = 1) %>%
  pull(b2)
mde_b3_prim <- power_b3 %>%
  filter(window == "Primary (Jan 2022)") %>%
  filter(power_b3 >= 0.80) %>%
  slice_min(b3, n = 1) %>%
  pull(b3)
mde_b3_sens <- power_b3 %>%
  filter(window == "Sensitivity (Jan 2023)") %>%
  filter(power_b3 >= 0.80) %>%
  slice_min(b3, n = 1) %>%
  pull(b3)

fmt_mde <- function(x, unit) {
  if (length(x) == 1) sprintf("%.1f %s", x, unit) else ">grid max"
}
cat(sprintf("    Level change b2 -- Primary   : %s\n",
            fmt_mde(mde_b2_prim, "%-pts")))
cat(sprintf("    Level change b2 -- Sensitivity: %s\n",
            fmt_mde(mde_b2_sens, "%-pts")))
cat(sprintf("    Slope change b3 -- Primary   : %s\n",
            fmt_mde(mde_b3_prim, "%-pts/month")))
cat(sprintf("    Slope change b3 -- Sensitivity: %s\n",
            fmt_mde(mde_b3_sens, "%-pts/month")))

if (!is.na(b2_obs))
  cat(sprintf("\n  Observed b2 from Script 04: %.2f %%-pts\n", b2_obs))
if (!is.na(b3_obs))
  cat(sprintf("  Observed b3 from Script 04: %.2f %%-pts/month\n", b3_obs))

cat(sprintf("\nRun completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=== Script 05 complete ===\n")

# -- Close log ----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
