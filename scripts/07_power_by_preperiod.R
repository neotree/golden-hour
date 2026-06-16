# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 07: Power by Pre-Intervention Window Length
# =============================================================================
# Purpose:
#   Sweep across every possible pre-intervention window length (n_pre = 4 to
#   34 months, counting back from October 2024) and estimate the power of the
#   ITS segmented regression at each length, for both key estimands:
#     b2 -- immediate level change at November 2024 (%-pts)
#     b3 -- change in monthly slope after November 2024 (%-pts/month)
#
# Motivation:
#   Script 05 showed a counterintuitive result: the sensitivity window
#   (22 pre-months, Jan 2023) has substantially higher power than the primary
#   window (34 pre-months, Jan 2022), despite fewer time points. This is
#   because the anomalous 2022 temperature recordings inflate the residual
#   sigma estimate, which overwhelms the benefit of the additional data.
#   This script sweeps all possible start dates to identify where power peaks.
#
# Method:
#   For each n_pre in 4:max_pre:
#     1. Subset the monthly data to the last n_pre historical months
#        (always ending October 2024) plus all 12 intervention months.
#     2. Fit a pre-intervention linear trend to estimate sigma (residual SD).
#     3. Run N_SIM simulations at each effect size using the same weighted
#        ITS model and noise structure as Script 05.
#     4. Record power(b2), power(b3), and power(either).
#
#   Effect sizes tested (chosen to bracket the observed effects from Sc. 04):
#     b2: 3, 5, 7, 10 %-pts  (observed b2 approx 5 %-pts)
#     b3: 0.5, 1.0, 1.65, 2.0 %-pts/month  (observed b3 approx 1.65 %-pts/mo)
#
# Outputs (all to 03-OUTPUTS/):
#   07a_power_by_preperiod.csv    -- full results table
#   07b_power_b2_by_preperiod.png -- power curves for b2 vs n_pre
#   07c_power_b3_by_preperiod.png -- power curves for b3 vs n_pre
#   07_preperiod_log.txt          -- analysis log
#
# Run time: approximately 4-8 minutes.
#
# DSH note: script is ASCII-only (no non-ASCII characters anywhere).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

set.seed(20241101)

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

LOG_FILE <- file.path(OUTPUT_DIR, "07_preperiod_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

cat("Golden Hour -- Power by Pre-Intervention Window Length\n")
cat(sprintf("Run started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# -- Fixed study boundaries ---------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2025-10-31")
HIST_START <- as.Date("2022-01-01")   # earliest available historical data
HIST_END   <- as.Date("2024-10-31")

MIN_N_TEMP <- 5

# -- Simulation parameters ----------------------------------------------------
N_SIM  <- 500     # per cell; reduced from Script 05 for tractable runtime
ALPHA  <- 0.05

# Effect sizes to test (brackets the observed effects from Script 04)
B2_TEST <- c(3, 5, 7, 10)          # %-pts  (observed ~5)
B3_TEST <- c(0.5, 1.0, 1.65, 2.0) # %-pts/month  (observed ~1.65)

# Reference lines for the two windows already analysed in Scripts 04-05
N_PRE_PRIMARY     <- 34   # Jan 2022 - Oct 2024
N_PRE_SENSITIVITY <- 22   # Jan 2023 - Oct 2024

# -- Colours ------------------------------------------------------------------
COL_HIST <- "#4A90C4"
COL_INT  <- "#E07B32"


# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)
cat(sprintf("  Rows loaded: %d\n", nrow(raw)))

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

smch_inborn <- df %>%
  filter(
    facility       == "SMCH",
    inborn         == TRUE,
    is_readmission == FALSE | is.na(is_readmission),
    period %in% c("Intervention", "Historical")
  )

# Full monthly series (all historical + intervention months)
monthly_all <- smch_inborn %>%
  group_by(adm_month, period) %>%
  summarise(
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

# Separate historical and intervention months
hist_months <- monthly_all %>%
  filter(period == "Historical") %>%
  arrange(adm_month_date)

int_months <- monthly_all %>%
  filter(period == "Intervention") %>%
  arrange(adm_month_date)

n_hist_max <- nrow(hist_months)
n_int      <- nrow(int_months)

cat(sprintf("\nHistorical months available: %d (max pre-intervention window)\n",
            n_hist_max))
cat(sprintf("Intervention months fixed  : %d\n", n_int))
cat(sprintf("Pre-intervention start dates: %s to %s\n",
            format(tail(hist_months$adm_month_date, N_PRE_PRIMARY)[1], "%b %Y"),
            format(tail(hist_months$adm_month_date, 4)[1], "%b %Y")))


# -----------------------------------------------------------------------------
# 2. SIMULATION HELPER (same logic as Script 05)
# -----------------------------------------------------------------------------

run_power_sim <- function(its_data, b0, b1, b2, b3, sigma, n_sim) {

  n        <- nrow(its_data)
  tv       <- its_data$time
  iv       <- its_data$intervention
  tav      <- its_data$time_after
  wv       <- its_data$n_temp_recorded
  noise_sd <- sigma / sqrt(wv)
  mu       <- b0 + b1 * tv + b2 * iv + b3 * tav

  sig_b2 <- logical(n_sim)
  sig_b3 <- logical(n_sim)

  for (i in seq_len(n_sim)) {
    y_sim <- pmin(pmax(mu + rnorm(n, 0, noise_sd), 0), 100)
    fit   <- lm(y_sim ~ tv + iv + tav, weights = wv)
    cf    <- summary(fit)$coefficients
    sig_b2[i] <- cf["iv",  "Pr(>|t|)"] < ALPHA
    sig_b3[i] <- cf["tav", "Pr(>|t|)"] < ALPHA
  }

  c(power_b2 = mean(sig_b2), power_b3 = mean(sig_b3),
    power_any = mean(sig_b2 | sig_b3))
}


# -----------------------------------------------------------------------------
# 3. SWEEP ACROSS PRE-INTERVENTION WINDOW LENGTHS
# -----------------------------------------------------------------------------
# n_pre = 4 is the minimum so the pre-intervention linear trend model
# (2 parameters) has df >= 2 for a stable sigma estimate.

MIN_PRE <- 4
MAX_PRE <- n_hist_max

cat(sprintf("\nSweeping n_pre = %d to %d (%d windows) x %d b2 x %d b3 effect sizes...\n",
            MIN_PRE, MAX_PRE,
            MAX_PRE - MIN_PRE + 1,
            length(B2_TEST), length(B3_TEST)))
cat(sprintf("N_SIM = %d per cell. Estimated run time: 4-8 min.\n\n", N_SIM))

results <- list()
win_idx <- 0

for (n_pre in MIN_PRE:MAX_PRE) {

  win_idx <- win_idx + 1

  # Take the last n_pre historical months (always ending Oct 2024)
  pre_subset <- tail(hist_months, n_pre)
  start_date <- min(pre_subset$adm_month_date)

  # Combine with fixed intervention months and build ITS variables
  window_data <- bind_rows(pre_subset, int_months) %>%
    arrange(adm_month_date) %>%
    mutate(
      time         = row_number(),
      intervention = as.integer(period == "Intervention"),
      time_after   = cumsum(intervention)
    )

  # Estimate sigma from pre-intervention linear trend
  pre_only <- window_data %>% filter(intervention == 0)

  sigma_ok <- tryCatch({
    pre_model <- lm(pct_normal ~ time,
                    data    = pre_only,
                    weights = n_temp_recorded)
    sigma_val <- summary(pre_model)$sigma
    b0_val    <- coef(pre_model)[["(Intercept)"]]
    b1_val    <- coef(pre_model)[["time"]]
    list(sigma = sigma_val, b0 = b0_val, b1 = b1_val, ok = TRUE)
  }, error = function(e) list(ok = FALSE))

  if (!sigma_ok$ok || is.na(sigma_ok$sigma) || sigma_ok$sigma <= 0) {
    cat(sprintf("  n_pre=%2d: sigma estimation failed -- skipping.\n", n_pre))
    next
  }

  sigma <- sigma_ok$sigma
  b0    <- sigma_ok$b0
  b1    <- sigma_ok$b1

  if (win_idx %% 5 == 1 || n_pre %in% c(MIN_PRE, N_PRE_SENSITIVITY, N_PRE_PRIMARY)) {
    cat(sprintf("  n_pre=%2d  start=%s  sigma=%.1f\n",
                n_pre, format(start_date, "%b %Y"), sigma))
  }

  # Run simulations for each b2 effect size (b3 = 0)
  for (eff_b2 in B2_TEST) {
    res <- run_power_sim(window_data, b0, b1,
                         b2 = eff_b2, b3 = 0,
                         sigma = sigma, n_sim = N_SIM)
    results[[length(results) + 1]] <- data.frame(
      n_pre      = n_pre,
      start_date = start_date,
      sigma      = sigma,
      estimand   = "b2",
      effect     = eff_b2,
      effect_label = sprintf("b2 = %g %%-pts", eff_b2),
      power      = res[["power_b2"]]
    )
  }

  # Run simulations for each b3 effect size (b2 = 0)
  for (eff_b3 in B3_TEST) {
    res <- run_power_sim(window_data, b0, b1,
                         b2 = 0, b3 = eff_b3,
                         sigma = sigma, n_sim = N_SIM)
    results[[length(results) + 1]] <- data.frame(
      n_pre        = n_pre,
      start_date   = start_date,
      sigma        = sigma,
      estimand     = "b3",
      effect       = eff_b3,
      effect_label = sprintf("b3 = %g %%-pts/mo", eff_b3),
      power        = res[["power_b3"]]
    )
  }
}

power_sweep <- bind_rows(results)

write_csv(power_sweep,
          file.path(OUTPUT_DIR, "07a_power_by_preperiod.csv"))
cat("\nSaved: 07a_power_by_preperiod.csv\n")


# -----------------------------------------------------------------------------
# 4. SUMMARY TABLE (PRINTED TO LOG)
# -----------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("  SIGMA AND POWER SUMMARY BY PRE-INTERVENTION WINDOW LENGTH\n")
cat("================================================================\n\n")

# Sigma by n_pre (one row per window)
sigma_tbl <- power_sweep %>%
  distinct(n_pre, start_date, sigma) %>%
  arrange(n_pre)

# Wide power table: rows = n_pre, columns = effect size
wide_b2 <- power_sweep %>%
  filter(estimand == "b2") %>%
  select(n_pre, start_date, sigma, effect_label, power) %>%
  pivot_wider(names_from = effect_label, values_from = power) %>%
  arrange(n_pre)

wide_b3 <- power_sweep %>%
  filter(estimand == "b3") %>%
  select(n_pre, start_date, sigma, effect_label, power) %>%
  pivot_wider(names_from = effect_label, values_from = power) %>%
  arrange(n_pre)

# Print b2 table
cat("  Power for level change (b2), b3 = 0:\n\n")
b2_cols <- paste(sprintf("b2 = %g", B2_TEST), "%-pts")
cat(sprintf("  %-6s  %-10s  %-8s  %s\n",
            "n_pre", "Start",
            "sigma",
            paste(sprintf("%-10s", b2_cols), collapse = "")))
cat(sprintf("  %s\n", strrep("-", 6 + 2 + 10 + 2 + 8 + 2 + 10 * length(B2_TEST))))

for (i in seq_len(nrow(wide_b2))) {
  row      <- wide_b2[i, ]
  marker   <- if (row$n_pre == N_PRE_PRIMARY)     " <- primary"     else
              if (row$n_pre == N_PRE_SENSITIVITY) " <- sensitivity" else ""
  pow_strs <- sapply(b2_cols, function(col) {
    v <- row[[col]]
    if (is.null(v) || is.na(v)) "  --      " else
      sprintf("%-10s", sprintf("%.0f%%", v * 100))
  })
  cat(sprintf("  %-6d  %-10s  %-8.1f  %s%s\n",
              row$n_pre,
              format(as.Date(row$start_date), "%b %Y"),
              row$sigma,
              paste(pow_strs, collapse = ""),
              marker))
}

# Print b3 table
cat(sprintf("\n\n  Power for slope change (b3), b2 = 0:\n\n"))
b3_cols <- paste(sprintf("b3 = %g", B3_TEST), "%-pts/mo")
cat(sprintf("  %-6s  %-10s  %-8s  %s\n",
            "n_pre", "Start",
            "sigma",
            paste(sprintf("%-12s", b3_cols), collapse = "")))
cat(sprintf("  %s\n", strrep("-", 6 + 2 + 10 + 2 + 8 + 2 + 12 * length(B3_TEST))))

for (i in seq_len(nrow(wide_b3))) {
  row      <- wide_b3[i, ]
  marker   <- if (row$n_pre == N_PRE_PRIMARY)     " <- primary"     else
              if (row$n_pre == N_PRE_SENSITIVITY) " <- sensitivity" else ""
  pow_strs <- sapply(b3_cols, function(col) {
    v <- row[[col]]
    if (is.null(v) || is.na(v)) "  --         " else
      sprintf("%-12s", sprintf("%.0f%%", v * 100))
  })
  cat(sprintf("  %-6d  %-10s  %-8.1f  %s%s\n",
              row$n_pre,
              format(as.Date(row$start_date), "%b %Y"),
              row$sigma,
              paste(pow_strs, collapse = ""),
              marker))
}

# Peak power rows
cat("\n  Peak power by estimand and effect size:\n")
for (est in c("b2", "b3")) {
  for (eff in if (est == "b2") B2_TEST else B3_TEST) {
    peak_row <- power_sweep %>%
      filter(estimand == est, abs(effect - eff) < 0.01) %>%
      slice_max(power, n = 1, with_ties = FALSE)
    if (nrow(peak_row) == 1) {
      unit <- if (est == "b2") "%-pts" else "%-pts/mo"
      cat(sprintf("    %s = %g %s: peak power %.0f%% at n_pre = %d (start %s, sigma = %.1f)\n",
                  est, eff, unit,
                  peak_row$power * 100,
                  peak_row$n_pre,
                  format(as.Date(peak_row$start_date), "%b %Y"),
                  peak_row$sigma))
    }
  }
}


# -----------------------------------------------------------------------------
# 5. PLOTS
# -----------------------------------------------------------------------------

cat("\nBuilding plots...\n")

# Reference band data for the two existing windows
ref_lines <- data.frame(
  n_pre = c(N_PRE_SENSITIVITY, N_PRE_PRIMARY),
  label = c("Sensitivity\n(Jan 2023)", "Primary\n(Jan 2022)"),
  col   = c(COL_INT, COL_HIST)
)

# Ordered factor for effect labels (so legend reads top to bottom)
power_sweep <- power_sweep %>%
  mutate(
    effect_fct = factor(effect_label,
                        levels = unique(effect_label[order(estimand, -effect)]))
  )

# Colour palette: 4 shades per estimand
pal_b2 <- c("#08519C", "#3182BD", "#6BAED6", "#BDD7E7")
pal_b3 <- c("#A63603", "#E6550D", "#FD8D3C", "#FDD0A2")

# -- Helper: add reference lines and 80% power band --------------------------
add_refs <- function(p) {
  p +
    geom_hline(yintercept = 0.80, linetype = "dashed",
               colour = "grey55", linewidth = 0.5) +
    annotate("text", x = MAX_PRE, y = 0.825,
             label = "80% power", hjust = 1, size = 2.9, colour = "grey45") +
    geom_vline(xintercept = N_PRE_SENSITIVITY,
               linetype = "dotted", colour = COL_INT, linewidth = 0.8) +
    annotate("text", x = N_PRE_SENSITIVITY - 0.4, y = 0.97,
             label = "Sensitivity\n(Jan 2023)", hjust = 1, size = 2.8,
             colour = COL_INT, fontface = "italic") +
    geom_vline(xintercept = N_PRE_PRIMARY,
               linetype = "dotted", colour = COL_HIST, linewidth = 0.8) +
    annotate("text", x = N_PRE_PRIMARY + 0.4, y = 0.97,
             label = "Primary\n(Jan 2022)", hjust = 0, size = 2.8,
             colour = COL_HIST, fontface = "italic")
}

common_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11.5),
    plot.subtitle   = element_text(colour = "grey40", size = 8.8),
    plot.caption    = element_text(colour = "grey50", size = 7.2, hjust = 0),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


# -- Panel b2: power for level change -----------------------------------------

b2_data <- power_sweep %>% filter(estimand == "b2") %>% droplevels()

# Sigma trace (right-axis proxy: normalise sigma to [0,1] for overlay)
sigma_trace <- b2_data %>%
  distinct(n_pre, sigma) %>%
  mutate(sigma_scaled = sigma / max(sigma, na.rm = TRUE))

p_b2 <- ggplot(b2_data,
               aes(x = n_pre, y = power, colour = effect_fct,
                   group = effect_fct)) +

  # Sigma trend (scaled to y-axis, plotted first so it sits behind)
  geom_line(data  = sigma_trace,
            aes(x = n_pre, y = sigma_scaled * 0.6),
            inherit.aes = FALSE,
            colour = "grey70", linetype = "solid", linewidth = 0.7) +
  annotate("text",
           x = MAX_PRE, y = sigma_trace$sigma_scaled[nrow(sigma_trace)] * 0.6 + 0.03,
           label = "Residual sigma\n(scaled)", hjust = 1,
           size = 2.6, colour = "grey55", fontface = "italic") +

  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 1.8, alpha = 0.8) +

  scale_colour_manual(
    name   = "Level change (b2)",
    values = setNames(pal_b2, levels(b2_data$effect_fct))
  ) +
  scale_x_continuous(
    breaks = c(seq(4, MAX_PRE, by = 2)),
    sec.axis = sec_axis(
      ~ .,
      name   = "Pre-intervention start month",
      breaks = c(N_PRE_SENSITIVITY, N_PRE_PRIMARY,
                 seq(4, MAX_PRE, by = 6)),
      labels = function(x) {
        dates <- as.Date(hist_months$adm_month_date)
        idx   <- nrow(hist_months) - x + 1
        idx   <- pmax(1, pmin(idx, nrow(hist_months)))
        format(dates[idx], "%b %Y")
      }
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +

  labs(
    title    = "Power to detect a level change (b2) by pre-intervention window length",
    subtitle = paste0(
      "Post-intervention period fixed at 12 months (Nov 2024 - Oct 2025). ",
      "Slope change b3 = 0.\n",
      "Grey line = scaled residual sigma (right scale); ",
      "vertical dotted lines mark the two existing analysis windows."
    ),
    x       = "Number of pre-intervention months",
    y       = "Power (alpha = 0.05, two-sided)",
    caption = paste0(
      "N_SIM = ", N_SIM, " per cell. Pre-intervention months counted back from October 2024.\n",
      "Observed b2 from Script 04 (primary window, Newey-West): 4.95 %-pts."
    )
  ) +
  common_theme

p_b2 <- add_refs(p_b2)

ggsave(
  file.path(OUTPUT_DIR, "07b_power_b2_by_preperiod.png"),
  plot = p_b2, width = 13, height = 6, dpi = 300, units = "in"
)
cat("Saved: 07b_power_b2_by_preperiod.png\n")


# -- Panel b3: power for slope change -----------------------------------------

b3_data <- power_sweep %>% filter(estimand == "b3") %>% droplevels()

p_b3 <- ggplot(b3_data,
               aes(x = n_pre, y = power, colour = effect_fct,
                   group = effect_fct)) +

  geom_line(data  = sigma_trace,
            aes(x = n_pre, y = sigma_scaled * 0.6),
            inherit.aes = FALSE,
            colour = "grey70", linetype = "solid", linewidth = 0.7) +
  annotate("text",
           x = MAX_PRE, y = sigma_trace$sigma_scaled[nrow(sigma_trace)] * 0.6 + 0.03,
           label = "Residual sigma\n(scaled)", hjust = 1,
           size = 2.6, colour = "grey55", fontface = "italic") +

  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 1.8, alpha = 0.8) +

  scale_colour_manual(
    name   = "Slope change (b3)",
    values = setNames(pal_b3, levels(b3_data$effect_fct))
  ) +
  scale_x_continuous(
    breaks = c(seq(4, MAX_PRE, by = 2)),
    sec.axis = sec_axis(
      ~ .,
      name   = "Pre-intervention start month",
      breaks = c(N_PRE_SENSITIVITY, N_PRE_PRIMARY,
                 seq(4, MAX_PRE, by = 6)),
      labels = function(x) {
        dates <- as.Date(hist_months$adm_month_date)
        idx   <- nrow(hist_months) - x + 1
        idx   <- pmax(1, pmin(idx, nrow(hist_months)))
        format(dates[idx], "%b %Y")
      }
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +

  labs(
    title    = "Power to detect a slope change (b3) by pre-intervention window length",
    subtitle = paste0(
      "Post-intervention period fixed at 12 months (Nov 2024 - Oct 2025). ",
      "Level change b2 = 0.\n",
      "Grey line = scaled residual sigma (right scale); ",
      "vertical dotted lines mark the two existing analysis windows."
    ),
    x       = "Number of pre-intervention months",
    y       = "Power (alpha = 0.05, two-sided)",
    caption = paste0(
      "N_SIM = ", N_SIM, " per cell. Pre-intervention months counted back from October 2024.\n",
      "Observed b3 from Script 04 (primary window, Newey-West): 1.65 %-pts/month."
    )
  ) +
  common_theme

p_b3 <- add_refs(p_b3)

ggsave(
  file.path(OUTPUT_DIR, "07c_power_b3_by_preperiod.png"),
  plot = p_b3, width = 13, height = 6, dpi = 300, units = "in"
)
cat("Saved: 07c_power_b3_by_preperiod.png\n")


# -- Sigma profile plot -------------------------------------------------------

sigma_profile <- power_sweep %>%
  distinct(n_pre, start_date, sigma) %>%
  arrange(n_pre)

p_sigma <- ggplot(sigma_profile, aes(x = n_pre, y = sigma)) +

  geom_line(colour = "grey50", linewidth = 1.0) +
  geom_point(colour = "grey35", size = 2.0) +

  geom_vline(xintercept = N_PRE_SENSITIVITY,
             linetype = "dotted", colour = COL_INT, linewidth = 0.8) +
  annotate("text", x = N_PRE_SENSITIVITY - 0.4, y = max(sigma_profile$sigma) * 0.97,
           label = "Sensitivity (Jan 2023)", hjust = 1, size = 2.8,
           colour = COL_INT, fontface = "italic") +
  geom_vline(xintercept = N_PRE_PRIMARY,
             linetype = "dotted", colour = COL_HIST, linewidth = 0.8) +
  annotate("text", x = N_PRE_PRIMARY + 0.4, y = max(sigma_profile$sigma) * 0.97,
           label = "Primary (Jan 2022)", hjust = 0, size = 2.8,
           colour = COL_HIST, fontface = "italic") +

  scale_x_continuous(breaks = seq(4, MAX_PRE, by = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +

  labs(
    title    = "Estimated residual sigma by pre-intervention window length",
    subtitle = paste0(
      "Sigma estimated from a weighted linear trend model fit to the pre-intervention months only.\n",
      "A jump in sigma as n_pre crosses 22 indicates the 2022 anomalous period is being included."
    ),
    x       = "Number of pre-intervention months",
    y       = "Residual sigma (%-pts, weighted OLS)",
    caption = paste0(
      "Sigma represents the baseline month-to-month variability in normal-temperature rates,\n",
      "after removing the secular pre-intervention trend. Higher sigma = lower power."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(colour = "grey40", size = 8.8),
    plot.caption  = element_text(colour = "grey50", size = 7.2, hjust = 0),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "07d_sigma_by_preperiod.png"),
  plot = p_sigma, width = 10, height = 5, dpi = 300, units = "in"
)
cat("Saved: 07d_sigma_by_preperiod.png\n")


cat(sprintf("\nRun completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=== Script 07 complete ===\n")

sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
