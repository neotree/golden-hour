# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 03: Admission Temperature Trends (Q10)
# =============================================================================
# Purpose:
#   Plot the monthly percentage of babies with normal admission temperature
#   (36.5-37.5 C, WHO definition) from January 2022 to October 2025.
#
# Definitions & methodology:
#   Data source: ZIM_db_master_joined_to_20260525.csv.
#   Shades the historical period (Jan 2022 - Oct 2024) and the intervention
#   period (Nov 2024 - Oct 2025) in distinct colours; shows individual
#   monthly data points sized by number of observations; overlays a loess
#   smoothed trend line for each period.
#
# Outputs (to 03-OUTPUTS/):
#   03a_monthly_normal_temp.csv  -- underlying monthly data
#   03_temperature_trends.png    -- the plot
#   03_temperature_trends_log.txt
#
# DSH note: script is ASCII-only.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)    # for date_breaks, percent, comma

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
            "ZIM_db_master_joined_to_20260525.csv"),
  mustWork = FALSE
)

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -- Redirect all output to a log file ----------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "03_temperature_trends_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

# -- Study period boundaries --------------------------------------------------
INT_START  <- as.Date("2024-11-01")
INT_END    <- as.Date("2026-03-31")   # extended from Oct 2025 to Mar 2026
HIST_START <- as.Date("2022-01-01")
HIST_END   <- as.Date("2024-10-31")

# -- Plot appearance constants ------------------------------------------------
# Colours: muted blue for historical, amber for intervention
COL_HIST   <- "#4A90C4"   # steel blue
COL_INT    <- "#E07B32"   # burnt amber
FILL_HIST  <- "#DDEEF8"   # very light blue
FILL_INT   <- "#FEF0DC"   # very light amber
MIN_N_TEMP <- 5           # exclude months with fewer than this many temp records


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

cat("Loading data...\n")
# guess_max = 100000: sample up to 100k rows for type-guessing so that
# mixed-type columns whose problem values fall after row 1000 are handled
# correctly rather than silently coerced to NA with a parse warning.
# DATE-PARSING FIX (Aug 2026, applied here to match Scripts 16/17/18/26):
# read_csv's type guesser converts "datetimeadmission" to POSIXct; ymd_hms()
# then re-coerces it via as.character(), which for a midnight timestamp
# yields a date-only string ("2025-08-08") that ymd_hms cannot parse -> NA ->
# the record is silently dropped. Read the column as character and take the
# date part directly instead.
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000,
                 col_types = cols(datetimeadmission = col_character()))

# Report any residual parse problems so they are visible in the log
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
# Reproduces the derivations from 01_descriptive_overview.R so this script
# is fully self-contained.

df <- raw %>%
  mutate(
    adm_date  = as.Date(substr(datetimeadmission, 1, 10)),
    adm_month = floor_date(adm_date, "month"),

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
# 3. DEFINE ANALYSIS POPULATION (Q10)
# -----------------------------------------------------------------------------

all_smch    <- df %>% filter(facility == "SMCH")
smch_inborn <- all_smch %>%
  filter(inborn == TRUE, is_readmission == FALSE | is.na(is_readmission))

# Q10 population: combined intervention + historical, SMCH inborn
q10_pop <- smch_inborn %>%
  filter(period %in% c("Intervention", "Historical"))

cat(sprintf("\nQ10 population (SMCH inborn, Jan 2022 - Mar 2026): %d records\n",
            nrow(q10_pop)))
cat(sprintf("  Historical  : %d\n", sum(q10_pop$period == "Historical")))
cat(sprintf("  Intervention: %d\n", sum(q10_pop$period == "Intervention")))
cat(sprintf("  Temperature recorded: %d / %d (%.1f%%)\n",
            sum(!is.na(q10_pop$temperature_num)),
            nrow(q10_pop),
            100 * mean(!is.na(q10_pop$temperature_num))))


# -----------------------------------------------------------------------------
# 4. BUILD MONTHLY SUMMARY TABLE
# -----------------------------------------------------------------------------

monthly_temp <- q10_pop %>%
  group_by(adm_month, period) %>%
  summarise(
    n_admissions     = n(),
    n_temp_recorded  = sum(!is.na(temperature_num)),
    n_normal         = sum(temp_normal == TRUE, na.rm = TRUE),
    n_hypothermic    = sum(temperature_num < 36.5, na.rm = TRUE),
    n_hyperthermic   = sum(temperature_num > 37.5, na.rm = TRUE),
    mean_temp        = round(mean(temperature_num, na.rm = TRUE), 2),
    pct_normal       = if_else(
      n_temp_recorded > 0,
      round(100 * n_normal / n_temp_recorded, 1),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(adm_month)

cat("\nMonthly normal temperature summary:\n")
print(monthly_temp, n = Inf)

write_csv(monthly_temp,
          file.path(OUTPUT_DIR, "03a_monthly_normal_temp.csv"))
cat("\nSaved: 03a_monthly_normal_temp.csv\n")

# Subset for plotting: exclude months with too few temperature recordings
plot_data <- monthly_temp %>%
  filter(n_temp_recorded >= MIN_N_TEMP)

cat(sprintf(
  "\nMonths excluded from plot (n_temp < %d): %d\n",
  MIN_N_TEMP,
  nrow(monthly_temp) - nrow(plot_data)
))
cat(sprintf("Months retained for plot: %d\n", nrow(plot_data)))


# -----------------------------------------------------------------------------
# 5. BUILD THE PLOT
# -----------------------------------------------------------------------------

cat("\nBuilding temperature trend plot...\n")

# -- Period shading as rectangles (one per period) ----------------------------
# xmin/xmax are the first day of the first and last months (+1 month end buffer)
hist_xmin <- HIST_START
hist_xmax <- HIST_END + days(31)    # end of October 2024
int_xmin  <- INT_START
int_xmax  <- INT_END  + days(31)    # end of March 2026

# -- Overall y-range for annotation placement ---------------------------------
y_max <- max(plot_data$pct_normal, na.rm = TRUE)
y_min <- min(plot_data$pct_normal, na.rm = TRUE)
y_label_pos <- max(plot_data$pct_normal, na.rm = TRUE) + 2   # just above data

# -- Factor for consistent period colours/legend order -------------------------
plot_data <- plot_data %>%
  mutate(period = factor(period, levels = c("Historical", "Intervention")))

p <- ggplot(plot_data, aes(x = adm_month, y = pct_normal)) +

  # -- Period shading --------------------------------------------------------
  annotate("rect",
           xmin  = hist_xmin, xmax = hist_xmax,
           ymin  = -Inf,      ymax = Inf,
           fill  = FILL_HIST, alpha = 1) +
  annotate("rect",
           xmin  = int_xmin, xmax = int_xmax,
           ymin  = -Inf,     ymax = Inf,
           fill  = FILL_INT, alpha = 1) +

  # -- Vertical line at intervention start ------------------------------------
  geom_vline(xintercept = INT_START,
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +

  # -- Data points (size proportional to n_temp_recorded) --------------------
  geom_point(aes(colour = period, size = n_temp_recorded),
             alpha = 0.85) +

  # -- Loess smoothed trend line per period -----------------------------------
  geom_smooth(aes(colour = period, fill = period, group = period),
              method  = "loess",
              formula = y ~ x,
              span    = 0.75,
              se      = TRUE,
              alpha   = 0.15,
              linewidth = 1.1) +

  # -- Period label annotations -----------------------------------------------
  annotate("text",
           x     = as.Date("2023-01-01"),
           y     = y_label_pos,
           label = "Historical period",
           colour = COL_HIST, hjust = 0.5, size = 3.2,
           fontface = "bold") +
  annotate("text",
           x     = as.Date("2025-07-01"),
           y     = y_label_pos,
           label = "Intervention\nperiod",
           colour = COL_INT, hjust = 0.5, size = 3.2,
           fontface = "bold") +

  # -- Scales and axes --------------------------------------------------------
  scale_colour_manual(
    name   = "Period",
    values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
  ) +
  scale_fill_manual(
    name   = "Period",
    values = c("Historical" = COL_HIST, "Intervention" = COL_INT)
  ) +
  scale_size_continuous(
    name   = "Temps\nrecorded",
    range  = c(1.5, 5),
    breaks = c(50, 100, 200, 400)
  ) +
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b\n%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.10))
  ) +

  # -- Labels -----------------------------------------------------------------
  labs(
    title    = "Monthly percentage of SMCH admissions with normal admission temperature",
    subtitle = "Normal temperature defined as 36.5-37.5 C (WHO). SMCH inborn babies, Jan 2022 - Mar 2026.",
    x        = NULL,
    y        = "% with normal admission temperature",
    caption  = paste0(
      "Loess smoothed trend lines shown per period (span = 0.75, 95% CI shaded).",
      "\nDashed vertical line marks intervention start (November 2024).",
      "\nMonths with fewer than ", MIN_N_TEMP, " temperature recordings excluded.",
      "\nPoint size proportional to number of temperatures recorded."
    )
  ) +

  # -- Theme ------------------------------------------------------------------
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 9.5),
    plot.caption     = element_text(colour = "grey50", size = 8, hjust = 0),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 9),
    axis.title.y     = element_text(size = 10),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


# -----------------------------------------------------------------------------
# 6. SAVE PLOT
# -----------------------------------------------------------------------------

PNG_PATH <- file.path(OUTPUT_DIR, "03_temperature_trends.png")

ggsave(
  filename = PNG_PATH,
  plot     = p,
  width    = 12,
  height   = 5.5,
  dpi      = 300,
  units    = "in"
)
cat(sprintf("Saved: 03_temperature_trends.png  (%dx%d px @ 300 dpi)\n",
            12 * 300, round(5.5 * 300)))


# -----------------------------------------------------------------------------
# 7. DESCRIPTIVE SUMMARY FOR THE MANUSCRIPT
# -----------------------------------------------------------------------------

cat("\n=== DESCRIPTIVE SUMMARY FOR Q10 ===\n")

# Overall normal temp rate by period
cat("\nNormal temperature rate by period (months with >=", MIN_N_TEMP, "recordings):\n")
plot_data %>%
  group_by(period) %>%
  summarise(
    n_months        = n(),
    total_records   = sum(n_temp_recorded),
    total_normal    = sum(n_normal),
    overall_pct_normal = round(100 * total_normal / total_records, 1),
    median_monthly_pct = round(median(pct_normal, na.rm = TRUE), 1),
    iqr_lo          = round(quantile(pct_normal, 0.25, na.rm = TRUE), 1),
    iqr_hi          = round(quantile(pct_normal, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(factor(period, levels = c("Historical", "Intervention"))) %>%
  print()

cat("\n=== Script 03 complete ===\n")

# -- Close log -----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
