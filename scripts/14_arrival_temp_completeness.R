# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 14: Completeness of arrival-to-neonatal-unit TEMPERATURE variables
# =============================================================================
# Purpose (Action 4, 2026-06-09 meeting):
#   Quantify the recording completeness of the "temperature on arrival to the
#   neonatal unit" variable(s), to decide whether a delivery-room-to-NNU
#   temperature analysis is feasible. The meeting recollection was that this
#   variable was introduced partway through the project (~Jan/Feb 2025) and had
#   poor completeness. Decision: include only if completeness is adequate;
#   otherwise leave out (it is a "bonus", not a paper priority).
#
#   This script reports, overall and BY MONTH, the proportion of SMCH inborn
#   admissions with each candidate temperature field populated, and applies a
#   simple usability rule of thumb (>= 60% complete across the intervention
#   period to be analysable; >= 80% to be robust).
#
# CANDIDATE VARIABLES (master NeoTree file):
#   nnuadmtemp  -- NNU admission temperature (closest to "arrival to unit")
#   tempnnu     -- alternative NNU temperature field
#   tempadm     -- admission temperature (delivery-room / early)
#   tempbirth   -- temperature at birth
#   temperature -- primary admission temperature (reference; ~99% complete)
#   dischtemp   -- discharge temperature (reference)
#
#   The script auto-detects which of these columns exist in the file, so it is
#   robust to schema changes.
#
# OUTPUTS (all to 03-OUTPUTS/):
#   14a_arrival_temp_completeness_overall.csv  -- one row per candidate variable
#   14b_arrival_temp_completeness_monthly.csv  -- variable x month completeness
#   14c_arrival_temp_completeness_plot.png     -- monthly completeness lines
#   14_arrival_temp_completeness_log.txt
#
# DSH note: ASCII-only.
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
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag),
                                      mustWork = FALSE))
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

LOG_FILE <- file.path(OUTPUT_DIR, "14_arrival_temp_completeness_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2026-03-31")

# Usability thresholds (% complete over the intervention period)
THRESH_ANALYSABLE <- 60
THRESH_ROBUST     <- 80

CANDIDATES <- c("nnuadmtemp", "tempnnu", "tempadm", "tempbirth",
                "temperature", "dischtemp")

cat("=============================================================\n")
cat("  Script 14: Arrival-to-NNU temperature completeness\n")
cat("=============================================================\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat(sprintf("  Data   : %s\n", DATA_PATH))
cat("=============================================================\n\n")


# -----------------------------------------------------------------------------
# 1. LOAD
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)

df <- raw %>%
  mutate(
    adm_dt   = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date = as.Date(adm_dt),
    month    = format(floor_date(adm_dt, "month"), "%Y-%m"),
    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                             ~ NA
    ),
    is_readmission = readmission == "Y"
  )

pop <- df %>%
  filter(facility == "SMCH", inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission),
         adm_date >= INT_START, adm_date <= INT_END)

cat(sprintf("\nSMCH inborn intervention admissions: %d\n", nrow(pop)))

present <- CANDIDATES[CANDIDATES %in% names(pop)]
missing_cols <- setdiff(CANDIDATES, present)
if (length(missing_cols) > 0) {
  cat(sprintf("  NOTE: candidate column(s) not in file: %s\n",
              paste(missing_cols, collapse = ", ")))
}

# A value counts as "recorded" if non-NA, non-blank, and (when numeric-coercible)
# a plausible temperature 25-45 C; otherwise just non-blank.
is_recorded <- function(x) {
  chr <- trimws(as.character(x))
  nonblank <- !is.na(x) & chr != "" & toupper(chr) != "U"
  num <- suppressWarnings(as.numeric(chr))
  plausible <- is.na(num) | (num >= 25 & num <= 45)
  nonblank & plausible
}


# -----------------------------------------------------------------------------
# 2. OVERALL COMPLETENESS
# -----------------------------------------------------------------------------

overall <- map_dfr(present, function(v) {
  rec <- is_recorded(pop[[v]])
  n   <- nrow(pop); k <- sum(rec, na.rm = TRUE)
  pct <- round(100 * k / n, 1)
  verdict <- if (pct >= THRESH_ROBUST) "robust"
             else if (pct >= THRESH_ANALYSABLE) "analysable (caution)"
             else "too incomplete -- not analysable"
  tibble(variable = v, n_total = n, n_recorded = k,
         pct_complete = pct, verdict = verdict)
}) %>% arrange(desc(pct_complete))

cat("\n--- Overall completeness over intervention period ---\n")
print(as.data.frame(overall), row.names = FALSE)
write_csv(overall, file.path(OUTPUT_DIR, "14a_arrival_temp_completeness_overall.csv"))
cat("\nSaved: 14a_arrival_temp_completeness_overall.csv\n")


# -----------------------------------------------------------------------------
# 3. MONTHLY COMPLETENESS (when did each variable start being recorded?)
# -----------------------------------------------------------------------------

MONTHS <- format(seq(INT_START, INT_END, by = "month"), "%Y-%m")

monthly <- map_dfr(present, function(v) {
  map_dfr(MONTHS, function(m) {
    d <- pop %>% filter(month == m)
    n <- nrow(d); k <- sum(is_recorded(d[[v]]), na.rm = TRUE)
    tibble(variable = v, month = m, n_total = n, n_recorded = k,
           pct_complete = if (n > 0) round(100 * k / n, 1) else NA_real_)
  })
})

cat("\n--- Monthly completeness (wide, % complete) ---\n")
monthly_wide <- monthly %>%
  select(variable, month, pct_complete) %>%
  pivot_wider(names_from = variable, values_from = pct_complete)
print(as.data.frame(monthly_wide), row.names = FALSE)
write_csv(monthly, file.path(OUTPUT_DIR, "14b_arrival_temp_completeness_monthly.csv"))
cat("\nSaved: 14b_arrival_temp_completeness_monthly.csv\n")


# -----------------------------------------------------------------------------
# 4. PLOT
# -----------------------------------------------------------------------------

plot_df <- monthly %>% mutate(adm_month = ym(month))

p <- ggplot(plot_df, aes(x = adm_month, y = pct_complete, colour = variable)) +
  geom_hline(yintercept = THRESH_ANALYSABLE, linetype = "dashed",
             colour = "grey55") +
  annotate("text", x = ym(MONTHS[2]), y = THRESH_ANALYSABLE + 3,
           label = paste0(THRESH_ANALYSABLE, "% analysable threshold"),
           hjust = 0, size = 3, colour = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(name = "% recorded", limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_x_date(name = NULL, date_breaks = "2 months", date_labels = "%b %Y") +
  scale_colour_brewer(name = "Temperature field", palette = "Set1") +
  labs(
    title    = "Completeness of temperature fields by month",
    subtitle = paste0("SMCH inborn admissions, ", INT_START, " to ", INT_END,
                      "\nFocus: arrival-to-NNU fields (nnuadmtemp, tempnnu) vs reference (temperature)"),
    caption  = "A field is 'recorded' if non-blank and a plausible 25-45 C value."
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title  = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "14c_arrival_temp_completeness_plot.png"),
       plot = p, width = 10, height = 5.5, dpi = 150)
cat("Saved: 14c_arrival_temp_completeness_plot.png\n")


# -----------------------------------------------------------------------------
# 5. VERDICT
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  RECOMMENDATION\n")
cat("=============================================================\n")
arrival_vars <- intersect(c("nnuadmtemp", "tempnnu"), present)
best <- overall %>% filter(variable %in% arrival_vars) %>% slice_max(pct_complete, n = 1)
if (nrow(best) == 1) {
  cat(sprintf("  Best arrival-to-NNU field: %s (%.1f%% complete) -> %s\n",
              best$variable, best$pct_complete, best$verdict))
  if (best$pct_complete < THRESH_ANALYSABLE) {
    cat("  CONCLUSION: completeness is below the analysable threshold.\n")
    cat("  Consistent with the meeting decision to deprioritise / leave out.\n")
  } else {
    cat("  CONCLUSION: completeness may support a sensitivity analysis;\n")
    cat("  review monthly ramp-up before committing.\n")
  }
} else {
  cat("  No arrival-to-NNU temperature field found in the file.\n")
}

cat("\n=== Script 14 complete ===\n")
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
