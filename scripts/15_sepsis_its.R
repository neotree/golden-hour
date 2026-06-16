# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 15: Interrupted Time Series (ITS) -- Neonatal SEPSIS incidence
# =============================================================================
# Purpose (Action 5, 2026-06-09 meeting -- OPTIONAL / interpret with caution):
#   Apply the same segmented-regression ITS framework as Script 12 (mortality)
#   to monthly neonatal SEPSIS incidence, motivated by published evidence that
#   DCC and skin-to-skin reduce sepsis risk.
#
#   HEAVY CONFOUNDING CAVEAT (must accompany any reporting):
#   A concurrent quality-improvement programme (Tinashe's PhD: blood-culture
#   and sepsis-diagnosis improvement) ran over the SAME period. Any change in
#   recorded sepsis incidence cannot be cleanly attributed to the Golden Hour
#   intervention -- improved diagnostics could RAISE recorded sepsis even if
#   true incidence fell. Treat this analysis as exploratory only.
#
# OUTCOME:
#   Monthly sepsis rate = sepsis cases / admissions * 100.
#   Sepsis definition is configurable (see SEPSIS_DEF below). Default uses
#   'suspectedneonatalsepsis' (the only field with sufficient monthly counts;
#   confirmed 'seps'/'neonatalsepsis' are too rare for monthly ITS, ~2-3/mo).
#
# MODEL (identical to Scripts 04b / 12):
#   pct_sepsis ~ time + intervention + time_after
#   Weighted by monthly admissions. Ljung-Box check; Newey-West HAC SEs (lag=2)
#   if autocorrelation detected.
#
# POPULATION:
#   Master NeoTree file. facility == "SMCH", inborn == TRUE, non-readmission.
#   (Sepsis is an admission-side diagnosis -> discharge match NOT required.)
#   Strata: All inborn / LBW (<2500 g) / NBW (>=2500 g).
#   Window : primary Jan 2023 - Mar 2026 (22 pre + 17 post); 2022 excluded
#            (consistent with 04b/12). Sensitivity incl. 2022 also run.
#
# OUTPUTS (all to 03-OUTPUTS/):
#   15a_sepsis_monthly.csv        -- monthly sepsis rates (all + subgroups)
#   15b_sepsis_its_coefficients.csv
#   15c_sepsis_its_plot.png       -- all-inborn ITS
#   15_sepsis_its_log.txt
#
# DSH note: ASCII-only.
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

LOG_FILE <- file.path(OUTPUT_DIR, "15_sepsis_its_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

INT_START       <- as.Date("2024-11-01")
INT_END         <- as.Date("2026-03-31")
HIST_START      <- as.Date("2023-01-01")
HIST_END        <- as.Date("2024-10-31")
HIST_START_SENS <- as.Date("2022-01-01")
MIN_N_DENOM     <- 10

# Sepsis case definition. Options:
#   "suspected" -> suspectedneonatalsepsis populated (default; analysable)
#   "confirmed" -> seps OR neonatalsepsis populated (rare; likely too sparse)
#   "any"       -> suspected OR confirmed
SEPSIS_DEF <- "suspected"

COL_HIST  <- "#4A90C4"; COL_INT <- "#E07B32"
FILL_HIST <- "#DDEEF8"; FILL_INT <- "#FEF0DC"; COL_CF <- "#7B9E87"

cat("Golden Hour Sepsis ITS Analysis -- SMCH Zimbabwe\n")
cat(sprintf("Run started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Sepsis definition: %s\n\n", SEPSIS_DEF))


# -----------------------------------------------------------------------------
# 1. LOAD AND DERIVE
# -----------------------------------------------------------------------------

cat("Loading data...\n")
raw <- read_csv(DATA_PATH, show_col_types = FALSE, guess_max = 100000)
cat(sprintf("  Rows loaded: %d\n", nrow(raw)))

# Helper: a diagnosis field is "positive" if non-NA and non-blank.
pos <- function(x) !is.na(x) & trimws(as.character(x)) != "" &
                   toupper(trimws(as.character(x))) != "U"

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
    is_readmission = readmission == "Y",

    bw_num = suppressWarnings(as.numeric(birthweight)),
    lbw    = !is.na(bw_num) & bw_num < 2500,
    nbw    = !is.na(bw_num) & bw_num >= 2500,

    susp_sepsis = if ("suspectedneonatalsepsis" %in% names(.))
                    pos(suspectedneonatalsepsis) else FALSE,
    conf_sepsis = (if ("seps" %in% names(.)) pos(seps) else FALSE) |
                  (if ("neonatalsepsis" %in% names(.)) pos(neonatalsepsis) else FALSE),

    period_primary = case_when(
      adm_date >= INT_START  & adm_date <= INT_END   ~ "Intervention",
      adm_date >= HIST_START & adm_date <= HIST_END  ~ "Historical",
      TRUE                                            ~ "Other"
    ),
    period_sens = case_when(
      adm_date >= INT_START       & adm_date <= INT_END  ~ "Intervention",
      adm_date >= HIST_START_SENS & adm_date <= HIST_END ~ "Historical",
      TRUE                                                ~ "Other"
    )
  ) %>%
  mutate(
    sepsis_case = case_when(
      SEPSIS_DEF == "suspected" ~ susp_sepsis,
      SEPSIS_DEF == "confirmed" ~ conf_sepsis,
      SEPSIS_DEF == "any"       ~ susp_sepsis | conf_sepsis,
      TRUE                      ~ susp_sepsis
    )
  )

pop <- df %>%
  filter(facility == "SMCH", inborn == TRUE,
         is_readmission == FALSE | is.na(is_readmission))

cat(sprintf("\nSMCH inborn (all time): %d\n", nrow(pop)))
cat(sprintf("Sepsis cases (def=%s): %d (%.1f%% of admissions)\n",
            SEPSIS_DEF, sum(pop$sepsis_case, na.rm = TRUE),
            100 * mean(pop$sepsis_case, na.rm = TRUE)))


# -----------------------------------------------------------------------------
# 2. MONTHLY TABLES
# -----------------------------------------------------------------------------

build_monthly <- function(data, period_var) {
  data %>%
    group_by(adm_month, period = .data[[period_var]]) %>%
    summarise(
      n_admissions = n(),
      n_sepsis     = sum(sepsis_case, na.rm = TRUE),
      pct_sepsis   = if_else(n_admissions > 0,
                             round(100 * n_sepsis / n_admissions, 2),
                             NA_real_),
      .groups = "drop"
    ) %>%
    filter(period %in% c("Intervention", "Historical")) %>%
    arrange(adm_month) %>%
    mutate(adm_month_date = as.Date(adm_month))
}

q_primary <- pop %>% filter(period_primary %in% c("Intervention", "Historical"))
monthly_all <- build_monthly(q_primary, "period_primary")
monthly_lbw <- build_monthly(q_primary %>% filter(lbw == TRUE), "period_primary")
monthly_nbw <- build_monthly(q_primary %>% filter(nbw == TRUE), "period_primary")

cat("\n--- Monthly sepsis rate (all inborn, primary window) ---\n")
print(monthly_all %>% select(adm_month_date, period, n_admissions,
                             n_sepsis, pct_sepsis), n = Inf)

monthly_combined <- bind_rows(
  monthly_all %>% mutate(subgroup = "All inborn"),
  monthly_lbw %>% mutate(subgroup = "LBW (<2500g)"),
  monthly_nbw %>% mutate(subgroup = "NBW (>=2500g)")
) %>%
  select(subgroup, adm_month_date, period, n_admissions, n_sepsis, pct_sepsis)
write_csv(monthly_combined, file.path(OUTPUT_DIR, "15a_sepsis_monthly.csv"))
cat("\nSaved: 15a_sepsis_monthly.csv\n")


# -----------------------------------------------------------------------------
# 3. ITS FIT (same engine as Script 12)
# -----------------------------------------------------------------------------

fit_its <- function(monthly, label, min_n = MIN_N_DENOM) {
  d <- monthly %>%
    filter(n_admissions >= min_n) %>%
    arrange(adm_month) %>%
    mutate(time = row_number(),
           intervention = as.integer(period == "Intervention"),
           time_after = cumsum(intervention))

  n_pre <- sum(d$intervention == 0); n_post <- sum(d$intervention == 1)
  cat(sprintf("\n%s\n  ITS: %s\n  Pre months: %d  Int months: %d\n%s\n",
              strrep("-", 60), label, n_pre, n_post, strrep("-", 60)))
  if (n_pre < 3 || n_post < 2) {
    cat("  WARNING: too few months. Skipping.\n"); return(NULL)
  }

  model <- lm(pct_sepsis ~ time + intervention + time_after,
              data = d, weights = n_admissions)
  cat("\n--- OLS summary ---\n"); print(summary(model))

  lb <- Box.test(residuals(model), lag = 1, type = "Ljung-Box")
  cat(sprintf("\nLjung-Box (lag 1): Q=%.3f p=%.4f%s\n",
              lb$statistic, lb$p.value,
              if (lb$p.value < 0.05) "  >> Newey-West applied" else ""))
  nw <- coeftest(model, vcov = NeweyWest(model, lag = 2, prewhite = FALSE))
  cat("\n--- Newey-West (lag=2) ---\n"); print(nw)

  extract <- function(mat, se_type) {
    data.frame(analysis = label, term = rownames(mat),
               estimate = mat[, 1], std_error = mat[, 2],
               t_value = mat[, 3], p_value = mat[, 4],
               se_type = se_type, stringsAsFactors = FALSE) %>%
      mutate(ci_lower = estimate - 1.96 * std_error,
             ci_upper = estimate + 1.96 * std_error,
             across(where(is.numeric), ~ round(.x, 4))) %>%
      select(analysis, term, estimate, std_error, ci_lower, ci_upper,
             t_value, p_value, se_type)
  }
  cf <- coef(model)
  d <- d %>% mutate(fitted_its = fitted(model),
                    counterfactual = cf["(Intercept)"] + cf["time"] * time)
  list(tidy_ols = extract(summary(model)$coefficients, "OLS"),
       tidy_nw  = extract(as.matrix(nw), "Newey-West (lag=2)"),
       lb = lb, d = d, b0 = cf["(Intercept)"], b1 = cf["time"],
       b2 = cf["intervention"], b3 = cf["time_after"])
}

cat("\n", strrep("=", 64), "\n  FITTING SEPSIS ITS MODELS\n", strrep("=", 64), "\n", sep = "")
its_all <- fit_its(monthly_all, "All inborn -- primary (Jan 2023 - Mar 2026)")
its_lbw <- fit_its(monthly_lbw, "LBW (<2500g) -- primary")
its_nbw <- fit_its(monthly_nbw, "NBW (>=2500g) -- primary")

q_sens <- pop %>% filter(period_sens %in% c("Intervention", "Historical"))
its_sens <- fit_its(build_monthly(q_sens, "period_sens"),
                    "All inborn -- sensitivity (Jan 2022 - Mar 2026)")

all_coefs <- bind_rows(
  if (!is.null(its_all))  bind_rows(its_all$tidy_ols,  its_all$tidy_nw),
  if (!is.null(its_lbw))  bind_rows(its_lbw$tidy_ols,  its_lbw$tidy_nw),
  if (!is.null(its_nbw))  bind_rows(its_nbw$tidy_ols,  its_nbw$tidy_nw),
  if (!is.null(its_sens)) bind_rows(its_sens$tidy_ols, its_sens$tidy_nw)
)
write_csv(all_coefs, file.path(OUTPUT_DIR, "15b_sepsis_its_coefficients.csv"))
cat("\nSaved: 15b_sepsis_its_coefficients.csv\n")


# -----------------------------------------------------------------------------
# 4. PLOT (all-inborn)
# -----------------------------------------------------------------------------

if (!is.null(its_all)) {
  d <- its_all$d
  d_pre  <- d %>% filter(intervention == 0)
  d_post <- d %>% filter(intervention == 1)
  p <- ggplot(d, aes(x = adm_month_date)) +
    annotate("rect", xmin = min(d$adm_month_date), xmax = INT_START,
             ymin = -Inf, ymax = Inf, fill = FILL_HIST) +
    annotate("rect", xmin = INT_START, xmax = max(d$adm_month_date) + days(31),
             ymin = -Inf, ymax = Inf, fill = FILL_INT) +
    geom_vline(xintercept = INT_START, linetype = "dashed", colour = "grey45") +
    geom_line(data = d_post, aes(y = counterfactual),
              colour = COL_CF, linetype = "dotted", linewidth = 1) +
    geom_line(data = d_pre,  aes(y = fitted_its), colour = COL_HIST, linewidth = 1.2) +
    geom_line(data = d_post, aes(y = fitted_its), colour = COL_INT,  linewidth = 1.2) +
    geom_point(aes(y = pct_sepsis, colour = period, size = n_admissions), alpha = 0.8) +
    scale_colour_manual(name = "Period",
                        values = c("Historical" = COL_HIST, "Intervention" = COL_INT)) +
    scale_size_continuous(name = "N admissions", range = c(1.5, 5.5)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title = "ITS -- Neonatal sepsis incidence, SMCH inborn",
         subtitle = paste0("Sepsis definition: ", SEPSIS_DEF,
                           ". Jan 2023 - Mar 2026 (22 pre + 17 post).",
                           "\nEXPLORATORY: confounded by concurrent blood-culture/sepsis-diagnosis QI work."),
         x = NULL, y = "% admissions with sepsis",
         caption = "pct_sepsis ~ time + intervention + time_after, weighted by admissions.") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(colour = "grey40", size = 8.8),
          legend.position = "right",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(file.path(OUTPUT_DIR, "15c_sepsis_its_plot.png"),
         plot = p, width = 13, height = 6, dpi = 300)
  cat("Saved: 15c_sepsis_its_plot.png\n")
}


# -----------------------------------------------------------------------------
# 5. SUMMARY
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 64), "\n  SEPSIS ITS RESULTS SUMMARY\n", strrep("=", 64), "\n", sep = "")
print_summary <- function(res, label) {
  if (is.null(res)) { cat(sprintf("\n%s: not fitted.\n", label)); return(invisible()) }
  coefs <- if (res$lb$p.value < 0.05) res$tidy_nw else res$tidy_ols
  cat(sprintf("\n--- %s (%s) ---\n", label,
              if (res$lb$p.value < 0.05) "Newey-West" else "OLS"))
  for (trm in c("time", "intervention", "time_after")) {
    r <- coefs %>% filter(term == trm)
    if (nrow(r) == 1)
      cat(sprintf("  %-13s est=%.3f [%.3f, %.3f] p=%.4f%s\n",
                  trm, r$estimate, r$ci_lower, r$ci_upper, r$p_value,
                  if (r$p_value < 0.05) " *" else ""))
  }
}
print_summary(its_all,  "All inborn")
print_summary(its_lbw,  "LBW (<2500g)")
print_summary(its_nbw,  "NBW (>=2500g)")
print_summary(its_sens, "All inborn (sensitivity incl. 2022)")

cat("\n  REMINDER: interpret with caution -- concurrent sepsis-diagnosis QI work\n")
cat("  (Tinashe PhD) is a major confounder for any change in recorded sepsis.\n")
cat("\n=== Script 15 complete ===\n")
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
