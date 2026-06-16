# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 09: DCC / ESSC Uptake STATUS BREAKDOWN (yes / no / unknown)
# =============================================================================
# Purpose (Action 1, 2026-06-09 meeting):
#   Verify the DCC / S2S(ESSC) uptake coding and report, by month, the FULL
#   status breakdown -- recorded vs not, and of all eligible babies the
#   proportion coded yes / no / unknown.
#
#   WHY THIS SCRIPT EXISTS:
#   Script 08 reports "uptake" as Y / (Y + N), i.e. it EXCLUDES the Unknown (U)
#   records from the denominator. That produces the ~98-99% DCC figure shown in
#   the June slides. Rachel flagged this looked high versus Metabase (~78-90%).
#   The explanation is the denominator: roughly 1 in 6 eligible babies are coded
#   U (unknown). When the percentage is expressed as a share of ALL eligible
#   babies, DCC "yes" is ~81%, "unknown" ~18%, "no" ~1% -- much closer to
#   Metabase. This script makes that breakdown explicit so the team can choose
#   the reporting denominator deliberately.
#
#   For each month it reports, as a percentage of ALL eligible babies:
#     pct_recorded  = (Y + N + U) / n_eligible        (status present at all)
#     pct_yes       = Y / n_eligible
#     pct_no        = N / n_eligible
#     pct_unknown   = U / n_eligible
#     pct_not_recorded = blank / n_eligible
#   plus, for comparison with Script 08 / the slides:
#     pct_yes_of_known = Y / (Y + N)                  (the "98-99%" definition)
#
#   Tables produced for:
#     (A) All eligible live births          -- DCC and ESSC
#     (C) LBW eligible live births (<2500 g) -- DCC and ESSC
#
# POPULATION / DEFINITIONS (identical to Script 08 Part B/C):
#   Source     : maternal outcome script (Prisca's delivery log).
#   Facility   : SMCH (all rows are inborn by construction).
#   Window     : 12-month intervention period, Nov 2024 - Oct 2025.
#                (Change INT_END to extend; later months are sparse.)
#   Live birth : neotreeoutcome in {LB, ENND}.
#   Eligible   : apgar1 > 6 OR resus == "N" (not requiring resuscitation).
#   LBW        : bwtdis < 2500 g.
#   Status     : dcc / s2s coded Y = received, N = not received,
#                U = unknown, "" / NA = not recorded.
#
# OUTPUTS (all to 03-OUTPUTS/Q1_uptake_status_tables/):
#   09a_dcc_all_eligible.csv
#   09b_essc_all_eligible.csv
#   09c_dcc_lbw_eligible.csv
#   09d_essc_lbw_eligible.csv
#   09e_status_breakdown_plot.png        -- stacked yes/no/unknown by month
#   09_uptake_status_breakdown_log.txt
#
# DSH note: ASCII-only (no non-ASCII characters anywhere).
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

MATERNAL_DATA <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA",
            "zim_db_maternal_outcomes_20260501_cleaned.csv"),
  mustWork = FALSE
)

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS", "Q1_uptake_status_tables"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(OUTPUT_DIR, "09_uptake_status_breakdown_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

# 12-month intervention window (Q1). Extend by changing INT_END.
INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")

COL_YES <- "#2166ac"   # blue
COL_NO  <- "#d6604d"   # orange-red
COL_UNK <- "#bdbdbd"   # grey
COL_NR  <- "#f0f0f0"   # near-white (not recorded)

cat("=============================================================\n")
cat("  Script 09: DCC/ESSC Uptake STATUS BREAKDOWN (yes/no/unknown)\n")
cat("=============================================================\n")
cat(sprintf("  Window           : %s to %s (12-month)\n", INT_START, INT_END))
cat(sprintf("  Source           : %s\n", MATERNAL_DATA))
cat(sprintf("  Output directory : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")


# -----------------------------------------------------------------------------
# 1. LOAD AND PREPARE
# -----------------------------------------------------------------------------

cat("Loading maternal outcomes data...\n")
mo_raw <- read_csv(MATERNAL_DATA, show_col_types = FALSE, guess_max = 100000)
cat(sprintf("  Rows loaded: %d\n", nrow(mo_raw)))

mo <- mo_raw %>%
  mutate(
    adm_dt      = ymd_hms(dateadmission, quiet = TRUE),
    adm_date    = as.Date(adm_dt),
    month_label = format(floor_date(adm_dt, "month"), "%Y-%m"),

    apgar1_num  = suppressWarnings(as.numeric(apgar1)),
    bwtdis_num  = suppressWarnings(as.numeric(bwtdis)),

    live_birth  = neotreeoutcome %in% c("LB", "ENND"),
    eligible    = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus == "N"                         ~ TRUE,
      TRUE                                 ~ FALSE
    ),
    lbw         = !is.na(bwtdis_num) & bwtdis_num < 2500,

    dcc_status  = case_when(
      dcc == "Y"                      ~ "yes",
      dcc == "N"                      ~ "no",
      dcc == "U"                      ~ "unknown",
      is.na(dcc) | trimws(dcc) == "" ~ "not_recorded",
      TRUE                            ~ "other"
    ),
    s2s_status  = case_when(
      s2s == "Y"                      ~ "yes",
      s2s == "N"                      ~ "no",
      s2s == "U"                      ~ "unknown",
      is.na(s2s) | trimws(s2s) == "" ~ "not_recorded",
      TRUE                            ~ "other"
    )
  )

mo_int <- mo %>%
  filter(facility == "SMCH",
         adm_date >= INT_START, adm_date <= INT_END,
         live_birth == TRUE)

cat(sprintf("\nLive births (SMCH, %s to %s): %d\n",
            INT_START, INT_END, nrow(mo_int)))
cat(sprintf("  Eligible (apgar1>6 or N-resus): %d\n",
            sum(mo_int$eligible, na.rm = TRUE)))
cat(sprintf("  LBW eligible (<2500 g)        : %d\n",
            sum(mo_int$eligible & mo_int$lbw, na.rm = TRUE)))


# -----------------------------------------------------------------------------
# 2. STATUS-BREAKDOWN BUILDER
# -----------------------------------------------------------------------------
# Denominator for all pct_* columns (except pct_yes_of_known) is n_eligible,
# i.e. ALL eligible babies in the month, so columns sum to ~100%.

MONTHS <- format(seq(INT_START, INT_END, by = "month"), "%Y-%m")

build_status <- function(data, status_col) {
  out <- lapply(MONTHS, function(m) {
    d  <- data %>% filter(month_label == m)
    n  <- nrow(d)
    s  <- d[[status_col]]
    ny <- sum(s == "yes",          na.rm = TRUE)
    nn <- sum(s == "no",           na.rm = TRUE)
    nu <- sum(s == "unknown",      na.rm = TRUE)
    nr <- sum(s == "not_recorded", na.rm = TRUE)
    nk <- ny + nn + nu                       # recorded (Y+N+U)
    pc <- function(x) if (n > 0) round(100 * x / n, 1) else NA_real_
    tibble(
      month            = m,
      n_eligible       = n,
      n_recorded       = nk,  pct_recorded     = pc(nk),
      n_yes            = ny,  pct_yes          = pc(ny),
      n_no             = nn,  pct_no           = pc(nn),
      n_unknown        = nu,  pct_unknown      = pc(nu),
      n_not_recorded   = nr,  pct_not_recorded = pc(nr),
      pct_yes_of_known = if ((ny + nn) > 0) round(100 * ny / (ny + nn), 1)
                         else NA_real_
    )
  })
  tab <- bind_rows(out)
  # TOTAL row across the window
  n  <- sum(tab$n_eligible)
  pc <- function(x) if (n > 0) round(100 * x / n, 1) else NA_real_
  ty <- sum(tab$n_yes); tn <- sum(tab$n_no)
  tot <- tibble(
    month            = "TOTAL (12mo)",
    n_eligible       = n,
    n_recorded       = sum(tab$n_recorded),     pct_recorded     = pc(sum(tab$n_recorded)),
    n_yes            = ty,                       pct_yes          = pc(ty),
    n_no             = tn,                       pct_no           = pc(tn),
    n_unknown        = sum(tab$n_unknown),       pct_unknown      = pc(sum(tab$n_unknown)),
    n_not_recorded   = sum(tab$n_not_recorded),  pct_not_recorded = pc(sum(tab$n_not_recorded)),
    pct_yes_of_known = if ((ty + tn) > 0) round(100 * ty / (ty + tn), 1) else NA_real_
  )
  bind_rows(tab, tot)
}

elig     <- mo_int %>% filter(eligible == TRUE)
elig_lbw <- mo_int %>% filter(eligible == TRUE, lbw == TRUE)

tab_dcc_all  <- build_status(elig,     "dcc_status")
tab_essc_all <- build_status(elig,     "s2s_status")
tab_dcc_lbw  <- build_status(elig_lbw, "dcc_status")
tab_essc_lbw <- build_status(elig_lbw, "s2s_status")

write_csv(tab_dcc_all,  file.path(OUTPUT_DIR, "09a_dcc_all_eligible.csv"))
write_csv(tab_essc_all, file.path(OUTPUT_DIR, "09b_essc_all_eligible.csv"))
write_csv(tab_dcc_lbw,  file.path(OUTPUT_DIR, "09c_dcc_lbw_eligible.csv"))
write_csv(tab_essc_lbw, file.path(OUTPUT_DIR, "09d_essc_lbw_eligible.csv"))

cat("\n--- DCC, all eligible babies ---\n");  print(as.data.frame(tab_dcc_all),  row.names = FALSE)
cat("\n--- ESSC, all eligible babies ---\n"); print(as.data.frame(tab_essc_all), row.names = FALSE)
cat("\n--- DCC, LBW eligible babies ---\n");  print(as.data.frame(tab_dcc_lbw),  row.names = FALSE)
cat("\n--- ESSC, LBW eligible babies ---\n"); print(as.data.frame(tab_essc_lbw), row.names = FALSE)

cat("\nSaved: 09a_dcc_all_eligible.csv, 09b_essc_all_eligible.csv,\n")
cat("       09c_dcc_lbw_eligible.csv, 09d_essc_lbw_eligible.csv\n")


# -----------------------------------------------------------------------------
# 3. KEY COMPARISON: two denominators side by side
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("  DENOMINATOR COMPARISON (overall, 12-month)\n")
cat("=============================================================\n")
cmp <- function(tab, lab) {
  t <- tab %>% filter(month == "TOTAL (12mo)")
  cat(sprintf("  %-18s yes/known = %.1f%%  |  yes/all-eligible = %.1f%%  (unknown %.1f%%)\n",
              lab, t$pct_yes_of_known, t$pct_yes, t$pct_unknown))
}
cmp(tab_dcc_all,  "DCC all")
cmp(tab_essc_all, "ESSC all")
cmp(tab_dcc_lbw,  "DCC LBW")
cmp(tab_essc_lbw, "ESSC LBW")
cat("\n  'yes/known' is the Script 08 / June-slide definition (excludes Unknown).\n")
cat("  'yes/all-eligible' counts Unknown in the denominator (recommended for\n")
cat("  population coverage; closer to Metabase).\n")


# -----------------------------------------------------------------------------
# 4. PLOT: stacked status by month (4 panels)
# -----------------------------------------------------------------------------

long_one <- function(tab, grp, intv) {
  # Use the raw COUNTS (not pre-rounded percentages). With geom_col(position =
  # "fill") ggplot normalises the counts per month to sum to exactly 100%, so
  # no segment is ever dropped by rounding (the earlier "blue Yes missing" bug).
  tab %>%
    filter(month != "TOTAL (12mo)") %>%
    select(month, n_yes, n_no, n_unknown, n_not_recorded) %>%
    pivot_longer(-month, names_to = "status", values_to = "n") %>%
    mutate(
      status = recode(status,
                      n_yes = "Yes", n_no = "No",
                      n_unknown = "Unknown", n_not_recorded = "Not recorded"),
      group = grp, intervention = intv, adm_month = ym(month)
    )
}

plot_df <- bind_rows(
  long_one(tab_dcc_all,  "All eligible", "DCC"),
  long_one(tab_essc_all, "All eligible", "ESSC"),
  long_one(tab_dcc_lbw,  "LBW eligible", "DCC"),
  long_one(tab_essc_lbw, "LBW eligible", "ESSC")
) %>%
  mutate(status = factor(status, levels = c("Yes", "No", "Unknown", "Not recorded")),
         panel  = paste0(intervention, " -- ", group))

p <- ggplot(plot_df, aes(x = adm_month, y = pct, fill = status)) +
  geom_col(width = 25) +
  facet_wrap(~ panel, ncol = 2) +
  scale_fill_manual(
    name   = "Recorded status",
    values = c("Yes" = COL_YES, "No" = COL_NO,
               "Unknown" = COL_UNK, "Not recorded" = COL_NR)
  ) +
  scale_y_continuous(name = "% of eligible babies",
                     breaks = seq(0, 100, 20)) +
  # use coord_cartesian (not scale limits) so rounding sums of 100.1% do not
  # drop the top stacked segment as "out of range" (ggplot oob clipping)
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_date(name = NULL, date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = "DCC / ESSC status by month -- share of all eligible babies",
    subtitle = paste0("Maternal outcome script, SMCH eligible live births, ",
                      INT_START, " to ", INT_END,
                      "\nUnknown is shown explicitly rather than dropped from the denominator"),
    caption  = "Denominator = all eligible babies each month. Bars sum to ~100%."
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
        plot.title       = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "09e_status_breakdown_plot.png"),
       plot = p, width = 11, height = 7, dpi = 150)
cat("\nSaved: 09e_status_breakdown_plot.png\n")

cat("\n=== Script 09 complete ===\n")
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
