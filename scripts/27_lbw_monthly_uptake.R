# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 27: DCC / ESSC UPTAKE BY MONTH, LOW BIRTH WEIGHT (<2500 g) ONLY
#            [Table 1a / 1b]
# =============================================================================
# Purpose:
#   Monthly uptake of DCC (Table 1a) and ESSC (Table 1b) among eligible
#   low-birthweight infants, in the same style as the Caesarean table
#   (Script 24).
#
# Definitions & methodology:
#   Data source : zim_db_maternal_outcomes_20260501_cleaned.csv.
#   Population  : maternal outcome script (delivery log), facility SMCH,
#                 live births (neotreeoutcome LB or ENND), birthweight at
#                 discharge (bwtdis) < 2500 g, range-filtered to 300-7000 g.
#   Eligibility : two definitions are produced --
#     (a) "table1_def" -- resus == "N" (DEFAULT, matches the Table 1 family)
#     (b) "table2_def" -- apgar1 > 6 OR resus == "N" (matches Table 2)
#     Both are written out so the intended definition can be confirmed
#     before use -- they will not give identical Ns.
#   Status      : dcc / s2s coded Y = received, N = not received,
#                 U = unknown, "" / NA = not recorded.
#   Percentages (per month, LBW eligible babies):
#     pct_yes          = Y / n_eligible_lbw          <- coverage figure
#     pct_yes_of_known = Y / (Y + N)                 <- "of documented" figure
#     pct_unknown      = (U + not recorded) / n_eligible_lbw
#   Window: 1 Nov 2024 - 31 Oct 2025 only (the key study period).
#
# Outputs (to ./outputs/maternal_august/):
#   27a_lbw_monthly_uptake_table1def.csv   -- default, resus == "N" eligibility
#   27b_lbw_monthly_uptake_table2def.csv   -- cross-check, apgar1>6 OR resus=N
#   27_lbw_monthly_uptake_log.txt
#   Golden_Hour_Table1a_1b_LBW_monthly_uptake.docx -- Table 1a/1b as Word tables
#
# Requires: tidyverse, lubridate, officer, flextable (docx output).
#
# DSH note: script is ASCII-only.
# =============================================================================

library(tidyverse)
library(lubridate)
library(officer)
library(flextable)

args <- commandArgs(trailingOnly = FALSE)
sf <- args[grep("--file=", args)]
SCRIPT_DIR <- if (length(sf) > 0)
  dirname(normalizePath(sub("--file=", "", sf), mustWork = FALSE)) else getwd()

MATERNAL <- normalizePath(file.path(SCRIPT_DIR, "..", "00-DATA",
  "zim_db_maternal_outcomes_20260501_cleaned.csv"), mustWork = FALSE)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs", "maternal_august")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(OUTPUT_DIR, "27_lbw_monthly_uptake_log.txt")
log_con <- file(LOG_FILE, open = "wt"); sink(log_con, split = TRUE, type = "output")

INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2025-10-31")
LBW_CUTOFF_G <- 2500

# Read every column as character: readr's type guesser otherwise mangles
# dateadmission for records outside its initial sample window.
raw <- read_csv(MATERNAL, show_col_types = FALSE,
                col_types = cols(.default = col_character()))

df <- raw %>%
  mutate(
    adm_date  = as.Date(dateadmission),
    month_label = format(adm_date, "%Y-%m"),
    live_birth = neotreeoutcome %in% c("LB", "ENND"),
    apgar1_num = suppressWarnings(as.numeric(apgar1)),

    bwtdis_num = suppressWarnings(as.numeric(bwtdis)),
    bwtdis_num = if_else(!is.na(bwtdis_num) & (bwtdis_num < 300 | bwtdis_num > 7000),
                          NA_real_, bwtdis_num),
    lbw = !is.na(bwtdis_num) & bwtdis_num < LBW_CUTOFF_G,

    # (a) Table 1 definition (Script 23): did NOT require immediate resuscitation.
    eligible_table1_def = resus == "N",
    # (b) Table 2 proxy definition (Script 24), for cross-checking.
    eligible_table2_def = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus == "N"                        ~ TRUE,
      TRUE                                ~ FALSE),

    dcc_Y = dcc == "Y", dcc_N = dcc == "N", dcc_U = dcc == "U",
    dcc_NR = is.na(dcc) | trimws(dcc) == "",
    s2s_Y = s2s == "Y", s2s_N = s2s == "N", s2s_U = s2s == "U",
    s2s_NR = is.na(s2s) | trimws(s2s) == "")

sum_na <- function(x) sum(x, na.rm = TRUE)

monthly_lbw <- function(elig_var) {
  d <- df %>% filter(facility == "SMCH", live_birth == TRUE, lbw == TRUE,
                     .data[[elig_var]] == TRUE,
                     adm_date >= INT_START, adm_date <= INT_END)
  body <- d %>% group_by(month_label) %>% summarise(
    n_eligible_lbw = n(),
    dcc_Y = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR),
    .groups = "drop")
  tot <- d %>% summarise(month_label = "TOTAL",
    n_eligible_lbw = n(),
    dcc_Y = sum_na(dcc_Y), dcc_N = sum_na(dcc_N),
    dcc_U = sum_na(dcc_U), dcc_not_recorded = sum_na(dcc_NR),
    essc_Y = sum_na(s2s_Y), essc_N = sum_na(s2s_N),
    essc_U = sum_na(s2s_U), essc_not_recorded = sum_na(s2s_NR))
  bind_rows(arrange(body, month_label), tot) %>%
    mutate(
      dcc_pct_yes          = round(100 * dcc_Y / n_eligible_lbw, 1),
      dcc_pct_yes_of_known = if_else(dcc_Y + dcc_N > 0,
                                     round(100 * dcc_Y / (dcc_Y + dcc_N), 1), NA_real_),
      dcc_pct_unknown      = round(100 * (dcc_U + dcc_not_recorded) / n_eligible_lbw, 1),
      essc_pct_yes          = round(100 * essc_Y / n_eligible_lbw, 1),
      essc_pct_yes_of_known = if_else(essc_Y + essc_N > 0,
                                      round(100 * essc_Y / (essc_Y + essc_N), 1), NA_real_),
      essc_pct_unknown      = round(100 * (essc_U + essc_not_recorded) / n_eligible_lbw, 1)) %>%
    relocate(month_label, n_eligible_lbw,
             dcc_Y, dcc_pct_yes, dcc_pct_yes_of_known, dcc_N, dcc_U, dcc_not_recorded, dcc_pct_unknown,
             essc_Y, essc_pct_yes, essc_pct_yes_of_known, essc_N, essc_U, essc_not_recorded, essc_pct_unknown)
}

cat("=============================================================\n")
cat("  Script 27: LBW (<2500g) DCC/ESSC uptake by month (Table 1a/1b)\n")
cat(sprintf("  Window : %s to %s\n", INT_START, INT_END))
cat("=============================================================\n\n")

tb_table1def <- monthly_lbw("eligible_table1_def")
cat("\n=========== (a) TABLE 1 eligibility def (resus == N) ===========\n")
print(as.data.frame(tb_table1def), row.names = FALSE)
write_csv(tb_table1def, file.path(OUTPUT_DIR, "27a_lbw_monthly_uptake_table1def.csv"))

tb_table2def <- monthly_lbw("eligible_table2_def")
cat("\n=========== (b) TABLE 2 eligibility def (apgar1>6 OR resus=N), for cross-check ===========\n")
print(as.data.frame(tb_table2def), row.names = FALSE)
write_csv(tb_table2def, file.path(OUTPUT_DIR, "27b_lbw_monthly_uptake_table2def.csv"))

# ---- Context for the footnotes ---------------------------------------------
cat("\n-------------------------------------------------------------\n")
cat("  CONTEXT / SENSITIVITY\n")
cat("-------------------------------------------------------------\n")
ctx <- df %>% filter(facility == "SMCH", live_birth == TRUE,
                     adm_date >= INT_START, adm_date <= INT_END)
cat(sprintf("12-month window, SMCH live births                     : %d\n", nrow(ctx)))
cat(sprintf("  of which LBW (bwtdis < %d g)                       : %d\n",
            LBW_CUTOFF_G, sum_na(ctx$lbw)))
cat(sprintf("  LBW AND eligible (table1_def: resus = N)             : %d\n",
            sum(ctx$lbw & ctx$eligible_table1_def, na.rm = TRUE)))
cat(sprintf("  LBW AND eligible (table2_def: apgar1>6 OR resus=N)   : %d\n",
            sum(ctx$lbw & ctx$eligible_table2_def, na.rm = TRUE)))
cat(sprintf("  Missing/out-of-range bwtdis in this window           : %d\n",
            sum(is.na(ctx$bwtdis_num) | is.na(as.numeric(ctx$bwtdis)) == FALSE &
                (as.numeric(ctx$bwtdis) < 300 | as.numeric(ctx$bwtdis) > 7000), na.rm = TRUE)))

# =============================================================================
# DOCX OUTPUT -- Table 1a (DCC) and Table 1b (ESSC), styled to match Table 2a/2b
# =============================================================================
# Column layout matches the existing Table 2a/2b exactly:
#   Month | Eligible [low birth weight] infants, n | Received, n (% of eligible)
#   | Not received, n | Unknown / not recorded, n (%) | % of documented (Y/N)
# Built from the DEFAULT eligibility definition (table1_def). To use the
# table2_def cross-check instead, swap tb_table1def for tb_table2def in the
# two format_intervention() calls below.

fmt_pct  <- function(x) if_else(is.na(x), "-", sprintf("%.1f%%", x))
fmt_cell <- function(n, pct) sprintf("%d (%s)", n, fmt_pct(pct))

format_intervention <- function(tb, prefix) {
  tb %>% transmute(
    Month = month_label,
    `Eligible low birth weight infants, n` = n_eligible_lbw,
    `Received, n (% of eligible)` = fmt_cell(.data[[paste0(prefix, "_Y")]],
                                              .data[[paste0(prefix, "_pct_yes")]]),
    `Not received, n` = .data[[paste0(prefix, "_N")]],
    `Unknown / not recorded, n (%)` = fmt_cell(
      .data[[paste0(prefix, "_U")]] + .data[[paste0(prefix, "_not_recorded")]],
      .data[[paste0(prefix, "_pct_unknown")]]),
    `% of documented (Y/N)` = fmt_pct(.data[[paste0(prefix, "_pct_yes_of_known")]]))
}

table1a_dcc  <- format_intervention(tb_table1def, "dcc")
table1b_essc <- format_intervention(tb_table1def, "essc")

style_table <- function(df_table) {
  ft <- flextable(df_table) %>%
    theme_box() %>%
    bold(part = "header") %>%
    bold(i = ~ Month == "TOTAL") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    autofit()
  ft
}

doc <- read_docx() %>%
  body_add_par("Table 1a. Delayed cord clamping among eligible low birth weight infants (<2500g), by month (1 November 2024 - 31 October 2025).",
               style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(style_table(table1a_dcc)) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table 1b. Early skin to skin care among eligible low birth weight infants (<2500g), by month (1 November 2024 - 31 October 2025).",
               style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(style_table(table1b_essc)) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Eligible = did not require immediate resuscitation (resus = No), per the Table 1 definition. Low birth weight = birthweight at discharge (bwtdis) < 2500 g, values <300 g or >7000 g treated as implausible/missing. If this should instead match Table 2's eligibility definition (1-minute Apgar >6 or resuscitation not required), re-run using tb_table2def -- see script header.",
               style = "Normal")

DOCX_PATH <- file.path(OUTPUT_DIR, "Golden_Hour_Table1a_1b_LBW_monthly_uptake.docx")
print(doc, target = DOCX_PATH)
cat(sprintf("\nDocx written to: %s\n", DOCX_PATH))

cat("\n=== Script 27 complete ===\n")
sink(type = "output"); close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
message(sprintf("Docx saved to: %s", DOCX_PATH))
