# =============================================================================
# Golden Hour Analysis -- SMCH Zimbabwe
# Script 08: DCC and S2S Uptake -- NeoTree vs Maternal Outcome Script
# =============================================================================
# Purpose:
#   Compare DCC and skin-to-skin (S2S/ESSC) uptake across two parallel data
#   streams over the analysis period (November 2024 - March 2026):
#     Part A -- NeoTree data (inborn neonatal admissions)
#     Part B -- Maternal outcome script (primary source)
#     Part C -- LBW subgroup (bwtdis < 2500 g), from the maternal outcome
#               script, eligible babies only
#
# Definitions & methodology:
#   Data sources: ZIM_db_master_joined_to_20260525.csv (NeoTree master) and
#   zim_db_maternal_outcomes_20260501_cleaned.csv (maternal outcomes).
#
#   Part A: monthly DCC and S2S uptake from the 'delivinter' multi-select
#   variable. Denominator: known delivinter (recorded and not coded Unknown).
#
#   Part B: monthly DCC and S2S uptake from the maternal delivery-log data.
#   dcc / s2s variables: Y = received, N = not received, U = unknown,
#   "" = not recorded. Denominator for uptake rate: known (Y + N only).
#   Secondary metric: unknown rate (U / all recorded) per month.
#
#   Eligibility (not requiring immediate resuscitation) differs by part,
#   since Part A and Part B/C read different source data:
#     Part A (NeoTree data): apgar1 > 6 OR resus %in% c("NONE","Norm").
#       Babies where eligibility cannot be determined are flagged separately
#       but retained in totals.
#     Part B/C (maternal outcome script): resus == "N" ("Did baby require
#       immediate resuscitation? = No") -- the confirmed primary/reported
#       definition, in step with Scripts 23/24/25. A broader apgar1>6 OR
#       resus=="N" proxy is also computed and written to a separate,
#       clearly-labelled cross-check output only.
#
#   Note on "inborn" in the maternal outcome script: it is completed
#   directly from the SMCH delivery log, so all entries represent babies
#   born at SMCH. No inborn/outborn filter variable is present. The
#   'externalsource' variable (ER = ~3% of rows) records that the MOTHER
#   was referred from an external facility -- the baby is still born at
#   SMCH and is included in all analyses.
#
# Outputs (to 03-OUTPUTS/):
#   08a_neotree_monthly_uptake.csv
#   08b_maternal_monthly_uptake.csv       (primary, resus=="N")
#   08c_lbw_monthly_uptake.csv
#   08d_neotree_monthly_uptake.png
#   08e_maternal_monthly_uptake.png
#   08f_maternal_unknown_rate.png
#   08g_lbw_monthly_uptake.png
#   08h_maternal_monthly_uptake_broadproxy_crosscheck.csv (cross-check only)
#   08_maternal_script_uptake_log.txt
#
# DSH note: script is ASCII-only.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(scales)
library(ggplot2)

# -- Locate script directory ---------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("--file=", args)]
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("--file=", "", script_flag),
                                      mustWork = FALSE))
} else {
  SCRIPT_DIR <- getwd()
}

# -- Path configuration -------------------------------------------------------
NEOTREE_DATA <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA",
            "ZIM_db_master_joined_to_20260525.csv"),
  mustWork = FALSE
)

MATERNAL_DATA <- normalizePath(
  file.path(SCRIPT_DIR, "..", "00-DATA",
            "zim_db_maternal_outcomes_20260501_cleaned.csv"),
  mustWork = FALSE
)

OUTPUT_DIR <- normalizePath(
  file.path(SCRIPT_DIR, "..", "03-OUTPUTS"),
  mustWork = FALSE
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -- Log file -----------------------------------------------------------------
LOG_FILE <- file.path(OUTPUT_DIR, "08_maternal_script_uptake_log.txt")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE, type = "output")

# -- Study period constants ---------------------------------------------------
INT_START <- as.Date("2024-11-01")
INT_END   <- as.Date("2026-03-31")

# -- Plot theme ---------------------------------------------------------------
# Consistent minimal theme used across all script 08 figures.
theme_gh <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      legend.position   = "bottom",
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 10, colour = "grey40"),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 9)
    )
}

# Colour palette (accessible, no red/green confusion)
COL_DCC <- "#2166ac"   # blue
COL_S2S <- "#d6604d"   # orange-red
COL_UNK <- "#999999"   # grey for unknown rate

cat("=============================================================\n")
cat("  Script 08: DCC/S2S Uptake -- NeoTree vs Maternal Script\n")
cat("=============================================================\n")
cat(sprintf("  Intervention period : %s to %s\n", INT_START, INT_END))
cat(sprintf("  NeoTree data        : %s\n", NEOTREE_DATA))
cat(sprintf("  Maternal outcomes   : %s\n", MATERNAL_DATA))
cat(sprintf("  Output directory    : %s\n", OUTPUT_DIR))
cat("=============================================================\n\n")


# =============================================================================
# PART A -- NEOTREE UPTAKE (inborn neonatal admissions, SMCH)
# =============================================================================

cat("=============================================================\n")
cat("  PART A: NeoTree DCC/S2S Uptake\n")
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# A1. Load and prepare NeoTree data
# -----------------------------------------------------------------------------

cat("Loading NeoTree master data...\n")
nt_raw <- read_csv(NEOTREE_DATA,
                   show_col_types = FALSE,
                   guess_max      = 100000)

parse_probs <- problems(nt_raw)
if (nrow(parse_probs) > 0) {
  cat(sprintf("  WARNING: %d parse problem(s) in NeoTree data.\n",
              nrow(parse_probs)))
  print(head(parse_probs, 10))
} else {
  cat("  NeoTree data: no parse problems.\n")
}
cat(sprintf("  Rows loaded: %d\n", nrow(nt_raw)))

# -- Derived variables --------------------------------------------------------
nt <- nt_raw %>%
  mutate(
    adm_dt   = ymd_hms(datetimeadmission, quiet = TRUE),
    adm_date = as.Date(adm_dt),
    adm_month = floor_date(adm_dt, "month"),
    month_label = format(adm_month, "%Y-%m"),

    inborn = case_when(
      inorout %in% c("Yes", "true", "True", "In")    ~ TRUE,
      inorout %in% c("No",  "false", "False", "Out") ~ FALSE,
      TRUE                                            ~ NA
    ),
    is_readmission = readmission == "Y",

    # -- delivinter parsing (multi-select; "U" = unknown) --------------------
    dcc_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                    ~ NA,
      grepl("DCC", delivinter, fixed = TRUE)          ~ TRUE,
      TRUE                                            ~ FALSE
    ),
    s2s_received = case_when(
      is.na(delivinter) | trimws(delivinter) == "" |
        trimws(delivinter) == "U"                    ~ NA,
      grepl("S2S", delivinter, fixed = TRUE)          ~ TRUE,
      TRUE                                            ~ FALSE
    ),
    # TRUE = delivinter recorded and not Unknown
    delivinter_known = (
      !is.na(delivinter) &
      trimws(delivinter) != "" &
      trimws(delivinter) != "U"
    ),
    delivinter_missing = is.na(delivinter) | trimws(delivinter) == "",
    delivinter_unknown = !is.na(delivinter) & trimws(delivinter) == "U",

    # -- Eligibility proxy: apgar1 > 6 OR no resuscitation needed ------------
    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    resus_none = resus %in% c("NONE", "Norm"),
    eligible = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus_none == TRUE                   ~ TRUE,
      TRUE                                 ~ FALSE   # conservative: unknown -> ineligible
    )
  )

# -- Intervention period, SMCH inborn, non-readmission -----------------------
int_pop <- nt %>%
  filter(
    facility      == "SMCH",
    inborn        == TRUE,
    is_readmission == FALSE | is.na(is_readmission),
    adm_date      >= INT_START,
    adm_date      <= INT_END
  )

cat(sprintf("\nNeoTree intervention population (SMCH inborn, non-readmission): %d\n",
            nrow(int_pop)))
cat(sprintf("  With delivinter known (excl. U/missing): %d (%.1f%%)\n",
            sum(int_pop$delivinter_known, na.rm = TRUE),
            100 * mean(int_pop$delivinter_known, na.rm = TRUE)))
cat(sprintf("  With delivinter missing                : %d (%.1f%%)\n",
            sum(int_pop$delivinter_missing, na.rm = TRUE),
            100 * mean(int_pop$delivinter_missing, na.rm = TRUE)))
cat(sprintf("  With delivinter Unknown (U)            : %d (%.1f%%)\n",
            sum(int_pop$delivinter_unknown, na.rm = TRUE),
            100 * mean(int_pop$delivinter_unknown, na.rm = TRUE)))
cat(sprintf("  Eligible proxy (apgar1>6 or no resus)  : %d (%.1f%%)\n",
            sum(int_pop$eligible, na.rm = TRUE),
            100 * mean(int_pop$eligible, na.rm = TRUE)))

# -- Overall NeoTree uptake (all admitted, not restricted to eligible) --------
# Denominator: delivinter_known. This mirrors Script 01 Section 5.
n_nt_known <- sum(int_pop$delivinter_known, na.rm = TRUE)
cat(sprintf("\nOverall NeoTree uptake (denominator = delivinter_known, n=%d):\n",
            n_nt_known))
nt_known <- int_pop %>% filter(delivinter_known == TRUE)
cat(sprintf("  DCC received    : %d / %d (%.1f%%)\n",
            sum(nt_known$dcc_received, na.rm = TRUE), n_nt_known,
            100 * mean(nt_known$dcc_received, na.rm = TRUE)))
cat(sprintf("  S2S received    : %d / %d (%.1f%%)\n",
            sum(nt_known$s2s_received, na.rm = TRUE), n_nt_known,
            100 * mean(nt_known$s2s_received, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# A2. NeoTree monthly uptake table
# -----------------------------------------------------------------------------

cat("\n--- NeoTree monthly DCC and S2S uptake (all admitted, delivinter_known denom) ---\n")

nt_monthly <- int_pop %>%
  group_by(month_label) %>%
  summarise(
    n_admitted        = n(),
    n_dcc_recorded    = sum(!is.na(dcc_received)),
    n_s2s_recorded    = sum(!is.na(s2s_received)),
    n_dcc_unknown     = sum(delivinter_unknown, na.rm = TRUE),
    n_dv_missing      = sum(delivinter_missing, na.rm = TRUE),
    n_delivinter_known = sum(delivinter_known, na.rm = TRUE),
    n_dcc_received    = sum(dcc_received == TRUE,  na.rm = TRUE),
    n_s2s_received    = sum(s2s_received == TRUE,  na.rm = TRUE),
    n_dcc_not_recv    = sum(dcc_received == FALSE, na.rm = TRUE),
    n_s2s_not_recv    = sum(s2s_received == FALSE, na.rm = TRUE),
    pct_dcc           = if_else(n_delivinter_known > 0,
                                round(100 * n_dcc_received / n_delivinter_known, 1),
                                NA_real_),
    pct_s2s           = if_else(n_delivinter_known > 0,
                                round(100 * n_s2s_received / n_delivinter_known, 1),
                                NA_real_),
    pct_dv_missing    = round(100 * n_dv_missing / n_admitted, 1),
    pct_dcc_unknown   = round(100 * n_dcc_unknown / n_admitted, 1),
    .groups = "drop"
  ) %>%
  arrange(month_label)

print(nt_monthly, n = Inf)
write_csv(nt_monthly, file.path(OUTPUT_DIR, "08a_neotree_monthly_uptake.csv"))
cat("\nSaved: 08a_neotree_monthly_uptake.csv\n")

# -----------------------------------------------------------------------------
# A3. NeoTree monthly uptake chart
# -----------------------------------------------------------------------------

# Reshape for plotting (DCC and S2S side by side)
nt_monthly_long <- nt_monthly %>%
  select(month_label, pct_dcc, pct_s2s) %>%
  pivot_longer(
    cols      = c(pct_dcc, pct_s2s),
    names_to  = "intervention",
    values_to = "pct"
  ) %>%
  mutate(
    intervention = recode(intervention,
                          pct_dcc = "DCC",
                          pct_s2s = "S2S (ESSC)"),
    adm_month = ym(month_label),
    col       = if_else(intervention == "DCC", COL_DCC, COL_S2S)
  )

nt_col_map <- c("DCC" = COL_DCC, "S2S (ESSC)" = COL_S2S)

p_nt <- ggplot(nt_monthly_long,
               aes(x = adm_month, y = pct, colour = intervention,
                   group = intervention)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  # Annotate n_delivinter_known per month
  geom_text(
    data = nt_monthly %>% mutate(adm_month = ym(month_label)),
    aes(x = adm_month, y = 2, label = paste0("n=", n_delivinter_known)),
    inherit.aes = FALSE,
    size = 2.8, colour = "grey40", vjust = 0
  ) +
  scale_y_continuous(
    name   = "Uptake rate (%, of known delivinter)",
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_x_date(
    name        = NULL,
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_colour_manual(
    name   = "Intervention",
    values = nt_col_map
  ) +
  labs(
    title    = "NeoTree DCC and S2S Uptake by Month",
    subtitle = paste0("SMCH inborn neonatal admissions, intervention period ",
                      INT_START, " to ", INT_END,
                      "\nDenominator: babies with known delivinter (excl. Unknown/missing)"),
    caption  = "Note: n = number with known delivinter each month; pct = n received / n known"
  ) +
  theme_gh()

ggsave(file.path(OUTPUT_DIR, "08d_neotree_monthly_uptake.png"),
       plot   = p_nt,
       width  = 10,
       height = 5.5,
       dpi    = 150)
cat("Saved: 08d_neotree_monthly_uptake.png\n")


# =============================================================================
# PART B -- MATERNAL OUTCOME SCRIPT UPTAKE (primary source)
# =============================================================================

cat("\n\n=============================================================\n")
cat("  PART B: Maternal Outcome Script DCC/S2S Uptake\n")
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# B1. Load and prepare maternal outcomes data
# -----------------------------------------------------------------------------

cat("Loading maternal outcomes data...\n")
mo_raw <- read_csv(MATERNAL_DATA,
                   show_col_types = FALSE,
                   guess_max      = 100000)

parse_probs_mo <- problems(mo_raw)
if (nrow(parse_probs_mo) > 0) {
  cat(sprintf("  WARNING: %d parse problem(s) in maternal outcomes data.\n",
              nrow(parse_probs_mo)))
  print(head(parse_probs_mo, 10))
} else {
  cat("  Maternal outcomes data: no parse problems.\n")
}
cat(sprintf("  Rows loaded: %d\n", nrow(mo_raw)))

# -- Derived variables --------------------------------------------------------
mo <- mo_raw %>%
  mutate(
    adm_dt    = ymd_hms(dateadmission, quiet = TRUE),
    adm_date  = as.Date(adm_dt),
    adm_month = floor_date(adm_dt, "month"),
    month_label = format(adm_month, "%Y-%m"),

    apgar1_num = suppressWarnings(as.numeric(apgar1)),
    bwtdis_num = suppressWarnings(as.numeric(bwtdis)),

    # -- Live birth flag -------------------------------------------------------
    # Restrict to live births (LB) and early neonatal deaths (ENND, who were
    # alive at birth). Exclude stillbirths (STBM, STBF) and miscarriages (MISC)
    # since DCC/S2S are not applicable to these outcomes.
    live_birth = neotreeoutcome %in% c("LB", "ENND"),

    # -- Eligibility: "Did baby require immediate resuscitation? = No" -------
    # Primary/reported definition, in step with Scripts 23/24/25.
    eligible = resus == "N",
    # Broader proxy, cross-check only -- NOT used for the primary tables below.
    eligible_broadproxy = case_when(
      !is.na(apgar1_num) & apgar1_num > 6 ~ TRUE,
      resus == "N"                         ~ TRUE,
      TRUE                                 ~ FALSE
    ),

    # -- DCC coding -----------------------------------------------------------
    # Y = DCC received, N = not received, U = unknown, "" = not recorded
    dcc_code = case_when(
      dcc == "Y"                        ~ "received",
      dcc == "N"                        ~ "not_received",
      dcc == "U"                        ~ "unknown",
      is.na(dcc) | trimws(dcc) == ""   ~ "not_recorded",
      TRUE                              ~ "other"
    ),
    dcc_received    = dcc == "Y",
    dcc_not_recv    = dcc == "N",
    dcc_unknown     = dcc == "U",
    dcc_not_recorded = is.na(dcc) | trimws(dcc) == "",
    dcc_known       = dcc %in% c("Y", "N"),   # known outcome (not U, not missing)

    # -- S2S coding -----------------------------------------------------------
    s2s_code = case_when(
      s2s == "Y"                        ~ "received",
      s2s == "N"                        ~ "not_received",
      s2s == "U"                        ~ "unknown",
      is.na(s2s) | trimws(s2s) == ""   ~ "not_recorded",
      TRUE                              ~ "other"
    ),
    s2s_received    = s2s == "Y",
    s2s_not_recv    = s2s == "N",
    s2s_unknown     = s2s == "U",
    s2s_not_recorded = is.na(s2s) | trimws(s2s) == "",
    s2s_known       = s2s %in% c("Y", "N"),

    # -- LBW flag (bwtdis < 2500 g) ------------------------------------------
    lbw = !is.na(bwtdis_num) & bwtdis_num < 2500
  )

# -- Filter: SMCH, intervention period, live births --------------------------
mo_int <- mo %>%
  filter(
    facility == "SMCH",
    adm_date >= INT_START,
    adm_date <= INT_END,
    live_birth == TRUE
  )

# Note: all records in the maternal outcomes dataset are from the SMCH
# delivery log, so all are inborn by definition.
# No separate inborn/outborn filter is applied.

cat(sprintf("\nMaternal outcomes -- SMCH, Nov 2024 - Mar 2026 (extended):\n"))
cat(sprintf("  Total rows (all neotreeoutcome)   : %d\n",
            nrow(mo %>% filter(facility == "SMCH",
                               adm_date >= INT_START, adm_date <= INT_END))))
cat(sprintf("  Live births (LB + ENND)           : %d\n", nrow(mo_int)))
cat(sprintf("  Eligible (resus == N)             : %d (%.1f%% of live births)\n",
            sum(mo_int$eligible, na.rm = TRUE),
            100 * mean(mo_int$eligible, na.rm = TRUE)))
cat(sprintf("  LBW (bwtdis < 2500 g)             : %d (%.1f%% of live births)\n",
            sum(mo_int$lbw, na.rm = TRUE),
            100 * mean(mo_int$lbw, na.rm = TRUE)))

# Raw value distributions
cat("\nDCC raw values (live births, int period):\n")
print(table(dcc = mo_int$dcc, useNA = "ifany"))
cat("\nS2S raw values (live births, int period):\n")
print(table(s2s = mo_int$s2s, useNA = "ifany"))

# Overall uptake (known denominator)
n_mo_dcc_known <- sum(mo_int$dcc_known, na.rm = TRUE)
n_mo_s2s_known <- sum(mo_int$s2s_known, na.rm = TRUE)
cat(sprintf("\nOverall maternal script DCC uptake (n known = %d): %d / %d (%.1f%%)\n",
            n_mo_dcc_known,
            sum(mo_int$dcc_received, na.rm = TRUE), n_mo_dcc_known,
            100 * sum(mo_int$dcc_received, na.rm = TRUE) / n_mo_dcc_known))
cat(sprintf("Overall maternal script S2S uptake (n known = %d): %d / %d (%.1f%%)\n",
            n_mo_s2s_known,
            sum(mo_int$s2s_received, na.rm = TRUE), n_mo_s2s_known,
            100 * sum(mo_int$s2s_received, na.rm = TRUE) / n_mo_s2s_known))

# Among eligible only
mo_elig <- mo_int %>% filter(eligible == TRUE)
n_elig_dcc <- sum(mo_elig$dcc_known, na.rm = TRUE)
n_elig_s2s <- sum(mo_elig$s2s_known, na.rm = TRUE)
cat(sprintf("\nAmong eligible (resus == N, n=%d):\n", nrow(mo_elig)))
cat(sprintf("  DCC uptake (n known = %d): %d / %d (%.1f%%)\n",
            n_elig_dcc,
            sum(mo_elig$dcc_received, na.rm = TRUE), n_elig_dcc,
            100 * sum(mo_elig$dcc_received, na.rm = TRUE) / n_elig_dcc))
cat(sprintf("  S2S uptake (n known = %d): %d / %d (%.1f%%)\n",
            n_elig_s2s,
            sum(mo_elig$s2s_received, na.rm = TRUE), n_elig_s2s,
            100 * sum(mo_elig$s2s_received, na.rm = TRUE) / n_elig_s2s))

# -----------------------------------------------------------------------------
# B2. Maternal outcomes monthly uptake table
# -----------------------------------------------------------------------------
# Primary denominators:
#   n_total    = all live births in month
#   n_eligible = live births meeting eligibility proxy
#   n_dcc_known = eligible with dcc Y or N (excl. U and not_recorded)
#   pct_dcc_uptake = n_dcc_Y / n_dcc_known * 100
#   pct_dcc_unknown = (dcc_U + dcc_not_recorded) / n_eligible * 100
#                      (documentation quality measure)
# The unknown + not_recorded rate captures the proportion of eligible babies
# where DCC/S2S status cannot be determined.
# -----------------------------------------------------------------------------

cat("\n--- Maternal outcomes monthly uptake table (eligible live births) ---\n")

mo_monthly <- mo_int %>%
  filter(eligible == TRUE) %>%
  group_by(month_label) %>%
  summarise(
    n_eligible         = n(),
    # DCC
    n_dcc_Y            = sum(dcc_received,    na.rm = TRUE),
    n_dcc_N            = sum(dcc_not_recv,    na.rm = TRUE),
    n_dcc_U            = sum(dcc_unknown,     na.rm = TRUE),
    n_dcc_not_recorded = sum(dcc_not_recorded, na.rm = TRUE),
    n_dcc_known        = sum(dcc_known,       na.rm = TRUE),
    pct_dcc_uptake     = if_else(n_dcc_known > 0,
                                 round(100 * n_dcc_Y / n_dcc_known, 1),
                                 NA_real_),
    pct_dcc_unknown    = round(100 * (n_dcc_U + n_dcc_not_recorded) / n_eligible, 1),
    # S2S
    n_s2s_Y            = sum(s2s_received,    na.rm = TRUE),
    n_s2s_N            = sum(s2s_not_recv,    na.rm = TRUE),
    n_s2s_U            = sum(s2s_unknown,     na.rm = TRUE),
    n_s2s_not_recorded = sum(s2s_not_recorded, na.rm = TRUE),
    n_s2s_known        = sum(s2s_known,       na.rm = TRUE),
    pct_s2s_uptake     = if_else(n_s2s_known > 0,
                                 round(100 * n_s2s_Y / n_s2s_known, 1),
                                 NA_real_),
    pct_s2s_unknown    = round(100 * (n_s2s_U + n_s2s_not_recorded) / n_eligible, 1),
    .groups = "drop"
  ) %>%
  arrange(month_label)

# Join total live-birth count per month from the unfiltered mo_int
# (mo_monthly above is filtered to eligible only, so n_total_livebirth
# must come from the pre-eligibility dataset).
mo_total_live <- mo_int %>%
  group_by(month_label) %>%
  summarise(n_total_livebirth = n(), .groups = "drop")

mo_monthly <- mo_monthly %>%
  left_join(mo_total_live, by = "month_label") %>%
  relocate(month_label, n_total_livebirth, n_eligible)

print(mo_monthly, n = Inf)
write_csv(mo_monthly, file.path(OUTPUT_DIR, "08b_maternal_monthly_uptake.csv"))
cat("\nSaved: 08b_maternal_monthly_uptake.csv\n")

# -- Cross-check only: broader apgar1>6 OR resus=="N" proxy ------------------
# Not used for the primary table above; kept so the effect of the broader
# proxy on the eligible count can be seen directly.
mo_monthly_broadproxy <- mo_int %>%
  filter(eligible_broadproxy == TRUE) %>%
  group_by(month_label) %>%
  summarise(
    n_eligible_broadproxy = n(),
    n_dcc_Y     = sum(dcc_received, na.rm = TRUE),
    n_dcc_known = sum(dcc_known,    na.rm = TRUE),
    pct_dcc_uptake = if_else(n_dcc_known > 0,
                             round(100 * n_dcc_Y / n_dcc_known, 1), NA_real_),
    n_s2s_Y     = sum(s2s_received, na.rm = TRUE),
    n_s2s_known = sum(s2s_known,    na.rm = TRUE),
    pct_s2s_uptake = if_else(n_s2s_known > 0,
                             round(100 * n_s2s_Y / n_s2s_known, 1), NA_real_),
    .groups = "drop"
  ) %>%
  arrange(month_label)
write_csv(mo_monthly_broadproxy,
          file.path(OUTPUT_DIR, "08h_maternal_monthly_uptake_broadproxy_crosscheck.csv"))
cat("Saved: 08h_maternal_monthly_uptake_broadproxy_crosscheck.csv (cross-check only)\n")

# -- Documentation quality: unknown rate trend --------------------------------
cat("\n--- Documentation quality: DCC unknown rate (U + not_recorded) by month ---\n")
cat("  (Lower % = better recording; Nov 2024 elevated due to intervention ramp-up)\n\n")

doc_quality <- mo_monthly %>%
  select(month_label, n_eligible,
         n_dcc_U, n_dcc_not_recorded, pct_dcc_unknown,
         n_s2s_U, n_s2s_not_recorded, pct_s2s_unknown)
print(doc_quality, n = Inf)

cat("\nDCC unknown/not-recorded rate range:\n")
cat(sprintf("  Min: %.1f%% (%s)\n",
            min(mo_monthly$pct_dcc_unknown, na.rm = TRUE),
            mo_monthly$month_label[which.min(mo_monthly$pct_dcc_unknown)]))
cat(sprintf("  Max: %.1f%% (%s)\n",
            max(mo_monthly$pct_dcc_unknown, na.rm = TRUE),
            mo_monthly$month_label[which.max(mo_monthly$pct_dcc_unknown)]))
cat(sprintf("  Mean across months: %.1f%%\n",
            mean(mo_monthly$pct_dcc_unknown, na.rm = TRUE)))

cat("\nS2S unknown/not-recorded rate range:\n")
cat(sprintf("  Min: %.1f%% (%s)\n",
            min(mo_monthly$pct_s2s_unknown, na.rm = TRUE),
            mo_monthly$month_label[which.min(mo_monthly$pct_s2s_unknown)]))
cat(sprintf("  Max: %.1f%% (%s)\n",
            max(mo_monthly$pct_s2s_unknown, na.rm = TRUE),
            mo_monthly$month_label[which.max(mo_monthly$pct_s2s_unknown)]))
cat(sprintf("  Mean across months: %.1f%%\n",
            mean(mo_monthly$pct_s2s_unknown, na.rm = TRUE)))

# Correlation between unknown rate and uptake rate (is poor documentation
# masking real uptake?)
if (sum(!is.na(mo_monthly$pct_dcc_unknown) &
        !is.na(mo_monthly$pct_dcc_uptake)) > 3) {
  r_dcc <- cor(mo_monthly$pct_dcc_unknown, mo_monthly$pct_dcc_uptake,
               use = "complete.obs")
  cat(sprintf("\nCorrelation (Pearson r) between DCC unknown rate and DCC uptake: %.3f\n",
              r_dcc))
}

# -----------------------------------------------------------------------------
# B3. Maternal outcomes uptake chart
# -----------------------------------------------------------------------------

# Uptake line plot (DCC and S2S)
mo_monthly_long <- mo_monthly %>%
  select(month_label, pct_dcc_uptake, pct_s2s_uptake) %>%
  pivot_longer(
    cols      = c(pct_dcc_uptake, pct_s2s_uptake),
    names_to  = "intervention",
    values_to = "pct"
  ) %>%
  mutate(
    intervention = recode(intervention,
                          pct_dcc_uptake = "DCC",
                          pct_s2s_uptake = "S2S (ESSC)"),
    adm_month    = ym(month_label)
  )

mo_col_map <- c("DCC" = COL_DCC, "S2S (ESSC)" = COL_S2S)

p_mo <- ggplot(mo_monthly_long,
               aes(x = adm_month, y = pct, colour = intervention,
                   group = intervention)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_text(
    data = mo_monthly %>% mutate(adm_month = ym(month_label)),
    aes(x = adm_month, y = 2, label = paste0("n=", n_dcc_known)),
    inherit.aes = FALSE,
    size = 2.8, colour = "grey40", vjust = 0
  ) +
  scale_y_continuous(
    name   = "Uptake rate (%, of known: Y or N)",
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_x_date(
    name        = NULL,
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_colour_manual(
    name   = "Intervention",
    values = mo_col_map
  ) +
  labs(
    title    = "Maternal Outcome Script: DCC and S2S Uptake by Month",
    subtitle = paste0("SMCH eligible live births (resus == N), ",
                      INT_START, " to ", INT_END,
                      "\nDenominator: babies with known status (Y or N; excludes Unknown and not recorded)"),
    caption  = "Note: n = DCC known denominator; uptake = Y / (Y+N)"
  ) +
  theme_gh()

ggsave(file.path(OUTPUT_DIR, "08e_maternal_monthly_uptake.png"),
       plot   = p_mo,
       width  = 10,
       height = 5.5,
       dpi    = 150)
cat("\nSaved: 08e_maternal_monthly_uptake.png\n")

# -- Documentation quality chart (unknown + not-recorded rate) ----------------
mo_unk_long <- mo_monthly %>%
  select(month_label, pct_dcc_unknown, pct_s2s_unknown) %>%
  pivot_longer(
    cols      = c(pct_dcc_unknown, pct_s2s_unknown),
    names_to  = "intervention",
    values_to = "pct_unknown"
  ) %>%
  mutate(
    intervention = recode(intervention,
                          pct_dcc_unknown = "DCC",
                          pct_s2s_unknown = "S2S (ESSC)"),
    adm_month    = ym(month_label)
  )

p_unk <- ggplot(mo_unk_long,
                aes(x = adm_month, y = pct_unknown, colour = intervention,
                    group = intervention)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    name   = "Unknown / not-recorded rate (% of eligible)",
    limits = c(0, max(mo_unk_long$pct_unknown, na.rm = TRUE) * 1.15),
    breaks = seq(0, 100, 10)
  ) +
  scale_x_date(
    name        = NULL,
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_colour_manual(
    name   = "Intervention",
    values = mo_col_map
  ) +
  geom_hline(yintercept = 20, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  annotate("text", x = ym("2025-03"), y = 21,
           label = "20% threshold", colour = "grey40", size = 3, hjust = 0) +
  labs(
    title    = "Documentation Quality: Unknown / Not-Recorded Rate by Month",
    subtitle = paste0("Maternal outcome script, SMCH eligible live births\n",
                      "Unknown (U) + not recorded, as % of eligible babies per month"),
    caption  = paste0("Note: Nov 2024 elevated due to intervention ramp-up (first month).\n",
                      "From Dec 2024 onwards, 'not recorded' = 0; Unknown (U) reflects genuine ambiguity in delivery log.")
  ) +
  theme_gh()

ggsave(file.path(OUTPUT_DIR, "08f_maternal_unknown_rate.png"),
       plot   = p_unk,
       width  = 10,
       height = 5.0,
       dpi    = 150)
cat("Saved: 08f_maternal_unknown_rate.png\n")


# =============================================================================
# PART C -- LBW SUBGROUP (maternal outcome script, bwtdis < 2500 g)
# =============================================================================

cat("\n\n=============================================================\n")
cat("  PART C: LBW Subgroup -- Maternal Outcome Script\n")
cat("=============================================================\n\n")

# -- LBW eligible babies from maternal outcomes script -----------------------
mo_lbw <- mo_int %>%
  filter(lbw == TRUE, eligible == TRUE)

cat(sprintf("LBW eligible babies (bwtdis < 2500 g, resus == N): %d\n",
            nrow(mo_lbw)))
cat(sprintf("  As %% of all eligible live births: %.1f%%\n",
            100 * nrow(mo_lbw) / nrow(mo_elig)))

# Overall LBW uptake
n_lbw_dcc_known <- sum(mo_lbw$dcc_known, na.rm = TRUE)
n_lbw_s2s_known <- sum(mo_lbw$s2s_known, na.rm = TRUE)
cat(sprintf("\nOverall LBW DCC uptake  (n known = %d): %d / %d (%.1f%%)\n",
            n_lbw_dcc_known,
            sum(mo_lbw$dcc_received, na.rm = TRUE), n_lbw_dcc_known,
            100 * sum(mo_lbw$dcc_received, na.rm = TRUE) / n_lbw_dcc_known))
cat(sprintf("Overall LBW S2S uptake  (n known = %d): %d / %d (%.1f%%)\n",
            n_lbw_s2s_known,
            sum(mo_lbw$s2s_received, na.rm = TRUE), n_lbw_s2s_known,
            100 * sum(mo_lbw$s2s_received, na.rm = TRUE) / n_lbw_s2s_known))

# -- LBW monthly uptake table ------------------------------------------------
cat("\n--- LBW monthly DCC and S2S uptake ---\n")

lbw_monthly <- mo_lbw %>%
  group_by(month_label) %>%
  summarise(
    n_lbw_eligible     = n(),
    n_dcc_Y            = sum(dcc_received,    na.rm = TRUE),
    n_dcc_N            = sum(dcc_not_recv,    na.rm = TRUE),
    n_dcc_U            = sum(dcc_unknown,     na.rm = TRUE),
    n_dcc_not_recorded = sum(dcc_not_recorded, na.rm = TRUE),
    n_dcc_known        = sum(dcc_known,       na.rm = TRUE),
    pct_dcc_uptake     = if_else(n_dcc_known > 0,
                                 round(100 * n_dcc_Y / n_dcc_known, 1),
                                 NA_real_),
    pct_dcc_unknown    = round(100 * (n_dcc_U + n_dcc_not_recorded) / n_lbw_eligible, 1),
    n_s2s_Y            = sum(s2s_received,    na.rm = TRUE),
    n_s2s_N            = sum(s2s_not_recv,    na.rm = TRUE),
    n_s2s_U            = sum(s2s_unknown,     na.rm = TRUE),
    n_s2s_not_recorded = sum(s2s_not_recorded, na.rm = TRUE),
    n_s2s_known        = sum(s2s_known,       na.rm = TRUE),
    pct_s2s_uptake     = if_else(n_s2s_known > 0,
                                 round(100 * n_s2s_Y / n_s2s_known, 1),
                                 NA_real_),
    pct_s2s_unknown    = round(100 * (n_s2s_U + n_s2s_not_recorded) / n_lbw_eligible, 1),
    .groups = "drop"
  ) %>%
  arrange(month_label)

print(lbw_monthly, n = Inf)
write_csv(lbw_monthly, file.path(OUTPUT_DIR, "08c_lbw_monthly_uptake.csv"))
cat("\nSaved: 08c_lbw_monthly_uptake.csv\n")

# -- LBW vs all-eligible comparison ------------------------------------------
cat("\n--- LBW vs all-eligible uptake comparison (overall) ---\n")
cat(sprintf("                    DCC uptake      S2S uptake\n"))
cat(sprintf("  All eligible    : %.1f%%            %.1f%%\n",
            100 * sum(mo_elig$dcc_received, na.rm = TRUE) / n_elig_dcc,
            100 * sum(mo_elig$s2s_received, na.rm = TRUE) / n_elig_s2s))
cat(sprintf("  LBW eligible    : %.1f%%            %.1f%%\n",
            100 * sum(mo_lbw$dcc_received, na.rm = TRUE) / n_lbw_dcc_known,
            100 * sum(mo_lbw$s2s_received, na.rm = TRUE) / n_lbw_s2s_known))

# -- LBW monthly chart -------------------------------------------------------
lbw_monthly_long <- lbw_monthly %>%
  select(month_label, pct_dcc_uptake, pct_s2s_uptake) %>%
  pivot_longer(
    cols      = c(pct_dcc_uptake, pct_s2s_uptake),
    names_to  = "intervention",
    values_to = "pct"
  ) %>%
  mutate(
    intervention = recode(intervention,
                          pct_dcc_uptake = "DCC",
                          pct_s2s_uptake = "S2S (ESSC)"),
    adm_month    = ym(month_label)
  )

p_lbw <- ggplot(lbw_monthly_long,
                aes(x = adm_month, y = pct, colour = intervention,
                    group = intervention)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_text(
    data = lbw_monthly %>% mutate(adm_month = ym(month_label)),
    aes(x = adm_month, y = 2, label = paste0("n=", n_dcc_known)),
    inherit.aes = FALSE,
    size = 2.8, colour = "grey40", vjust = 0
  ) +
  scale_y_continuous(
    name   = "Uptake rate (%, of known: Y or N)",
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_x_date(
    name        = NULL,
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_colour_manual(
    name   = "Intervention",
    values = mo_col_map
  ) +
  labs(
    title    = "LBW Subgroup: DCC and S2S Uptake by Month",
    subtitle = paste0("Maternal outcome script, SMCH; LBW eligible babies (bwtdis < 2500 g,\n",
                      "resus == N), ", INT_START, " to ", INT_END,
                      "\nDenominator: LBW babies with known status (Y or N)"),
    caption  = "Note: n = DCC known denominator (LBW subset)"
  ) +
  theme_gh()

ggsave(file.path(OUTPUT_DIR, "08g_lbw_monthly_uptake.png"),
       plot   = p_lbw,
       width  = 10,
       height = 5.5,
       dpi    = 150)
cat("Saved: 08g_lbw_monthly_uptake.png\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n=============================================================\n")
cat("  SUMMARY\n")
cat("=============================================================\n\n")

cat("Part A -- NeoTree (inborn neonatal admissions, SMCH, Nov 2024 - Mar 2026)\n")
cat(sprintf("  Total admissions           : %d\n", nrow(int_pop)))
cat(sprintf("  delivinter known           : %d (%.1f%%)\n",
            sum(int_pop$delivinter_known, na.rm = TRUE),
            100 * mean(int_pop$delivinter_known, na.rm = TRUE)))
cat(sprintf("  delivinter missing         : %d (%.1f%%)\n",
            sum(int_pop$delivinter_missing, na.rm = TRUE),
            100 * mean(int_pop$delivinter_missing, na.rm = TRUE)))
cat(sprintf("  DCC received (of known)    : %d (%.1f%%)\n",
            sum(nt_known$dcc_received, na.rm = TRUE),
            100 * mean(nt_known$dcc_received, na.rm = TRUE)))
cat(sprintf("  S2S received (of known)    : %d (%.1f%%)\n",
            sum(nt_known$s2s_received, na.rm = TRUE),
            100 * mean(nt_known$s2s_received, na.rm = TRUE)))

cat("\nPart B -- Maternal outcome script (SMCH deliveries, Nov 2024 - Mar 2026)\n")
cat(sprintf("  Live births                : %d\n", nrow(mo_int)))
cat(sprintf("  Eligible live births       : %d (%.1f%%)\n",
            nrow(mo_elig),
            100 * nrow(mo_elig) / nrow(mo_int)))
cat(sprintf("  DCC received (eligible, known denom %d): %.1f%%\n",
            n_elig_dcc,
            100 * sum(mo_elig$dcc_received, na.rm = TRUE) / n_elig_dcc))
cat(sprintf("  S2S received (eligible, known denom %d): %.1f%%\n",
            n_elig_s2s,
            100 * sum(mo_elig$s2s_received, na.rm = TRUE) / n_elig_s2s))

cat("\nPart C -- LBW subgroup (maternal script, bwtdis < 2500 g, eligible)\n")
cat(sprintf("  LBW eligible babies        : %d\n", nrow(mo_lbw)))
cat(sprintf("  DCC received (known denom %d): %.1f%%\n",
            n_lbw_dcc_known,
            100 * sum(mo_lbw$dcc_received, na.rm = TRUE) / n_lbw_dcc_known))
cat(sprintf("  S2S received (known denom %d): %.1f%%\n",
            n_lbw_s2s_known,
            100 * sum(mo_lbw$s2s_received, na.rm = TRUE) / n_lbw_s2s_known))

cat("\nOutputs saved to:", OUTPUT_DIR, "\n")
cat("  08a_neotree_monthly_uptake.csv\n")
cat("  08b_maternal_monthly_uptake.csv\n")
cat("  08c_lbw_monthly_uptake.csv\n")
cat("  08d_neotree_monthly_uptake.png\n")
cat("  08e_maternal_monthly_uptake.png\n")
cat("  08f_maternal_unknown_rate.png\n")
cat("  08g_lbw_monthly_uptake.png\n")

cat("\n=== Script 08 complete ===\n")

# -- Close log ----------------------------------------------------------------
sink(type = "output")
close(log_con)
message(sprintf("Log saved to: %s", LOG_FILE))
