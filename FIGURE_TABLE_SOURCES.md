# Provenance map — manuscript elements to source code

Each manuscript table and figure below is linked to the R script that produces
it and the key output file(s) it writes under `output/`. This lets co-authors
and reviewers trace every reported number to the exact code that generated it.
Mapped against `The_Golden_Hour_20260821.docx`.

Script numbers are **not renumbered/sequential** — they keep their original
numbers from the analysis working folder (e.g. `01`, `03`, `08`, `16`–`27`) so
they stay traceable to the request emails and provenance notes each was built
from.

## Tables

| Manuscript item | Source script | Key output file(s) |
|---|---|---|
| Table 1 — Documentation by month of DCC/ESSC status, all live births | `25_documentation_all_babies.R` | `25a_documentation_all_vs_eligible_12mo.csv` |
| Table 1a — DCC among eligible LBW infants, by month | `27_lbw_monthly_uptake.R` | `27a_lbw_monthly_uptake_table1def.csv` |
| Table 1b — ESSC among eligible LBW infants, by month | `27_lbw_monthly_uptake.R` | `27b_lbw_monthly_uptake_table2def.csv` |
| Table 2a — DCC among eligible caesarean births, by month | `24_csection_monthly_uptake.R` | `24a_csection_monthly_uptake_12mo.csv` |
| Table 2b — ESSC among eligible caesarean births, by month | `24_csection_monthly_uptake.R` | `24a_csection_monthly_uptake_12mo.csv` |
| Table 3a — Odds ratios for death before discharge, by intervention group | `18_outcome_or_and_sensitivity.R` | `18a_or_full_models.csv`, `18b_or_table4.csv` |
| Table 4 — Uptake by mode of delivery, birth weight, gestational age | `23_table1_maternal_uptake_12mo.R` | `23a_table1_maternal_uptake_12mo.csv` |
| Table 5 — Characteristics of the eligible NICU cohort, by intervention group | `17_eligible_cohort_tables.R` | `17a_table2_characteristics.csv`, `17c_table3_outcomes.csv` |
| Table 6 — Sensitivity analysis, 1500-2499g subgroup | `18_outcome_or_and_sensitivity.R` | `18c_or_counts.csv` |
| Appendix — cohort-flow numbers (22,036 -> 4,860 -> 3,090 -> 1,794 -> 1,654) | `16_data_readiness_checks.R`, `17_eligible_cohort_tables.R` | `16a_completeness.csv`, `16b_readiness_summary.csv` |

## Figures

| Manuscript item | Source script | Key output file(s) |
|---|---|---|
| Figure 2 — DCC/ESSC uptake trend (all / LBW / theatre) | `08_maternal_script_uptake.R` (all-eligible line); `27_lbw_monthly_uptake.R` (LBW line); `24_csection_monthly_uptake.R` (theatre/caesarean line) | `08b_maternal_monthly_uptake.csv`, `27a/27b_*.csv`, `24a_*.csv` |
| Figure 3 — Admission temperature, monthly %, descriptive | `03_temperature_trends.R` | `03a_monthly_normal_temp.csv`, `03_temperature_trends.png` |
| Figure 3b — Forest plot of odds ratios for death before discharge | `20_forest_table4.R` | `20_forest_table4.png` |

Figure 2 is a composite built from three scripts' monthly series (all-eligible,
LBW subgroup, caesarean subgroup) — there is no single script that produces
the combined 3-line figure directly.

## Supporting (not directly cited by a table/figure number, but needed to reproduce reported numbers)

| Item | Source script | Notes |
|---|---|---|
| Descriptive baseline overview | `01_descriptive_overview.R` | Included for general documentation of the source data; not tied to a specific manuscript table |
| Data-readiness / cohort-eligibility pre-check | `16_data_readiness_checks.R` | Confirms the variables the eligible-cohort spec depends on are present; also the source of the Results cohort-flow numbers (see Table above) |
