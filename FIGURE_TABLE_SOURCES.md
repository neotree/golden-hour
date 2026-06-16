# Provenance map — manuscript elements to source code

Each manuscript table and figure below is linked to the R script that produces it
and the key output file(s) it writes to `03-OUTPUTS/`. This lets co-authors and
reviewers trace every reported number to the exact code that generated it.
Output filenames refer to files created when the scripts are run against the
source data.

## Tables

| Manuscript item | Source script | Key output file(s) |
|---|---|---|
| Table 1 — Analysis populations / cohort counts | `01_descriptive_overview.R`; `12_mortality_its.R`; `13_mortality_or_by_bw.R` | `01c_key_counts_summary.csv`; `12a_mortality_monthly.csv`; `13c_mortality_or_counts.csv` |
| Table 2 — DCC / ESSC status (yes / no / unknown), all + LBW | `09_uptake_status_breakdown.R` | `09a_dcc_all_eligible.csv`, `09b_essc_all_eligible.csv`, `09c_dcc_lbw_eligible.csv`, `09d_essc_lbw_eligible.csv` |
| Table 3 — Temperature ITS coefficients (+ seasonality) | `04b_its_analysis_excl2022.R`; `04c_its_seasonality_correction.R`; `04d_seasonal_comparison.R` | `04b_its_coefficients.csv`; `04c_its_seasonality_correction_coefficients.csv`; `04d_seasonal_comparison_monthly.csv` |
| Table 4 — Mortality ITS coefficients (+ Fourier) | `12_mortality_its.R`; `12b_mortality_its_fourier.R` | `12b_its_coefficients.csv`; `12b_mortality_fourier_coefficients.csv` |
| Table 5 — Per-baby mortality OR, unadjusted vs adjusted, by BW | `13_mortality_or_by_bw.R` | `13b_mortality_or_summary.csv` (headline); `13a_mortality_or_full_models.csv` (all coefficients); `13c_mortality_or_counts.csv` |

## Figures

| Manuscript item | Source script | Key output file(s) |
|---|---|---|
| Uptake status by month (stacked yes/no/unknown, all + LBW) | `09_uptake_status_breakdown.R` | `09e_status_breakdown_plot.png` |
| Uptake / coverage trend over time | `08_maternal_script_uptake.R` | `08e_maternal_monthly_uptake.png` |
| Temperature ITS time series | `04b_its_analysis_excl2022.R` (or `06_its_results_figure.R` for the publication figure) | `04b_its_plot.png`; `06a_its_timeseries.png`, `06b_its_forest.png` |
| Temperature seasonality (Fourier / same-season) | `04c_its_seasonality_correction.R`; `04d_seasonal_comparison.R` | `04c_its_seasonality_correction_plot.png`; `04d_seasonal_comparison_plot.png` |
| Mortality ITS — all inborn and LBW | `12_mortality_its.R` | `12c_mortality_its_plot.png`; `12d_lbw_mortality_its_plot.png` |
| Mortality Fourier sensitivity | `12b_mortality_its_fourier.R` | `12b_mortality_fourier_plot_all.png`, `12b_mortality_fourier_plot_lbw.png` |
| Arrival-to-unit temperature completeness | `14_arrival_temp_completeness.R` | `14c_arrival_temp_completeness_plot.png`; `14a/14b_*.csv` |

## Supporting / not reported

| Item | Source script | Notes |
|---|---|---|
| Descriptive baseline & variable completeness | `01_descriptive_overview.R` | `01a`–`01e_*.csv` |
| Statistical power (ITS) | `05_power_analysis.R`; `07_power_by_preperiod.R` | Justifies the Jan-2023 primary window |
| Sepsis ITS | `15_sepsis_its.R` | **Exploratory; not reported** (incomplete recent data + confounding) |
| Original exposure/outcome logistic models | `02_outcome_analysis.R` | **Retired**; superseded by `13` |
