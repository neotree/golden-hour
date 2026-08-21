# Input data

Place the following 3 files in this folder. They are confidential
patient-level extracts and are **not included in this repository** — see the
main [README.md](../README.md#data-availability) for how to request them.

| Filename | Used by | Description |
|---|---|---|
| `ZIM_db_master_joined_to_20260401.csv` | `01`, `03` | Neotree master file, all SMCH neonatal-unit admissions joined to discharge record where matched |
| `ZIM_db_master_joined_to_20260525.csv` | `16`, `17`, `18`, `19` (also cross-reads the 20260401 file), `20` (via `18`'s output) | Later Neotree master extract — used for the eligible-cohort mortality analysis. `19_reference_group_sensitivity.R` quantifies and explains why two dated extracts are both needed (see the manuscript's Limitations paragraph) |
| `zim_db_maternal_outcomes_20260501_cleaned.csv` | `08`, `23`, `24`, `25`, `27` | Maternal delivery-log extract (all SMCH births), with DCC/skin-to-skin status |

`run_all.R` checks for all 3 files before running anything, and stops with a
clear error listing exactly what's missing if any are absent.
