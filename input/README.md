# Input data

Place the following 3 files in this folder. They are confidential
patient-level extracts and are **not included in this repository** — see the
main [README.md](../README.md#data-availability) for how to request them.

| Filename | Used by | Description |
|---|---|---|
| `ZIM_db_master_joined_to_20260525.csv` | `01`, `03`, `16`, `17`, `18`, `19` (also cross-reads the 20260401 file), `20` (via `18`'s output) | Neotree master file, all SMCH neonatal-unit admissions joined to discharge record where matched. Every script that reads the master file now uses this later extract — see next row for why `19` still needs the older one too |
| `ZIM_db_master_joined_to_20260401.csv` | `19` only | An earlier snapshot of the master file, kept solely for `19_reference_group_sensitivity.R`, which deliberately compares this extract against the 20260525 one to quantify the reference-group sensitivity quoted in the manuscript's Limitations paragraph |
| `zim_db_maternal_outcomes_20260501_cleaned.csv` | `08`, `23`, `24`, `25`, `27` | Maternal delivery-log extract (all SMCH births), with DCC/skin-to-skin status |

`run_all.R` checks for all 3 files before running anything, and stops with a
clear error listing exactly what's missing if any are absent.
