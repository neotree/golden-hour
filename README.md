# Golden Hour intervention at Sally Mugabe Central Hospital — analysis code

Analysis scripts for the evaluation of a **Golden Hour** quality-improvement
intervention (delayed cord clamping and early skin-to-skin care) in the neonatal
service at **Sally Mugabe Central Hospital (SMCH)**, Harare, Zimbabwe, using
routinely collected **NeoTree** and maternal delivery-log data.

This repository is provided **for transparency and reproducibility** — it shows
exactly how every result, table and figure in the manuscript was produced. It
contains **code only**. No patient data are included (see
[Data availability](#data-availability)).

---

## What the analysis does

The intervention began in **November 2024**. It is evaluated mainly with an
**interrupted time series (ITS)** design, supported by seasonality sensitivity
analyses and an individual-level (per-baby) comparison. The scripts cover:

- **Intervention coverage / uptake** of DCC and skin-to-skin among eligible
  live births, with an explicit yes / no / unknown breakdown.
- **Primary outcome** — monthly % of admissions with a normal admission
  temperature (36.5–37.5 °C), by segmented-regression ITS, plus two
  **seasonality** checks (Fourier-adjusted ITS and same-season year-on-year).
- **Neonatal mortality** over time (ITS), overall and by birthweight, with a
  Fourier-adjusted sensitivity analysis.
- **Per-baby mortality odds ratios** by intervention receipt, **unadjusted and
  adjusted**, stratified by birthweight.
- **Supporting analyses** — descriptive overview, power analysis, and a
  data-completeness check of the arrival-to-unit temperature field.

A line-by-line mapping of each manuscript table and figure to the script and
output file that produced it is in **[FIGURE_TABLE_SOURCES.md](FIGURE_TABLE_SOURCES.md)**.

---

## Repository structure

```
neotree-golden-hour/
├── README.md                  # this file
├── LICENSE                    # MIT licence (code)
├── CITATION.cff               # how to cite this repository
├── FIGURE_TABLE_SOURCES.md    # provenance map: manuscript element -> script -> output
├── .gitignore                 # excludes data, outputs, R session files
└── scripts/                   # all R analysis scripts
    ├── 01_descriptive_overview.R
    ├── 02_outcome_analysis.R          # RETIRED (superseded by 13); kept for the record
    ├── 03_temperature_trends.R
    ├── 04_its_analysis.R
    ├── 04b_its_analysis_excl2022.R    # primary temperature ITS
    ├── 04c_its_seasonality_correction.R
    ├── 04d_seasonal_comparison.R
    ├── 05_power_analysis.R
    ├── 06_its_results_figure.R
    ├── 07_power_by_preperiod.R
    ├── 08_maternal_script_uptake.R
    ├── 09_uptake_status_breakdown.R   # uptake yes/no/unknown
    ├── 12_mortality_its.R
    ├── 12b_mortality_its_fourier.R
    ├── 13_mortality_or_by_bw.R        # per-baby OR, unadjusted + adjusted
    ├── 14_arrival_temp_completeness.R
    ├── 15_sepsis_its.R                # exploratory; not reported (see manuscript)
    └── SCRIPT_INDEX.txt               # detailed, per-script documentation
```

`scripts/SCRIPT_INDEX.txt` gives a full description of every script (purpose,
population, model, and outputs). Start there for detail.

---

## Data availability

**No data are stored in this repository.** The analyses use two routinely
collected, patient-level datasets that cannot be shared publicly for
confidentiality and data-governance reasons:

| Expected file (place in `00-DATA/`) | Description |
|---|---|
| `ZIM_db_master_joined_to_20260401.csv` | **All** NeoTree neonatal-unit admissions; each admission joined to its discharge record where one exists. Admissions and discharges are linked **deterministically** (exact match on unique identifier + facility), not by probabilistic matching. ~91% of SMCH inborn admissions have a matched discharge. |
| `zim_db_maternal_outcomes_20260501_cleaned.csv` | Maternal delivery-log extract (all SMCH births), with DCC/skin-to-skin status. |

The master file therefore contains every admission, whether or not a discharge
record was found. Analyses use it accordingly: **admission-side outcomes**
(admission temperature) use **all** admissions, while **discharge-side outcomes**
(neonatal mortality) are necessarily restricted to **matched** records
(`match_type == "direct_match"`), because the discharge outcome exists only for
matched admissions.

These data are held by the NeoTree programme / SMCH and are **available from the
data custodians on reasonable request, subject to data-sharing agreement and the
relevant ethics/governance approvals**. The variable definitions used by the
scripts are documented inline and in `SCRIPT_INDEX.txt`.

> A ready-to-paste **Data and code availability statement** for the manuscript is
> at the bottom of this file.

---

## How to run

1. **Install R** (the analyses were run on **R 4.6**) and the required packages:

   ```r
   install.packages(c("tidyverse", "lubridate", "scales", "broom",
                      "lmtest", "sandwich", "patchwork", "cowplot"))
   ```

2. **Place the two data files** (above) in a folder named `00-DATA/` at the
   repository root, i.e. alongside `scripts/`. Outputs (CSVs, figures, logs) are
   written to `03-OUTPUTS/`, which is created automatically. The scripts locate
   their own directory at run time and use paths relative to it, so no editing is
   needed if this layout is kept:

   ```
   neotree-golden-hour/
   ├── scripts/        # the R files
   ├── 00-DATA/        # you provide the two CSVs here (not in the repo)
   └── 03-OUTPUTS/     # created automatically when scripts run
   ```

3. **Run the scripts** from a terminal, e.g.:

   ```bash
   cd scripts
   Rscript 01_descriptive_overview.R
   ```

   Suggested order (only script 06 has a dependency — run 04 first):

   ```
   01_descriptive_overview.R
   03_temperature_trends.R
   04_its_analysis.R
   04b_its_analysis_excl2022.R       # primary temperature ITS
   04c_its_seasonality_correction.R
   04d_seasonal_comparison.R
   05_power_analysis.R
   06_its_results_figure.R           # requires 04 to have run
   07_power_by_preperiod.R
   08_maternal_script_uptake.R
   09_uptake_status_breakdown.R
   12_mortality_its.R
   12b_mortality_its_fourier.R
   13_mortality_or_by_bw.R
   14_arrival_temp_completeness.R
   15_sepsis_its.R                   # exploratory; not reported
   ```

   `02_outcome_analysis.R` is **retired** (superseded by `13`) and is not part of
   the reported pipeline; it is retained only for completeness.

Each script is self-contained, ASCII-only, and writes a `*_log.txt` alongside its
numeric outputs so a full record of every run is captured.

---

## Notes on methods

- **Eligibility** for the intervention is a proxy for not needing immediate
  resuscitation (1-minute Apgar > 6, or resuscitation recorded as none).
- **Primary temperature window** is January 2023 – March 2026 (2022 excluded for
  data-quality reasons; a January-2022 sensitivity analysis is included).
- The **per-baby odds ratios are observational** and subject to indication bias;
  they are reported as supportive, with the ITS as the primary evaluation.
- A **sepsis** ITS (`15`) was explored but is **not reported** in the manuscript
  (incomplete recent data and confounding by a concurrent sepsis-diagnosis QI
  programme).

---

## How to cite

See [CITATION.cff](CITATION.cff). Please cite the repository (and its release/DOI,
once minted — see below) in the manuscript's code-availability statement.

---

## Data and code availability statement (ready to paste)

> **Code availability.** All analysis code is openly available at
> https://github.com/<your-org-or-user>/neotree-golden-hour (archived at
> [Zenodo DOI — to be added on release]). The repository contains the complete set
> of R scripts used to generate every result, table and figure, together with a
> provenance map linking each manuscript item to its source script.
>
> **Data availability.** The individual-level NeoTree and maternal delivery-log
> data that support these findings are not publicly available owing to patient
> confidentiality and data-governance restrictions. De-identified data are
> available from the NeoTree programme / Sally Mugabe Central Hospital data
> custodians on reasonable request and subject to a data-sharing agreement and
> the relevant ethical approvals.

*(Replace the GitHub URL, add the Zenodo DOI after creating a release, and adjust
the custodian wording to match your governance arrangements.)*
