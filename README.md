# Golden Hour intervention at Sally Mugabe Central Hospital — analysis code

Analysis scripts for the evaluation of a **Golden Hour** quality-improvement
intervention (delayed cord clamping and early skin-to-skin care) in the neonatal
service at **Sally Mugabe Central Hospital (SMCH)**, Harare, Zimbabwe, using
routinely collected **Neotree** and maternal delivery-log data.

This repository is provided **for transparency and reproducibility** — it shows
exactly how every result, table and figure in the manuscript was produced. It
contains **code only**. No patient data are included (see
[Data availability](#data-availability)).

---

## What the analysis does

The intervention began in **November 2024** and ran for a 12-month key study
period (1 November 2024 – 31 October 2025). The scripts cover:

- **Intervention coverage / uptake** of DCC and skin-to-skin care among
  eligible live births — overall, monthly, and by low-birth-weight and
  caesarean-birth subgroups.
- **Admission temperature** — descriptive monthly trend.
- **Neonatal mortality** for a subset of infants admitted to the neonatal
  unit — logistic-regression odds ratios by intervention group (unadjusted,
  adjusted for birth weight, adjusted for birth weight and sex), with a
  1500–2499g sensitivity analysis and a reference-group sensitivity check.

Eligibility for both the uptake analysis and the mortality analysis is a
1-minute Apgar score above 6 (resuscitation type was not used, as it was too
incompletely recorded). Full method detail is in the manuscript's Statistical
Analysis Methods appendix.

A line-by-line mapping of each manuscript table and figure to the script and
output file that produced it is in
**[FIGURE_TABLE_SOURCES.md](FIGURE_TABLE_SOURCES.md)**.

---

## Repository structure

```
golden-hour/
├── README.md                  # this file
├── LICENSE                    # MIT licence (code)
├── CITATION.cff               # how to cite this repository
├── FIGURE_TABLE_SOURCES.md    # provenance map: manuscript element -> script -> output
├── run_all.R                  # runs every script end-to-end
├── .gitignore                 # excludes data, outputs, R session files
├── scripts/                   # all 11 R analysis scripts
│   ├── 01_descriptive_overview.R
│   ├── 03_temperature_trends.R
│   ├── 08_maternal_script_uptake.R
│   ├── 16_data_readiness_checks.R
│   ├── 17_eligible_cohort_tables.R
│   ├── 18_outcome_or_and_sensitivity.R
│   ├── 20_forest_table4.R
│   ├── 23_table1_maternal_uptake_12mo.R
│   ├── 24_csection_monthly_uptake.R
│   ├── 25_documentation_all_babies.R
│   └── 27_lbw_monthly_uptake.R
├── input/                     # you provide the 2 confidential data files here
│   └── README.md              # exact filenames expected
└── output/                    # created automatically when run_all.R runs
```

Script numbers keep their original numbering from the analysis working
folder (not renumbered sequentially) so they stay traceable to request
emails and provenance notes — see `FIGURE_TABLE_SOURCES.md` for why the gaps
exist.

---

## Data availability

**No data are stored in this repository.** The analyses use three routinely
collected, patient-level extracts that cannot be shared publicly for
confidentiality and data-governance reasons — see
**[input/README.md](input/README.md)** for the exact filenames expected.

These data are held by the Neotree programme / SMCH and are **available from
the data custodians on reasonable request, subject to data-sharing agreement
and the relevant ethics/governance approvals**.

> A ready-to-paste **Data and code availability statement** for the manuscript
> is at the bottom of this file.

---

## How to run

1. **Install R** (the analyses were run on **R 4.6**, tidyverse 2.0.0).
2. **Place the 2 required data files** in `input/` — see
   [input/README.md](input/README.md) for exact filenames.
3. From the repository root:

   ```bash
   Rscript run_all.R
   ```

   `run_all.R` installs any missing R packages automatically, checks `input/`
   has everything it needs (and stops with a clear message listing exactly
   what's missing if not), runs all 11 scripts in dependency order, and
   writes every output under `output/`. It continues past a script that
   errors rather than aborting the whole run, and prints a succeeded/failed
   summary at the end.

   Each script also remains fully runnable on its own — see the header
   comment in `run_all.R` for how the `input`/`output` redirection works if
   you want to run a single script directly from `scripts/`.

Each script is self-contained, ASCII-only, and writes a `*_log.txt` alongside
its numeric outputs so a full record of every run is captured.

---

## Notes on methods

- **Eligibility** for both the uptake and mortality analyses is a 1-minute
  Apgar score above 6. Resuscitation type was assessed as a possible
  additional exclusion but not used, as it was incompletely recorded
  (present for 16.0% of the cohort).
- All scripts read a single Neotree master extract,
  `ZIM_db_master_joined_to_20260525.csv` — see `input/README.md`.
- The **mortality odds ratios are observational** and the study was not
  powered to detect a mortality effect; all reported confidence intervals
  cross 1.

---

## How to cite

See [CITATION.cff](CITATION.cff).

---

## Data and code availability statement (ready to paste)

> **Code availability.** All analysis code is openly available at
> https://github.com/neotree/golden-hour. The repository contains the complete
> set of R scripts used to generate every result, table and figure, together
> with a provenance map linking each manuscript item to its source script.
>
> **Data availability.** The individual-level Neotree and maternal delivery-log
> data that support these findings are not publicly available owing to patient
> confidentiality and data-governance restrictions. De-identified data are
> available from the Neotree programme / Sally Mugabe Central Hospital data
> custodians on reasonable request and subject to a data-sharing agreement and
> the relevant ethical approvals.
