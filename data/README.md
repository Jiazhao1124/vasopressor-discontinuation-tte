# Data

No patient-level data are distributed with this repository. MIMIC-IV and eICU-CRD are
available to credentialed users on PhysioNet under the PhysioNet Credentialed Health Data
Use Agreement, which prohibits redistribution.

Run the SQL in `sql/` against your own PostgreSQL build of the databases and save the
outputs as:

| File | Query |
|---|---|
| `raw/mimic_vs.csv` | `sql/mimiciv/01_cohort_features.sql` |
| `raw/mimic_vx.csv` | `sql/mimiciv/02_extra_covariates.sql` |
| `raw/mimic_cs.csv` | `sql/mimiciv/03_code_status.sql` |
| `raw/mimic_note.csv` | `sql/mimiciv/04_discharge_note_keywords.sql` |
| `raw/eicu_vs.csv` | `sql/eicu/01_cohort_features.sql` |

`derived/` is written by `analysis/01_prepare_mimic.py` and `analysis/02_prepare_eicu.py`.
