# Abrupt versus tapered vasopressor discontinuation: target trial emulation

Analysis code for:

> Jing S, Wen L, Song J. Abrupt versus tapered discontinuation of vasopressors in critically ill adults: a target trial emulation.

The study compares discontinuation of vasopressors at a norepinephrine-equivalent dose (NEE) above 0.05 µg/kg/min (abrupt) with discontinuation at 0.05 µg/kg/min or less (tapered). The primary outcome is vasopressor re-initiation within 24 hours. MIMIC-IV is the derivation cohort and eICU-CRD the replication cohort.

## Data access

No patient-level data are included. Both databases are available to credentialed users on [PhysioNet](https://physionet.org):

- MIMIC-IV v3.0, including the `mimiciv_derived` concept tables (norepinephrine-equivalent dose, SOFA, Sepsis-3) and MIMIC-IV-Note v2.2
- eICU Collaborative Research Database v2.0 (schema `eicuii` in the queries; change it to match your build)

Run the queries in `sql/` and save the outputs in `data/raw/` as listed in [`data/README.md`](data/README.md).

## Reproducing the results

```bash
pip install -r requirements.txt
bash run_all.sh
```

The scripts reproduce the published estimates exactly when they are run with the package versions in `requirements.txt`. All random steps use fixed seeds.

| Script | Purpose | Output in `results/` | Manuscript |
|---|---|---|---|
| `analysis/01_prepare_mimic.py` | MIMIC-IV eligibility and analytic cohort | `flow_mimic.csv` | Figure S1 |
| `analysis/02_prepare_eicu.py` | eICU-CRD eligibility and analytic cohort | `flow_eicu.csv` | Figure S1 |
| `analysis/03_main_analysis.py` | Overlap-weighted primary and secondary outcomes, sensitivity analyses, subgroups, eICU-CRD replication | `main_results.csv` | Table 2, Table S3, Figure S2 |
| `analysis/04_competing_risk.py` | Aalen–Johansen cumulative incidence with death as a competing event | `cumulative_incidence_mimic.csv` | Results |
| `analysis/05_dose_response.py` | Spline g-computation of risk by NEE at discontinuation | `dose_response_mimic.csv`, `dose_response_eicu.csv` | Figure 1 |
| `analysis/06_descriptives.py` | Crude risks, timing of re-initiation, hospital-level consistency | `descriptives.csv` | Results |
| `analysis/07_table1.py` | Baseline characteristics of both cohorts | `table1.csv`, `table1_missing.csv` | Table 1 |
| `analysis/08_psm_balance.py` | 1:1 propensity matching and covariate balance | `tableS4_psm_balance.csv`, `psm_matched_effects.csv`, `tableS2_missing_by_group.csv` | Tables S2, S4 |
| `figures/*` | Figures (R and matplotlib) | `figures/output/` | Figures 1, S1, S2 |

## Methods in brief

- **Time zero:** end of the first vasopressor episode lasting at least 6 hours (infusion gaps ≤60 min in MIMIC-IV, ≤120 min in eICU-CRD), in patients alive and in the ICU 1 hour later.
- **Exclusions:** peak NEE below 0.05 µg/kg/min; comfort care charted up to 2 hours after discontinuation (MIMIC-IV) or care limitation records up to 24 hours after (eICU-CRD); comfort care charted within 24 hours in MIMIC-IV patients who died within 24 hours without re-initiation; comfort care or hospice in the MIMIC-IV discharge summary; NEE at discontinuation above 0.2 µg/kg/min.
- **Confounding:** overlap weights from a cross-fitted gradient-boosted propensity score (5-fold, grouped by patient in MIMIC-IV and by hospital in eICU-CRD, clipped at 0.01–0.99).
- **Inference:** cluster bootstrap with 500 resamples (patient in MIMIC-IV; patient and hospital in eICU-CRD).
- **Sensitivity analyses:** logistic overlap weights, trimmed IPTW, AIPW, 1:1 propensity matching (caliper 0.2 SD of the logit), first stay per patient, alternative thresholds of 0.02 to 0.10 µg/kg/min, re-initiation beyond 2 hours, and E-value.

## Repository layout

```
sql/        extraction queries (PostgreSQL)
analysis/   cohort construction and statistical analysis (Python)
figures/    figure scripts (R, matplotlib)
results/    aggregate results only
data/       empty; place query outputs in data/raw/
```

## License

Code is released under the MIT License. The data remain subject to the PhysioNet Credentialed Health Data Use Agreement.

## Contact

Jiazhao Song, Cixi People's Hospital, Wenzhou Medical University, and University Hospital Erlangen.
