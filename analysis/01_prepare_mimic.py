"""Build the MIMIC-IV analytic cohort from the four query outputs.

Input : data/raw/mimic_vs.csv, mimic_vx.csv, mimic_cs.csv, mimic_note.csv
Output: data/derived/mimic_analytic.csv, results/flow_mimic.csv
"""
import numpy as np
import pandas as pd
from common import RAW, DERIVED, RESULTS, THRESHOLD, MAX_STOP_DOSE

d = pd.read_csv(RAW / "mimic_vs.csv").merge(pd.read_csv(RAW / "mimic_vx.csv"), on="stay_id")
sed = ["propofol_on", "fentanyl_on", "dexmed_on", "midaz_on"]
d[sed] = d[sed].fillna(0)
d["sedation_any"] = d[sed].max(axis=1)
d["cu"] = d.first_careunit.astype("category").cat.codes
d = d.merge(pd.read_csv(RAW / "mimic_cs.csv"), on="stay_id")

d["abrupt"] = (d.last_d > THRESHOLD).astype(int)
# death within 24 h without vasopressor re-initiation
d["d24nr"] = ((d.death_h <= 24) & ~(d.restart_gap_h <= 24)).astype(int)

flow = [("ICU stays with a qualifying first vasopressor episode", len(d))]
# Exclude comfort care charted up to 2 h after discontinuation, peak NEE < 0.05, and
# death within 24 h without re-initiation after comfort care charted within 24 h
# (treated as withdrawal of support rather than weaning).
el = ((d.cmo_by_stop2h == 0) & (d.max_d >= 0.05)
      & ~((d.d24nr == 1) & (d.cmo_by_24h == 1)))
d = d[el].copy()
flow.append(("Peak NEE >= 0.05 and no comfort care at discontinuation", len(d)))

d = d.merge(pd.read_csv(RAW / "mimic_note.csv"), on="stay_id")
d = d[(d.cmo_note == 0) & (d.hospice_note == 0)]
flow.append(("No comfort care or hospice in discharge summary", len(d)))
d = d[d.last_d <= MAX_STOP_DOSE].reset_index(drop=True)
flow.append(("NEE at discontinuation <= 0.2 (analytic cohort)", len(d)))
flow += [("  tapered (NEE <= 0.05)", int((d.abrupt == 0).sum())),
         ("  abrupt (NEE > 0.05)", int(d.abrupt.sum())),
         ("  unique patients", d.subject_id.nunique())]

d.to_csv(DERIVED / "mimic_analytic.csv", index=False, float_format="%.17g")
pd.DataFrame(flow, columns=["step", "n"]).to_csv(RESULTS / "flow_mimic.csv", index=False)
for s, n in flow:
    print(f"{s:60s} {n:>7,}")
