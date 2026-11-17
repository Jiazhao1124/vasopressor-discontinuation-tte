"""Build the eICU-CRD analytic cohort.

Input : data/raw/eicu_vs.csv
Output: data/derived/eicu_analytic.csv, results/flow_eicu.csv
"""
import pandas as pd
from common import RAW, DERIVED, RESULTS, THRESHOLD, MAX_STOP_DOSE

e = pd.read_csv(RAW / "eicu_vs.csv")
flow = [("First-unit vasopressor episodes >= 6 h", len(e))]
e = e[(e.last_unk == 0) & e.last_d.notna()]
flow.append(("Dose at discontinuation computable", len(e)))
e = e[e.max_d >= 0.05]
flow.append(("Peak NEE >= 0.05", len(e)))
e = e[e.cmo_by_24h == 0]
flow.append(("No care limitation for comfort <= 24 h after discontinuation", len(e)))
e = e[e.last_d <= MAX_STOP_DOSE].reset_index(drop=True)
flow.append(("NEE at discontinuation <= 0.2 (analytic cohort)", len(e)))

# arterial MAP when available, otherwise non-invasive
e["map_m2"] = e.map_m2_art.fillna(e.map_m3_nibp)
e["map_24"] = e.map_24_art.fillna(e.map_24_nibp)
e["ut"] = e.unittype.astype("category").cat.codes
e["abrupt"] = (e.last_d > THRESHOLD).astype(int)
e["d24nr"] = ((e.hosp_death_h <= 24) & ~(e.restart_gap_h <= 24)).astype(int)
flow += [("  tapered (NEE <= 0.05)", int((e.abrupt == 0).sum())),
         ("  abrupt (NEE > 0.05)", int(e.abrupt.sum())),
         ("  hospitals", e.hospitalid.nunique())]

e.to_csv(DERIVED / "eicu_analytic.csv", index=False, float_format="%.17g")
pd.DataFrame(flow, columns=["step", "n"]).to_csv(RESULTS / "flow_eicu.csv", index=False)
for s, n in flow:
    print(f"{s:60s} {n:>7,}")
