"""Crude descriptive results quoted in the Results section.

Output: results/descriptives.csv
"""
import numpy as np
import pandas as pd
from common import RESULTS, load_mimic, load_eicu

rows = []
def add(db, item, abrupt, tapered=None):
    rows.append(dict(database=db, item=item, abrupt=abrupt, tapered=tapered))
    print(f"{db:9s} {item:55s} {abrupt}" + ("" if tapered is None else f" | {tapered}"))

d, e = load_mimic(), load_eicu()
for db, D in (("MIMIC-IV", d), ("eICU-CRD", e)):
    a, tp = D[D.abrupt == 1], D[D.abrupt == 0]
    add(db, "n", len(a), len(tp))
    add(db, "abrupt stops, %", round(100 * D.abrupt.mean(), 1))
    add(db, "median NEE at discontinuation", round(a.last_d.median(), 3), round(tp.last_d.median(), 3))
    add(db, "vasopressin at discontinuation, %", round(100 * a.vaso_late.mean(), 1), round(100 * tp.vaso_late.mean(), 1))
    add(db, "crude 24-h re-initiation, %", round(100 * (a.restart_gap_h <= 24).mean(), 1),
        round(100 * (tp.restart_gap_h <= 24).mean(), 1))

# timing of re-initiation (MIMIC-IV)
gap = d.restart_gap_h[d.restart_gap_h <= 24]
add("MIMIC-IV", "median h to re-initiation (<=24 h)", round(gap.median(), 1))
add("MIMIC-IV", "IQR h to re-initiation", f"{gap.quantile(.25):.1f}-{gap.quantile(.75):.1f}")
add("MIMIC-IV", "re-initiations within 6 h, % of those <=24 h", round(100 * (gap <= 6).mean(), 0))

# state at 72 h (crude, mutually exclusive): re-initiation within 72 h;
# death within 72 h without re-initiation within 72 h; neither
r = d.restart_gap_h.fillna(np.inf); de = d.death_h.fillna(np.inf)
re72 = r <= 72
dd72 = (de <= 72) & (r > 72)
for lab, s in (("re-initiation by 72 h, %", re72), ("death without re-initiation by 72 h, %", dd72)):
    add("MIMIC-IV", lab, round(100 * s[d.abrupt == 1].mean(), 1), round(100 * s[d.abrupt == 0].mean(), 1))

# hospital-level consistency (eICU-CRD): hospitals with >=40 patients and >=8 per group
h = e.groupby("hospitalid").filter(lambda s: len(s) >= 40 and s.abrupt.sum() >= 8 and (1 - s.abrupt).sum() >= 8)
pos = sum((s[s.abrupt == 1].restart_gap_h <= 24).mean() > (s[s.abrupt == 0].restart_gap_h <= 24).mean()
          for _, s in h.groupby("hospitalid"))
add("eICU-CRD", "hospitals evaluable", h.hospitalid.nunique())
add("eICU-CRD", "hospitals with higher re-initiation after abrupt stops", pos)

pd.DataFrame(rows).to_csv(RESULTS / "descriptives.csv", index=False)
