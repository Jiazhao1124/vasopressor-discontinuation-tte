"""Overlap-weighted Aalen-Johansen cumulative incidence of re-initiation, death as competing event.

Output: results/cumulative_incidence_mimic.csv
"""
import numpy as np
import pandas as pd
from common import CONF_MIMIC, RESULTS, load_mimic, psfit, cboot

rng = np.random.default_rng(2026)
d = load_mimic()
g = d.subject_id.values
t, ps = psfit(d, 0.05, g, CONF_MIMIC, "cu")
w = np.where(t == 1, 1 - ps, ps)

HORIZON = 72.0
r = np.nan_to_num(d.restart_gap_h.values, nan=np.inf)
de = np.nan_to_num(d.death_h.values, nan=np.inf)
first = np.minimum(r, de)
C = np.zeros(len(d), int)                    # 0 none, 1 re-initiation, 2 death
C[(first <= HORIZON) & (r <= de)] = 1
C[(first <= HORIZON) & (de < r)] = 2
T = np.where(C > 0, first, HORIZON)
times = [6, 12, 24, 48, 72]


def aj(T, C, w):
    out, S, c1, c2 = {}, 1.0, 0.0, 0.0
    for tt in sorted(set(times) | set(np.unique(T[C > 0]))):
        m = T == tt
        n = w[T >= tt].sum()
        d1, d2 = w[m & (C == 1)].sum(), w[m & (C == 2)].sum()
        if n > 0:
            c1 += S * d1 / n
            c2 += S * d2 / n
            S *= 1 - (d1 + d2) / n
        if tt in times:
            out[tt] = (c1, c2)
    return out


def both(i):
    a, b = t[i] == 1, t[i] == 0
    return aj(T[i][a], C[i][a], w[i][a]), aj(T[i][b], C[i][b], w[i][b])


a, b = both(np.arange(len(d)))
bs = [both(i) for i in cboot(g, rng, 200)]
rows = []
for tt in times:
    df = [x[0][tt][0] - x[1][tt][0] for x in bs]
    rows.append(dict(hour=tt, cif_reinit_abrupt=a[tt][0], cif_reinit_tapered=b[tt][0],
                     diff=a[tt][0] - b[tt][0], diff_lo=np.percentile(df, 2.5),
                     diff_hi=np.percentile(df, 97.5),
                     cif_death_abrupt=a[tt][1], cif_death_tapered=b[tt][1]))
    print(f"{tt:>2} h  {100*a[tt][0]:.1f} vs {100*b[tt][0]:.1f}  diff {100*rows[-1]['diff']:.1f} "
          f"({100*rows[-1]['diff_lo']:.1f} to {100*rows[-1]['diff_hi']:.1f})")
pd.DataFrame(rows).to_csv(RESULTS / "cumulative_incidence_mimic.csv", index=False)
