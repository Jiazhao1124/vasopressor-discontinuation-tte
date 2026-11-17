"""Figure S1: cohort flow, MIMIC-IV and eICU-CRD (matplotlib)."""
from fig_common import *
import pandas as pd
from pathlib import Path
RES = Path(__file__).resolve().parents[1] / 'results'
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

W, H = 174, 135
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis('off')
K = 'black'; LW = 0.5; FS = 7; LH = 3.4   # line height (mm)


def box(x, y, w, h, lines, bold=False, align='center'):
    ax.add_patch(Rectangle((x, y), w, h, fc='white', ec=K, lw=LW))
    n = len(lines)
    for i, t in enumerate(lines):
        tx = x + w / 2 if align == 'center' else x + 1.6
        ax.text(tx, y + h / 2 + (i - (n - 1) / 2) * LH, t, ha=align, va='center', fontsize=FS,
                fontweight='bold' if (bold and i == 0) else 'normal')


def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', lw=LW, color=K, mutation_scale=6, shrinkA=0, shrinkB=0))


def hline(x1, x2, y):
    ax.plot([x1, x2], [y, y], color=K, lw=LW)


# Counts are read from the flow tables written by analysis/01 and 02.
fm_ = pd.read_csv(RES / "flow_mimic.csv").n.tolist()
fe_ = pd.read_csv(RES / "flow_eicu.csv").n.tolist()
n = lambda x: f"n = {x:,}"
D = {
    'MIMIC-IV (derivation)': dict(
        steps=[['ICU stays with a first', 'vasopressor episode ≥6 h', n(fm_[0])],
               [n(fm_[1])], [n(fm_[2])],
               ['Analytic cohort', f'{fm_[3]:,} stays', f'({fm_[6]:,} patients)']],
        excl=[['Peak NEE <0.05 µg/kg/min or', 'comfort care ≤2 h after stop', n(fm_[0] - fm_[1])],
              ['Comfort care or hospice in', 'discharge summary', n(fm_[1] - fm_[2])],
              ['NEE at stop >0.2 µg/kg/min', n(fm_[2] - fm_[3])]],
        arms=[['Tapered stop', 'NEE ≤0.05', n(fm_[4])], ['Abrupt stop', 'NEE >0.05', n(fm_[5])]]),
    'eICU-CRD (replication)': dict(
        steps=[['First-unit vasopressor', 'episodes ≥6 h', n(fe_[0])],
               [n(fe_[1])], [n(fe_[3])],
               ['Analytic cohort', f'{fe_[4]:,} patients', f'({fe_[7]} hospitals)']],
        excl=[['Dose at stop not computable', n(fe_[0] - fe_[1])],
              ['Peak NEE <0.05 µg/kg/min', f'(n = {fe_[1] - fe_[2]:,}) or comfort care', f'≤24 h (n = {fe_[2] - fe_[3]:,})'],
              ['NEE at stop >0.2 µg/kg/min', n(fe_[3] - fe_[4])]],
        arms=[['Tapered stop', 'NEE ≤0.05', n(fe_[5])], ['Abrupt stop', 'NEE >0.05', n(fe_[6])]])}

x0s = [3, 89]; bw = 36; ew = 41; gap = 4
ys = [10, 40, 66, 92]; hs = [14, 7, 7, 14]
for x0, (title, d) in zip(x0s, D.items()):
    ax.text(x0, 5, title, ha='left', va='center', fontsize=7.5, fontweight='bold')
    cx = x0 + bw / 2
    for i, (y, h, lines) in enumerate(zip(ys, hs, d['steps'])):
        box(x0, y, bw, h, lines, bold=(i == 3))
        if i < 3:
            arrow(cx, y + h, cx, ys[i + 1])
            ey = (y + h + ys[i + 1]) / 2
            eh = LH * len(d['excl'][i]) + 2.4
            box(x0 + bw + gap, ey - eh / 2, ew, eh, d['excl'][i], align='left')
            hline(cx, x0 + bw + gap, ey)
    # split into strategies
    ya = ys[3] + hs[3]; yb = ya + 6; aw = 30
    acx = [cx - aw / 2, cx - aw / 2 + 44]
    ax.plot([cx, cx], [ya, yb], color=K, lw=LW)
    hline(acx[0] + aw / 2, acx[1] + aw / 2, yb)
    for ax_, lines in zip(acx, d['arms']):
        arrow(ax_ + aw / 2, yb, ax_ + aw / 2, yb + 6)
        box(ax_, yb + 6, aw, 13, lines)

save(fig, 'output/FigureS1')
