#!/usr/bin/env python3
"""Figure 3: dual-signature connectivity mapping.

Panel a is the HCC-versus-heart-failure reversal-score plane; panels b and c
summarise the prioritised candidates. Every target and mechanism label in b and c
now comes from the Broad Drug Repurposing Hub (release 2020-03-24) rather than
from a hand-assembled annotation table, so the figure and the manuscript use the
same annotation source. Outputs PDF and TIF only.
"""
import os
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from config import OUT, RAW, FIG

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

BLUE = '#5B9BD5'; ORANGE = '#D08A6A'; PURPLE = '#8E7CC3'; RED = '#C00000'

drugs = pd.read_csv(f'{RAW}/drug_scores_all3.csv')
final = pd.read_csv(f'{OUT}/A12_primary_final_118.csv')
dual = drugs[(drugs.hcc_score > 0) & (drugs.cvd_score > 0) & (drugs.has_name)
             & (drugs.pubchem.astype(str) != '-666')]   # the 867 dual-reversal node
final = final.assign(combined=final.hcc_score + final.cvd_score).sort_values(
    'combined', ascending=False)


def mclass(moa):
    m = str(moa).lower()
    if 'hdac' in m:
        return ORANGE
    if 'cdk' in m:
        return BLUE
    return PURPLE


fig = plt.figure(figsize=(7.8, 3.0))
# manual placement: a wide, generous gap a->b (panel b's labels sit in it), tight b->c
ax_a = fig.add_axes([0.090, 0.17, 0.240, 0.68])
ax_b = fig.add_axes([0.545, 0.17, 0.205, 0.68])
ax_c = fig.add_axes([0.845, 0.17, 0.130, 0.68])
axes = [ax_a, ax_b, ax_c]


fig.canvas.draw()      # needed before measuring tick labels


def panel_label(ax, s):
    """Panel letter at the top-left of the panel, aligned with the left edge of
    its tick labels (not with the axes spine)."""
    xs = [l.get_window_extent().x0 for l in ax.get_yticklabels() if l.get_text()]
    if ax.get_ylabel():
        xs.append(ax.yaxis.label.get_window_extent().x0)   # include the axis label
    if xs:
        x0 = fig.transFigure.inverted().transform((min(xs), 0))[0]
    else:
        x0 = ax.get_position().x0
    fig.text(x0, ax.get_position().y1 + 0.05, s, fontsize=9, va='bottom', ha='left')


# ------------------------------------------------- a: the score plane
ax = axes[0]
ax.scatter(drugs.hcc_score, drugs.cvd_score, s=2.5, c='#C0C0C0', alpha=0.5,
           linewidths=0, rasterized=True)
ax.scatter(dual.hcc_score, dual.cvd_score, s=12, c=RED, alpha=0.7,
           linewidths=0, rasterized=True)
# per-compound label offsets, as in the submitted figure, so labels do not collide
_cfg = {'actinomycin-d': ('left', 'bottom', 0.06, 0.06), 'plumbagin': ('right', 'bottom', -0.06, 0.06),
        'ER-27319': ('left', 'bottom', 0.06, 0.10), 'MW-STK33-23': ('left', 'top', 0.06, -0.10),
        'JNJ-7706621': ('right', 'bottom', -0.06, 0.06), 'thiomersal': ('right', 'top', -0.05, -0.07)}
for _, r in dual.sort_values('joint_cvd', ascending=False).head(6).iterrows():
    ha, va, dx, dy = _cfg.get(r.true_name, ('left', 'bottom', 0.06, 0.06))
    ax.text(r.hcc_score + dx, r.cvd_score + dy, r.true_name, fontsize=6, ha=ha, va=va)
ax.axhline(0, color='#999', ls='--', lw=0.5, alpha=0.4, zorder=0)
ax.axvline(0, color='#999', ls='--', lw=0.5, alpha=0.4, zorder=0)
ax.set_xlabel('HCC reversal score'); ax.set_ylabel('HF reversal score')
for _s in ('top', 'right', 'left', 'bottom'):
    ax.spines[_s].set_visible(False)      # no frame around the scatter
panel_label(ax, 'a')

# ------------------------------------------- b: top candidates (Hub moa)
ax = axes[1]
top15 = final.head(15)
sns.heatmap(top15[['hcc_score', 'cvd_score']].values, annot=True, fmt='.2f',
            cmap='Greens', square=False, cbar=False,
            xticklabels=['HCC', 'HF'], yticklabels=top15.true_name.tolist(),
            linewidths=0.5, linecolor='white', ax=ax)
ax.set_yticklabels(top15.true_name.tolist(), fontsize=6, rotation=0)
_n = len(top15)
fig.canvas.draw()
_labs = [l for l in ax.get_yticklabels() if l.get_text()]
_x0fig = fig.transFigure.inverted().transform((min(l.get_window_extent().x0 for l in _labs), 0))[0]
_pos = ax.get_position()
_ov = fig.add_axes([0, 0, 1, 1], zorder=5)          # overlay in figure fractions
_ov.set_xlim(0, 1); _ov.set_ylim(0, 1); _ov.axis('off')
for i, moa in enumerate(top15.moa):
    _ov.plot([_x0fig - 0.016], [_pos.y0 + _pos.height * (1 - (i + 0.5) / _n)],  # row centre
             marker='o', ms=4.5, color=mclass(moa), clip_on=False, mec='none')
panel_label(ax, 'b')

# ------------------------------------- c: target families (Hub annotation)
ax = axes[2]
import re as _re
_fams = []
for _t in final.head(25).target:
    for _x in str(_t).split('|'):
        _x = _x.strip()
        if _x and _x.lower() not in ('nan', 'none'):
            _fams.append(_re.sub(r'\d+$', '', _x))     # HDAC1 -> HDAC, NR3C1 -> NR3C
_top = Counter(_fams).most_common(8)[::-1]
ax.barh(range(len(_top)), [v for _, v in _top], color=BLUE, height=0.6)
ax.set_yticks(range(len(_top)))
ax.set_yticklabels([k for k, _ in _top], fontsize=6)
ax.set_xlabel('Targets per family')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
panel_label(ax, 'c')

# panels are placed manually, so the canvas is used as-is (no tight crop)
fig.savefig(os.path.join(FIG, 'Figure3_connectivity.tif'), dpi=600)
plt.close(fig)
print('written', os.path.join(FIG, 'Figure3_connectivity.tif'))
