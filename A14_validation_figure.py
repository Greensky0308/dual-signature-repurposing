#!/usr/bin/env python3
"""A14: revised Figure 4 (three-source prioritisation of the candidates).

Panel a  waterfall of the primary chain (dual-reversal -> touchstone ->
         documented indication), with the published values marked.
Panel b  candidate counts by mechanism class in the final node, annotated with
         the Fisher exact enrichment result (HDAC and tubulin inhibitors,
         which score identically).
Panel c  HCC and HF reversal scores of the cardiac-glycoside candidates: the
         cardenolides that reach the touchstone layer, plus digitoxin, which
         enters under the sensitivity chain (PubChem-CID requirement removed).

Counts come from A12_revised_chain.json; compound scores come from the
manuscript-aggregation table raw/data/drug_scores_all3.csv.
Outputs: figures/Figure4_validation.tif
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from config import RAW as DATA, OUT, FIG, HUB

chain = json.load(open(f'{OUT}/A12_revised_chain.json'))
P, S = chain['primary'], chain['sensitivity']
enr = {e['mechanism_class']: e for e in P['enrichment']}   # primary chain
scores = pd.read_csv(f'{DATA}/drug_scores_all3.csv').set_index('true_name')

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 8, 'pdf.fonttype': 42})
INK, GREY, BLUE, RED, GREEN = '#1a1a1a', '#8f8f8f', '#2c6fbb', '#c0392b', '#1e7a4c'

fig = plt.figure(figsize=(7.2, 5.2))

# ---------------- panel a: funnel ----------------------------------------
ax = fig.add_axes([0.02, 0.60, 0.96, 0.36])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')

stages = [('Dual-reversal screening', P['dual_cid']),
          ('Touchstone annotation', P['touchstone']),
          ('Documented indication', P['final'])]
top = max(v for _, v in stages)
funnel_blue = ['#b9d0e8', '#9dbddg' if False else '#9dbdd8', '#84aacd']
cx, maxw = 0.42, 0.33
band_h = 0.31
y0 = 0.97

widths = [max(w * (v / top), 0.075) for w, (_, v) in zip([maxw] * 3, stages)]
widths.append(widths[-1] * 0.72)          # slight closing taper at the bottom

for i, (name, v) in enumerate(stages):
    yt, yb = y0 - i * band_h, y0 - (i + 1) * band_h
    wt, wb = widths[i], widths[i + 1]
    ax.add_patch(plt.Polygon(
        [(cx - wt / 2, yt), (cx + wt / 2, yt), (cx + wb / 2, yb), (cx - wb / 2, yb)],
        closed=True, facecolor=funnel_blue[i], edgecolor='white', linewidth=1.4, zorder=2))
    yc = (yt + yb) / 2
    ax.text(cx - wt / 2 - 0.015, yc, name, ha='right', va='center',
            fontsize=8.2, color=INK)
    ax.text(cx + wt / 2 + 0.015, yc, f'{v:,}', ha='left', va='center',
            fontsize=9.6, color=INK, fontweight='bold')
    ax.text(cx + wt / 2 + 0.015, yc - 0.075, f'({100 * v / P["dual_cid"]:.0f}%)',
            ha='left', va='center', fontsize=7.0, color=GREY)

fig.text(0.030, 0.965, 'a', fontsize=9, color=INK, va='bottom')

# ---------------- panel b: mechanism classes -----------------------------
ax2 = fig.add_axes([0.115, 0.10, 0.36, 0.44])
cls = P['classes']
names = [c['mechanism_class'] for c in cls][::-1]
counts = [c['n'] for c in cls][::-1]
# tubulin scores identically to HDAC in the enrichment test, so both are highlighted
ENRICHED = {'HDAC inhibitor', 'Tubulin inhibitor'}
colors = [RED if n in ENRICHED else BLUE for n in names]
ax2.barh(np.arange(len(names)), counts, color=colors, alpha=0.85, zorder=3)
ax2.set_yticks(np.arange(len(names)))
SHORT = {'HDAC inhibitor': 'HDAC', 'Topoisomerase inhibitor': 'Topoisom.',
         'Tubulin inhibitor': 'Tubulin', 'CDK inhibitor': 'CDK',
         'Proteasome inhibitor': 'Proteasome', 'ATP1A1 target': 'ATP1A1'}
ax2.set_yticklabels([SHORT.get(n, n) for n in names], fontsize=7)
for i, v in enumerate(counts):
    ax2.text(v + 0.06, i, str(v), va='center', fontsize=6.8, color=INK)
ax2.set_xlabel('candidates', fontsize=7.2)
ax2.set_xlim(0, max(counts) * 1.35)
ax2.tick_params(labelsize=7)
for s in ('top', 'right'):
    ax2.spines[s].set_visible(False)
h = enr.get('HDAC inhibitor', {})
ax2.text(0.99, 1.03,
         f"HDAC & tubulin inhibitors: OR {h.get('odds_ratio', float('nan')):.1f}, "
         f"P = {h.get('p', float('nan')):.4f}",
         transform=ax2.transAxes, ha='right', va='bottom',
         fontsize=6.6, color=RED)
fig.text(0.030, 0.548, 'b', fontsize=9, color=INK, va='bottom')

# ---------------- panel c: cardiac glycoside scores ----------------------
ax3 = fig.add_axes([0.615, 0.10, 0.345, 0.44])
cg = ['peruvoside', 'sarmentogenin', 'strophanthidin', 'cinobufagin', 'digitoxin']
labels, hcc, hf = [], [], []
for d in cg:
    if d in scores.index:
        labels.append(d)
        hcc.append(float(scores.loc[d, 'hcc_score']))
        hf.append(float(scores.loc[d, 'cvd_score']))
y = np.arange(len(labels))
ax3.barh(y - 0.19, hcc, height=0.36, color=BLUE, alpha=0.85, label='HCC', zorder=3)
ax3.barh(y + 0.19, hf, height=0.36, color=GREEN, alpha=0.85, label='HF', zorder=3)
ax3.axvline(0, color=INK, lw=0.7)
ax3.set_yticks(y)
ax3.set_yticklabels(labels, fontsize=6.8)
ax3.set_xlabel('reversal score', fontsize=7.2)
ax3.tick_params(labelsize=7)
ax3.legend(fontsize=6.2, frameon=False, loc='upper right', ncol=1)
for s in ('top', 'right'):
    ax3.spines[s].set_visible(False)
fig.text(0.548, 0.548, 'c', fontsize=9, color=INK, va='bottom')


fig.savefig(f'{FIG}/Figure4_validation.tif', dpi=400)
print(f'wrote {FIG}/Figure4_validation.tif')
print(f'  funnel: {P["dual_cid"]} -> {P["touchstone"]} -> {P["final"]}')
print(f'  cardiac-glycoside scores: {dict(zip(labels, zip(hcc, hf)))}')
