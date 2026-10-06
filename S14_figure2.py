#!/usr/bin/env python3
"""Figure 2: differential expression and functional enrichment.

The four panels carry no titles; the two volcano panels identify their disease
through a legend (HCC / HF) instead, so the same wording is used everywhere in
the manuscript. Outputs PDF and TIF only.
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RAW, FIG, ROOT

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'legend.handletextpad': 0.4,
})

BLUE = '#5B9BD5'; RED = '#C00000'; LGRAY = '#C0C0C0'; DGRAY = '#303030'


def panel_label(ax, s):
    ax.text(-0.22, 1.05, s, transform=ax.transAxes, fontsize=9, va='top')


def volcano(ax, deg, label):
    up = (deg.log2FC > 1) & (deg.padj < 0.05)
    down = (deg.log2FC < -1) & (deg.padj < 0.05)
    ns = ~(up | down)
    ax.scatter(deg.log2FC[ns], -np.log10(deg.padj[ns]), s=2, c=LGRAY, alpha=0.45,
               linewidths=0, rasterized=True, label='NS')
    ax.scatter(deg.log2FC[up], -np.log10(deg.padj[up]), s=5, c=RED, alpha=0.7,
               linewidths=0, rasterized=True, label='Up')
    ax.scatter(deg.log2FC[down], -np.log10(deg.padj[down]), s=5, c=BLUE, alpha=0.7,
               linewidths=0, rasterized=True, label='Down')
    ax.axhline(-np.log10(0.05), color='#999', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(1, color='#999', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(-1, color='#999', ls='--', lw=0.5, alpha=0.4)
    # the disease is carried by the x-axis label rather than by a title or legend
    ax.set_xlabel(f'log$_2$ FC in {label}')
    for _s in ('top', 'right'):
        ax.spines[_s].set_visible(False)      # no frame; the axes stay
    _h, _l = ax.get_legend_handles_labels()
    ax.legend([_h[1], _h[2], _h[0]], [_l[1], _l[2], _l[0]], frameon=False, fontsize=6,
              loc='upper right', ncol=1, handletextpad=0.3, labelspacing=0.25,
              borderaxespad=0.3)


fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))

volcano(axes[0, 0], pd.read_csv(f'{RAW}/DEG_HCC.csv'), 'HCC')
axes[0, 0].set_ylabel('-log$_{10}$ adjusted P')
panel_label(axes[0, 0], 'a')

volcano(axes[0, 1], pd.read_csv(f'{RAW}/DEG_CVD.csv'), 'HF')
axes[0, 1].set_ylabel('-log$_{10}$ adjusted P')
panel_label(axes[0, 1], 'b')

# the enrichment table the paper reports lives in results/, not next to the raw inputs:
# the raw/data copy is a stale artefact whose GO terms match neither the submitted
# figure nor the manuscript text
enr = json.load(open(os.path.join(ROOT, 'results', 'enrichment_results.json')))
# panel c: warm/purple hues only -- no blue or green, which belong to panel d
WARM = ['#D55E00', '#E69F00', '#CC79A7', '#8E7CC3', '#999999']
COOL = ['#0072B2', '#56B4E9', '#009E73', '#4E7A8C', '#7A9C6B']


def enrich_bar(ax, terms, palette):
    terms = sorted(terms, key=lambda x: -np.log10(x['Adjusted P-value']))[:5][::-1]
    names = [t['Term'].split(' (GO')[0] for t in terms]
    vals = [-np.log10(t['Adjusted P-value']) for t in terms]
    ax.barh(range(len(vals)), vals, color=palette[:len(vals)], height=0.6)
    ax.set_yticks(range(len(vals))); ax.set_yticklabels(names, fontsize=5.5)
    ax.set_xlabel('-log$_{10}$ adjusted P')
    for _s in ('top', 'right'):
        ax.spines[_s].set_visible(False)


if 'HCC_down' in enr:
    enrich_bar(axes[1, 0], enr['HCC_down'], WARM)
    panel_label(axes[1, 0], 'c')
if 'CVD_HF_up' in enr:
    enrich_bar(axes[1, 1], enr['CVD_HF_up'], COOL)
    panel_label(axes[1, 1], 'd')

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'Figure2_de_enrichment_v2.tif'),
                dpi=600, bbox_inches='tight')
plt.close(fig)
print('written', os.path.join(FIG, 'Figure2_de_enrichment_v2.tif'))
