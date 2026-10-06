#!/usr/bin/env python3
"""Supplementary Figure 1: generalisability to coronary artery disease.

Same figure as submitted, with the x-axis label shortened to "log2 FC" so that
it matches Figure 2. Outputs PDF and TIF only.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RAW, FIG

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

BLUE = '#5B9BD5'; RED = '#C00000'; LGRAY = '#C0C0C0'; ORANGE = '#D08A6A'


def panel_label(ax, s):
    ax.text(-0.22, 1.05, s, transform=ax.transAxes, fontsize=9, va='top')


fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))


def volcano(ax, deg, title, n_sig):
    up = (deg.log2FC > 1) & (deg.padj < 0.05)
    down = (deg.log2FC < -1) & (deg.padj < 0.05)
    ns = ~(up | down)
    ax.scatter(deg.log2FC[ns], -np.log10(deg.padj[ns]), s=2, c=LGRAY, alpha=0.45,
               linewidths=0, rasterized=True)
    ax.scatter(deg.log2FC[up], -np.log10(deg.padj[up]), s=5, c=RED, alpha=0.7,
               linewidths=0, rasterized=True)
    ax.scatter(deg.log2FC[down], -np.log10(deg.padj[down]), s=5, c=BLUE, alpha=0.7,
               linewidths=0, rasterized=True)
    ax.axhline(-np.log10(0.05), color='#999', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(1, color='#999', ls='--', lw=0.5, alpha=0.4)
    ax.axvline(-1, color='#999', ls='--', lw=0.5, alpha=0.4)
    ax.set_xlabel('log$_2$ FC'); ax.set_ylabel('-log$_{10}$ adjusted P')
    ax.set_title(f'{title} (n={n_sig})', fontsize=8)
    for _s in ('top', 'right'):
        ax.spines[_s].set_visible(False)      # no frame; the axes stay


volcano(axes[0], pd.read_csv(f'{RAW}/DEG_CAD.csv'), 'CAD', 2434)
panel_label(axes[0], 'a')

cad_cand = pd.read_csv(f'{RAW}/candidate_drugs_HCC+CAD.csv')
ax = axes[1]
ax.scatter(cad_cand.hcc_score, cad_cand.cad_score, s=8, c=ORANGE, alpha=0.6,
           linewidths=0, rasterized=True)
ax.axhline(0, color='#999', ls='--', lw=0.5, alpha=0.4)
ax.axvline(0, color='#999', ls='--', lw=0.5, alpha=0.4)
ax.set_xlabel('HCC reversal score'); ax.set_ylabel('CAD reversal score')
ax.set_title('HCC + CAD dual-reversal', fontsize=8)
for _s in ('top', 'right'):
    ax.spines[_s].set_visible(False)
panel_label(ax, 'b')

fig.savefig(os.path.join(FIG, 'FigureS1_cad.tif'), dpi=600, bbox_inches='tight')
plt.close(fig)
print('written', os.path.join(FIG, 'FigureS1_cad.tif'))
