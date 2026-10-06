#!/usr/bin/env python3
"""A13: the filtering-progression panel (Figure 1a).

GEO datasets -> differential expression -> L1000 screening -> dual-reversal
candidates -> touchstone compounds -> documented indication, drawn as a vertical
stack of rounded boxes with arrows. Counts come from A12_revised_chain.json.

The same drawing backs panel a of Figure 1 (S12 imports `draw_workflow`), so the
published workflow and this standalone panel can never drift apart.

Outputs: figures/Figure_workflow.tif
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from config import OUT, FIG

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'pdf.fonttype': 42})

BLUE1, BLUE2, CREAM = '#d9e4ef', '#ebf1f7', '#faecd3'
SALMON, MAUVE, PINK = '#f5dcd9', '#e6daed', '#fbe3e0'
DARKRED = '#8c2f22'


def draw_workflow(ax, P, label_letter=True):
    """Draw the progression on `ax`. `P` holds the chain counts."""
    boxes = [
        ('HCC + heart failure datasets', 'GSE14520  -  GSE57345', BLUE2, 'black'),
        ('Differential expression', 'disease signatures', BLUE1, 'black'),
        ('L1000 connectivity map', f"{473647:,} perturbation profiles", CREAM, 'black'),
        ('Dual-reversal screening', f"{P['dual_cid']} candidates", SALMON, 'black'),
        ('Touchstone annotation', f"{P['touchstone']} compounds", MAUVE, 'black'),
        ('Documented indication', f"{P['final']} candidates", BLUE1, 'black'),
        ('HDAC inhibitors + cardiotonic steroids', '', PINK, DARKRED),
    ]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    if label_letter:
        ax.text(1, 99, 'a', fontsize=9, color='black', va='top')

    top, gap, h = 94.0, 2.6, 10.5
    for i, (title, detail, fill, tcol) in enumerate(boxes):
        y = top - i * (h + gap) - h
        ax.add_patch(FancyBboxPatch((17, y), 66, h,
                                    boxstyle='round,pad=0.4,rounding_size=1.2',
                                    linewidth=0, facecolor=fill, zorder=2))
        if detail:
            ax.text(50, y + h * 0.62, title, ha='center', va='center', fontsize=7.0,
                    color=tcol, zorder=3)
            ax.text(50, y + h * 0.26, detail, ha='center', va='center', fontsize=6.5,
                    color=tcol, zorder=3)
        else:
            ax.text(50, y + h * 0.5, title, ha='center', va='center', fontsize=7.0,
                    color=tcol, zorder=3)
        if i < len(boxes) - 1:
            ax.add_patch(FancyArrowPatch((50, y - 0.6), (50, y - gap + 0.6),
                                         arrowstyle='-|>', mutation_scale=11,
                                         linewidth=1.1, color='black', zorder=1))


if __name__ == '__main__':
    chain = json.load(open(f'{OUT}/A12_revised_chain.json'))['primary']
    fig = plt.figure(figsize=(3.9, 7.2))
    draw_workflow(fig.add_axes([0.02, 0.02, 0.96, 0.96]), chain)
    fig.savefig(f'{FIG}/Figure_workflow.tif', dpi=400)
    print(f'wrote {FIG}/Figure_workflow.tif')
    print(f"nodes: 473,647 -> {chain['dual_cid']} -> {chain['touchstone']} -> {chain['final']}")
