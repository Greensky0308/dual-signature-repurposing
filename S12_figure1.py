#!/usr/bin/env python3
"""Figure 1: the complete filtering pathway, with the cohort and signature panels.

Panel a carries the whole chain — datasets, differential expression, the L1000
screen, dual-reversal screening (867), touchstone annotation (384) and documented
indication (118) — which is what Reviewer 1 asked for in comment 8. Panels b
(cohort composition) and c (signature size) are unchanged from the submitted
figure. Panel a is drawn by A13_workflow_figure.draw_workflow, so the whole
figure is reproducible from code.

Outputs a single TIF.
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from config import FIG, OUT
from A13_workflow_figure import draw_workflow

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

BLUE = '#5B9BD5'; ORANGE = '#D08A6A'; RED = '#C00000'; DGRAY = '#303030'

fig = plt.figure(figsize=(7.2, 4.4))
# a spans the full left column; b and c are stacked on the right so the bars stay
# short rather than being stretched to the height of the funnel
gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1], height_ratios=[1, 1],
                      wspace=0.38, hspace=0.45)

# ------------------------------------------------------------ a: the workflow
chain = json.load(open(os.path.join(OUT, 'A12_revised_chain.json')))['primary']
ax = fig.add_subplot(gs[:, 0])
draw_workflow(ax, chain, label_letter=False)

# --------------------------------------------------- b: cohort composition
ax = fig.add_subplot(gs[0, 1])
categories = ['HCC', 'Heart failure']
disease = [225, 177]; control = [220, 136]
x = np.arange(2)
_ymax_b = max(np.array(disease) + np.array(control)) * 1.18
ax.set_ylim(0, _ymax_b)
ax.bar(x, disease, 0.32, label='Disease', color=ORANGE)
ax.bar(x, control, 0.32, bottom=disease, label='Control', color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(categories, fontsize=6)
ax.tick_params(axis='y', labelsize=6)
ax.set_ylabel('Samples', fontsize=6)
ax.legend(frameon=False, fontsize=5.5, loc='upper right', ncol=1, labelspacing=0.25, handletextpad=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)


def seg_labels(ax, i, pairs, ymax):
    """Label each stacked segment; a segment too short to hold its label is put
    just outside the bar so the text never overruns the bar."""
    base = 0.0
    for val, _c in pairs:
        ycen = base + val / 2
        if val < 0.03 * ymax:
            ax.text(i + 0.30, ycen, str(val), ha='left', va='center', fontsize=5, color='black')
        else:
            ax.text(i, ycen, str(val), ha='center', va='center',
                    fontsize=5, color='white', fontweight='bold')
        base += val


seg_labels(ax, 0, [(disease[0], ORANGE), (control[0], BLUE)], _ymax_b)
seg_labels(ax, 1, [(disease[1], ORANGE), (control[1], BLUE)], _ymax_b)

# ------------------------------------------------------- c: signature size
ax = fig.add_subplot(gs[1, 1])
up = [150, 25]; down = [150, 23]
_ymax_c = max(np.array(up) + np.array(down)) * 1.18
ax.set_ylim(0, _ymax_c)
ax.bar(x, up, 0.32, label='Up-regulated', color=RED)
ax.bar(x, down, 0.32, bottom=up, label='Down-regulated', color=BLUE)
ax.set_xticks(x); ax.set_xticklabels(categories, fontsize=6)
ax.tick_params(axis='y', labelsize=6)
ax.set_ylabel('Genes', fontsize=6)
ax.legend(frameon=False, fontsize=5.5, loc='upper right', ncol=1, labelspacing=0.25, handletextpad=0.3)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
seg_labels(ax, 0, [(up[0], RED), (down[0], BLUE)], _ymax_c)
seg_labels(ax, 1, [(up[1], RED), (down[1], BLUE)], _ymax_c)

# panel letters: a and b share one y (top of the upper row), c sits at its own top
_ax_a, _ax_b, _ax_c = fig.axes[0], fig.axes[1], fig.axes[2]
_y_top = _ax_b.get_position().y1 + 0.045
_ap, _bp, _cp = _ax_a.get_position(), _ax_b.get_position(), _ax_c.get_position()
fig.text(_ap.x0 - 0.022, max(_ap.y1, _y_top) - 0.005, 'a', fontsize=9, va='top')
fig.text(_bp.x0 - 0.085, _y_top, 'b', fontsize=9, va='top')
fig.text(_cp.x0 - 0.085, _cp.y1 + 0.045, 'c', fontsize=9, va='top')

fig.savefig(os.path.join(FIG, 'Figure1_overview.tif'), dpi=600, bbox_inches='tight')
plt.close(fig)
print('written', os.path.join(FIG, 'Figure1_overview.tif'))
