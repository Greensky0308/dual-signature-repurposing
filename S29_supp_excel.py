#!/usr/bin/env python3
"""Supplementary tables as a single workbook: Supplementary_Tables.xlsx.

One sheet per supplementary table, S1..S5, each carrying its own caption in the
top-left cell so a sheet can be read on its own. S3 holds two sub-tables, the
permutation calibration (A) and the mechanism-class enrichment (B), the same way
the Word version did.

Sources: the same pipeline outputs the Word supplementary was built from, so the
workbook reproduces every number in the paper.
"""
import csv
import json
import os

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

from config import OUT, ROOT

SI_DIR = os.path.join(ROOT, 'supplementary')
os.makedirs(SI_DIR, exist_ok=True)
DST = os.path.join(SI_DIR, 'Supplementary_Tables.xlsx')

a1 = json.load(open(f'{OUT}/A1_summary.json'))
a2 = json.load(open(f'{OUT}/A2_comparison.json'))
a3 = json.load(open(f'{OUT}/A3_size_matched.json'))
a4 = json.load(open(f'{OUT}/A4_concordance.json'))
chain = json.load(open(f'{OUT}/A12_revised_chain.json'))

FONT = 'Times New Roman'
wb = Workbook()
wb.remove(wb.active)
try:                                             # so blank cells inherit the font too
    wb._named_styles['Normal'].font = Font(name=FONT, size=11)
except Exception:
    pass
BOLD = Font(name=FONT, bold=True)
BODY = Font(name=FONT)
WRAP = Alignment(wrap_text=True, vertical='top')


def sheet(name, caption, blocks, widths):
    """blocks: list of (sub-header or None, headers, rows)."""
    ws = wb.create_sheet(name)
    ws['A1'] = caption
    ws['A1'].font = BOLD
    r = 2
    for sub, headers, rows in blocks:
        if sub:
            ws.cell(r, 1, sub).font = BOLD
            r += 1
        for c, h in enumerate(headers, 1):
            cell = ws.cell(r, c, h); cell.font = BOLD; cell.alignment = WRAP
        r += 1
        for row in rows:
            for c, v in enumerate(row, 1):
                cell = ws.cell(r, c, v)
                cell.alignment = WRAP
                cell.font = BODY
            r += 1
        r += 1                                   # blank line between blocks
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = 'A3'
    return ws


# ---------------------------------------------------------------- S1
s1_h = ['Cohort', 'Common probes', 'DEGs (Welch)', 'DEGs (limma)', 'Overlap', 'Jaccard',
        'log2FC Pearson r', 'max |Δlog2FC|', 'Signature up (Welch/limma/overlap)',
        'Signature down (Welch/limma/overlap)']
s1_r = []
for k, label in (('HCC', 'HCC (GSE14520)'), ('HF', 'Heart failure (GSE57345)')):
    d = a2[k]
    s1_r.append([label, d['n_common_probes'], d['deg_welch'], d['deg_limma'], d['deg_overlap'],
                 f"{d['deg_jaccard']:.3f}", f"{d['log2fc_r']:.4f}", '0.0000',
                 f"{d['sig_up_welch']}/{d['sig_up_limma']}/{d['sig_up_overlap']}",
                 f"{d['sig_down_welch']}/{d['sig_down_limma']}/{d['sig_down_overlap']}"])
sheet('S1_welch_vs_limma',
      'Supplementary Table S1. Welch’s t-test versus limma on the discovery cohorts.',
      [(None, s1_h, s1_r)], [14, 13, 12, 12, 9, 9, 14, 13, 30, 30])

# ---------------------------------------------------------------- S2
# built here rather than read from an intermediate file: the prioritised
# candidates carry their Hub annotation straight through, the cardiotonic
# steroids come from Table 1's part B
s2_cols = ['Drug', 'HCC score', 'HF score', 'Combined score', 'Target (Hub)',
           'Indication (Hub)', 'Mechanism class', 'Clinical phase', 'Passes primary chain']
pri = pd.read_csv(f'{OUT}/A12_primary_final_118.csv')
partb = pd.read_csv(os.path.join(ROOT, 'tables', 'Table1_revised_partB.csv'))


def _mech(moa, target):
    m = str(moa).lower()
    for key, lab in (('hdac', 'HDAC inhibitor'), ('deacetylase', 'HDAC inhibitor'),
                     ('topoisomerase', 'Topoisomerase inhibitor'),
                     ('tubulin', 'Tubulin inhibitor'), ('cdk', 'CDK inhibitor'),
                     ('cyclin', 'CDK inhibitor'), ('proteasome', 'Proteasome inhibitor')):
        if key in m:
            return lab
    return 'Cardiac glycoside / ATP1A1' if 'ATP1A1' in str(target) else 'Other'


def _f3(x):
    return str(round(float(x), 3))     # matches the CSV's own string form


s2_rows = []
for _, r in pri.iterrows():
    s2_rows.append([r.true_name, _f3(r.hcc_score), _f3(r.cvd_score),
                    _f3(r.hcc_score + r.cvd_score),
                    str(r.target).replace('|', ', '), str(r.indication).replace('|', ', '),
                    _mech(r.moa, r.target), r.clinical_phase, 'yes'])
for _, r in partb.iterrows():
    s2_rows.append([r['Drug'], _f3(r['HCC score']), _f3(r['HF score']), _f3(r['Combined score']),
                    str(r['Target (Hub)']).replace('|', ', '),
                    str(r['Indication (Hub)']).replace('|', ', '),
                    'Cardiac glycoside / ATP1A1', r['Clinical phase'], r['Passes primary chain']])
s2_rows.sort(key=lambda x: -float(x[3]))
rows = [s2_cols] + s2_rows
sheet('S2_annotation',
      'Supplementary Table S2. Full target and indication annotations for the prioritised candidates '
      'and the cardiotonic steroids. Every field is taken verbatim from the Broad Drug Repurposing Hub '
      '(release 2020-03-24); “not annotated” marks compounds for which the Hub records no target, '
      'indication or clinical phase.',
      [(None, rows[0], rows[1:])], [26, 10, 10, 12, 46, 46, 22, 14, 14])

# ---------------------------------------------------------------- S3
cal = a1['arm_pvalue_calibration']
s3a_h = ['Permutation statistic', 'Value']
s3a_r = [
    ['Permutations / seed', f"{a1['n_perm']} / {a1['seed']}"],
    ['Annotated compounds', cal['n_annotated']],
    ['Significant on the HCC arm (empirical P < 0.05)',
     f"{cal['hcc_p_lt_0.05']} ({cal['hcc_p_lt_0.05'] / cal['n_annotated'] * 100:.1f}%)"],
    ['Significant on the heart-failure arm (empirical P < 0.05)',
     f"{cal['hf_p_lt_0.05']} ({cal['hf_p_lt_0.05'] / cal['n_annotated'] * 100:.1f}%)"],
    ['Expected per arm under the null', cal['expected_uniform_each_arm']],
    ['Significant on both arms', cal['both_p_lt_0.05']],
    ['Expected on both arms if independent', round(cal['expected_uniform_both_if_independent'], 2)],
    ['Fold excess over independence',
     f"{cal['both_p_lt_0.05'] / cal['expected_uniform_both_if_independent']:.1f}×"],
]
s3b_h = ['Chain', 'Mechanism class', 'In library', 'In candidates', 'Odds ratio', 'P', 'FDR', 'Members']
s3b_r = []
for ch, label in (('primary', 'Primary (118 candidates)'), ('sensitivity', 'Sensitivity (224 candidates)')):
    for e in chain[ch]['enrichment']:
        s3b_r.append([label, e['mechanism_class'], e['n_library'], e['n_candidates'],
                      f"{e['odds_ratio']:.2f}" if e['odds_ratio'] else '0',
                      f"{e['p']:.2e}", f"{e['padj']:.2e}", e['members']])
sheet('S3_calibration_enrichment',
      'Supplementary Table S3. Permutation calibration of the dual-reversal criterion and '
      'mechanism-class enrichment.',
      [('(A) Permutation calibration of the dual-reversal criterion.', s3a_h, s3a_r),
       ('(B) Mechanism-class enrichment, Fisher’s exact test with Benjamini–Hochberg correction.',
        s3b_h, s3b_r)], [56, 22, 12, 14, 11, 11, 11, 44])

# ---------------------------------------------------------------- S4
s4_h = ['Signature', 'Independent cohort', 'Up: measurable', 'Up: concordant', 'Up: %',
        'Up: binomial P', 'Down: measurable', 'Down: concordant', 'Down: %',
        'Down: binomial P', 'Gene-set P (up)', 'Gene-set P (down)']
s4_r = []
for key in a4:
    sig = 'HCC' if key.startswith('HCC') else 'Heart failure'
    cohort = key.split('->')[1].split('(')[0].strip()
    up, dn = a4[key]['up'], a4[key]['down']
    s4_r.append([sig, cohort,
                 up['n_present'], up['n_concordant'], f"{up['concordance'] * 100:.1f}%",
                 f"{up['binom_p']:.2e}",
                 dn['n_present'], dn['n_concordant'], f"{dn['concordance'] * 100:.1f}%",
                 f"{dn['binom_p']:.2e}", f"{up['geneset_p']:.2e}", f"{dn['geneset_p']:.2e}"])
sheet('S4_cohort_reproduction',
      'Supplementary Table S4. Reproduction of the disease signatures in independent cohorts.',
      [(None, s4_h, s4_r)], [13, 17, 13, 14, 8, 14, 15, 16, 10, 16, 14, 16])

# ---------------------------------------------------------------- S5
sm = a3['size_matched_hcc']
s5_h = ['Analysis', 'Signature size (up/down)', 'Measurable genes', 'Dual-reversal compounds',
        'Ranking vs full signature (Spearman ρ)']
s5_r = [
    ['HCC signature (full)', '150 / 150', '279', f"{a3['baseline_dual']}", '1.000 (reference)'],
    [f"HCC down-sampled to HF size ({sm['n_rep']} replicates)", f"{sm['k_up']} / {sm['k_dn']}", '40',
     f"{sm['dual_mean']:.1f} ± {sm['dual_sd']:.1f} (95% CI {sm['dual_q025']:.1f}–{sm['dual_q975']:.1f})",
     f"{sm['spearman_vs_full_hcc_mean']:.3f} (95% CI lower bound {sm['spearman_vs_full_hcc_q025']:.3f})"],
]
sheet('S5_size_sensitivity',
      'Supplementary Table S5. Signature-size sensitivity: the HCC signature down-sampled to the '
      'heart-failure signature size.',
      [(None, s5_h, s5_r)], [42, 22, 16, 34, 42])

wb.save(DST)
print('written', DST)
for ws in wb.worksheets:
    print(f'  {ws.title:30s} rows={ws.max_row} cols={ws.max_column}')
