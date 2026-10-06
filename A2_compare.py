#!/usr/bin/env python3
"""A2 step 3: compare the limma results against the manuscript's Welch's t-test
results, at probe level and at disease-signature level.

Then rebuild the 150/150 and 25/23 signatures from limma and score them against
L1000, to test whether the reported dual-reversal candidates are robust to the
differential-expression method.
"""
import json
import numpy as np
import pandas as pd

from config import RAW as DATA, OUT

JOBS = [('HCC', 'GSE14520', 'DEG_HCC.csv', 'GPL3921_probe2gene.csv'),
        ('HF', 'GSE57345', 'DEG_CVD.csv', 'GPL11532_probe2gene.csv')]


def gene_signature(deg, p2g, n):
    deg = deg.copy()
    deg['probe'] = deg['probe'].astype(str).str.strip('"')
    p2g = p2g.copy()
    p2g['ID'] = p2g['ID'].astype(str)
    p2g = p2g.dropna(subset=['symbol'])
    p2g = p2g[p2g['symbol'] != '']
    m = deg.merge(p2g[['ID', 'symbol']], left_on='probe', right_on='ID', how='inner')
    m['absfc'] = m.log2FC.abs()
    m = m.sort_values('absfc', ascending=False).drop_duplicates('symbol')
    up = m[(m.log2FC > 1) & (m.padj < 0.05)].sort_values('log2FC', ascending=False)
    dn = m[(m.log2FC < -1) & (m.padj < 0.05)].sort_values('log2FC', ascending=True)
    return up['symbol'].tolist()[:n], dn['symbol'].tolist()[:n]


report = {}
old_sigs = json.load(open(f'{DATA}/disease_signatures.json'))
keymap = {'HCC': 'HCC', 'HF': 'CVD_HF'}

for tag, geo, welch_file, p2g_file in JOBS:
    welch = pd.read_csv(f'{DATA}/{welch_file}')
    welch['probe'] = welch['probe'].astype(str).str.strip('"')
    lim = pd.read_csv(f'{OUT}/A2_limma_{geo}.csv')
    lim['probe'] = lim['probe'].astype(str).str.strip('"')
    m = welch.merge(lim, on='probe', suffixes=('_w', '_l'))
    print(f'\n===== {tag} ({geo}) =====')
    print(f'common probes: {len(m)}')

    deg_w = set(m.loc[(m.log2FC_w.abs() > 1) & (m.padj_w < 0.05), 'probe'])
    deg_l = set(m.loc[(m.log2FC_l.abs() > 1) & (m.padj_l < 0.05), 'probe'])
    inter = deg_w & deg_l
    print(f'DEGs (Welch) {len(deg_w)}  DEGs (limma) {len(deg_l)}  overlap {len(inter)}  '
          f'{100*len(inter)/max(len(deg_w),1):.1f}% of Welch   '
          f'{100*len(inter)/max(len(deg_l),1):.1f}% of limma   Jaccard {len(inter)/len(deg_w|deg_l):.3f}')
    conc = np.mean(np.sign(m.log2FC_w) == np.sign(m.log2FC_l))
    print(f'log2FC correlation r = {m.log2FC_w.corr(m.log2FC_l):.4f}   '
          f'sign agreement = {100*conc:.2f}%   '
          f'max absolute difference = {(m.log2FC_w - m.log2FC_l).abs().max():.4f}')
    print(f'padj correlation (Spearman) = '
          f'{pd.Series(m.padj_w).rank().corr(pd.Series(m.padj_l).rank()):.4f}')

    # signature level
    n_up = 150 if tag == 'HCC' else 25
    n_dn = 150 if tag == 'HCC' else 23
    p2g = pd.read_csv(f'{DATA}/{p2g_file}')
    up_l, dn_l = gene_signature(lim, p2g, n_up)
    up_w, dn_w = old_sigs[keymap[tag]]['up'], old_sigs[keymap[tag]]['down']
    for lab, a, b in [('up', set(up_w), set(up_l)), ('down', set(dn_w), set(dn_l))]:
        print(f'  signature {lab}: Welch {len(a)}  limma {len(b)}  overlap {len(a&b)}  '
              f'Jaccard {len(a&b)/len(a|b):.3f}  Welch-only {len(a-b)}  limma-only {len(b-a)}')
    report[tag] = {
        'geo': geo, 'n_common_probes': len(m),
        'deg_welch': len(deg_w), 'deg_limma': len(deg_l), 'deg_overlap': len(inter),
        'deg_jaccard': len(inter) / len(deg_w | deg_l),
        'log2fc_r': float(m.log2FC_w.corr(m.log2FC_l)),
        'max_abs_diff': float((m.log2FC_w - m.log2FC_l).abs().max()),
        # the adjusted-P agreement carries the "identical DE sets" claim just as much as the
        # fold changes do, so it is exported rather than only printed
        'padj_spearman': float(pd.Series(m.padj_w).rank().corr(pd.Series(m.padj_l).rank())),
        'sign_concordance': float(conc),
        'sig_up_overlap': len(set(up_w) & set(up_l)), 'sig_up_welch': len(set(up_w)),
        'sig_up_limma': len(set(up_l)),
        'sig_down_overlap': len(set(dn_w) & set(dn_l)), 'sig_down_welch': len(set(dn_w)),
        'sig_down_limma': len(set(dn_l)),
    }
    # save limma-derived signature for the reversal-score sensitivity test
    json.dump({'up': up_l, 'down': dn_l},
              open(f'{OUT}/A2_limma_signature_{tag}.json', 'w'))

with open(f'{OUT}/A2_comparison.json', 'w') as fh:
    json.dump(report, fh, indent=2)
print(f'\nwrote {OUT}/A2_comparison.json')
