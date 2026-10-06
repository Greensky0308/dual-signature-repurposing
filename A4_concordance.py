#!/usr/bin/env python3
"""A4 step 3: are the discovery signatures reproduced in independent cohorts?

For each (signature, cohort) pair this reports four things:
  1. how many signature genes are measurable in the cohort;
  2. the directional concordance rate with a binomial test against 50%;
  3. a competitive gene-set test (one-sided Mann-Whitney on limma t-statistics,
     signature genes vs all other genes) -- the analogue of limma::geneSetTest;
  4. how well the signature is recovered de novo (Jaccard against the top gene
     list rebuilt from the validation cohort itself).
"""
import json
import os
import numpy as np
import pandas as pd
from scipy import stats

from config import RAW as DATA, OUT

PAIRS = [
    dict(name='HCC signature -> GSE76427 (tumour vs adjacent)',
         sig_key='HCC', deg='A4_limma_GSE76427.csv', p2g='A4_GPL10558_probe2gene.csv',
         n_up=150, n_dn=150),
    dict(name='HF signature -> GSE141910 (DCM vs non-failing donor)',
         sig_key='CVD_HF', deg='A4_limma_GSE141910.csv', p2g='A4_GSE141910_probe2gene.csv',
         n_up=25, n_dn=23),
]


def load_gene_level(deg_file, p2g_file):
    deg = pd.read_csv(f'{OUT}/{deg_file}')
    deg['probe'] = deg['probe'].astype(str)
    p2g = pd.read_csv(f'{OUT}/{p2g_file}', dtype=str).dropna()
    m = deg.merge(p2g, on='probe', how='inner')
    m = m[m.symbol.astype(str).str.strip() != '']
    m['absfc'] = m.log2FC.abs()
    m = m.sort_values('absfc', ascending=False).drop_duplicates('symbol')
    return m


def gene_signature(m, n_up, n_dn):
    up = m[(m.log2FC > 1) & (m.padj < 0.05)].sort_values('log2FC', ascending=False)
    dn = m[(m.log2FC < -1) & (m.padj < 0.05)].sort_values('log2FC', ascending=True)
    return up.symbol.tolist()[:n_up], dn.symbol.tolist()[:n_dn]


def main():
    # the RNA-seq cohort's rows are already gene symbols -> an identity mapping file
    for tag in ('GSE141910',):
        ident = f'{OUT}/A4_{tag}_probe2gene.csv'
        if not os.path.exists(ident):
            d = pd.read_csv(f'{OUT}/A4_limma_{tag}.csv')
            pd.DataFrame({'probe': d.probe.astype(str),
                          'symbol': d.probe.astype(str)}).to_csv(ident, index=False)
            print(f'wrote an identity probe map: {ident}')

    sigs = json.load(open(f'{DATA}/disease_signatures.json'))
    report = {}
    for p in PAIRS:
        m = load_gene_level(p['deg'], p['p2g'])
        tmap = dict(zip(m.symbol, m.t_limma))
        fc = dict(zip(m.symbol, m.log2FC))
        allt = m.t_limma.values
        sig = sigs[p['sig_key']]
        up, dn = sig['up'], sig['down']
        print(f"\n===== {p['name']} =====")
        print(f'probes mapping to a gene: {len(m)} ({m.symbol.nunique()} genes)')

        res = {}
        for lab, genes, direction in [('up', up, 'greater'), ('down', dn, 'less')]:
            present = [g for g in genes if g in tmap]
            exp_sign = 1 if lab == 'up' else -1
            conc = [g for g in present if np.sign(fc[g]) == exp_sign]
            k, n = len(conc), len(present)
            pbin = stats.binomtest(k, n, 0.5, alternative='greater').pvalue
            tvals = np.array([tmap[g] for g in present])
            mw = stats.mannwhitneyu(tvals, allt, alternative=direction)
            print(f'\n{lab}: signature of {len(genes)} genes, {n} measurable here')
            print(f'  concordant direction: {k}/{n} = {100*k/max(n,1):.1f}%   '
                  f'binomial p = {pbin:.3e}')
            print(f'  median limma t {np.median(tvals):+.3f} (all genes '
                  f'{np.median(allt):+.3f})')
            print(f'  competitive gene-set test (Mann-Whitney, {direction}): p = {mw.pvalue:.3e}')
            res[lab] = dict(n_present=n, n_concordant=k, concordance=k/max(n, 1),
                            binom_p=float(pbin), median_t=float(np.median(tvals)),
                            geneset_p=float(mw.pvalue))
        # de novo recovery
        up2, dn2 = gene_signature(m, p['n_up'], p['n_dn'])
        for lab, a, b in [('up', set(up), set(up2)), ('down', set(dn), set(dn2))]:
            jac = len(a & b) / max(len(a | b), 1)
            print(f'\nrebuilt {lab} signature de novo: {p["n_up"] if lab=="up" else p["n_dn"]} genes, '
                  f'{len(a&b)} shared with the discovery signature, Jaccard {jac:.3f}')
            res[f'denovo_{lab}'] = dict(overlap=len(a & b), jaccard=jac)
        report[p['name']] = res

    with open(f'{OUT}/A4_concordance.json', 'w') as fh:
        json.dump(report, fh, indent=2, default=float)
    print(f'\nwrote {OUT}/A4_concordance.json')


if __name__ == '__main__':
    main()
