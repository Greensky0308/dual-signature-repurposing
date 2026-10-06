#!/usr/bin/env python3
"""A5: quantitative enrichment analysis of mechanism classes.

R1 asked for "additional independent resources or quantitative enrichment
analyses" to reduce false-positive concerns, and for enrichment analysis of the
scoring framework (R1-3). This script tests, with Fisher's exact test, whether
each mechanism class is over-represented among the dual-reversal candidates
relative to the annotated L1000 library, and reports Benjamini-Hochberg adjusted
p-values.

Inputs are the manuscript-pipeline compound scores (output/A1_observed_drug_scores.csv,
2,947 annotated compounds) and the Broad Drug Repurposing Hub annotations.
"""
import numpy as np
import pandas as pd
from scipy import stats

from config import OUT, HUB

scores = pd.read_csv(f'{OUT}/A1_observed_drug_scores.csv')      # compound, hcc_score, hf_score
hub = pd.read_csv(HUB, sep='\t', comment='!', low_memory=False)
hub['has_ind'] = hub.indication.notna()
hub = hub.drop_duplicates('pert_iname')

df = scores.merge(hub[['pert_iname', 'clinical_phase', 'moa', 'target',
                       'disease_area', 'indication', 'has_ind']],
                  left_on='compound', right_on='pert_iname', how='left')
df['dual'] = (df.hcc_score > 0) & (df.hf_score > 0)
df['annotated'] = df.has_ind.fillna(False).astype(bool)
df['launched'] = df.clinical_phase.eq('Launched')
print(f'compounds {len(df)}; dual-positive {int(df.dual.sum())}; '
      f'Hub-annotated indication {int(df.annotated.sum())}; dual-positive with annotation {int((df.dual & df.annotated).sum())}')

# mechanism classes from the Hub: moa keywords and target gene symbols
MOA_CLASSES = {
    'HDAC inhibitor': r'HDAC|deacetylase',
    'CDK inhibitor': r'CDK|cyclin-dependent',
    'Topoisomerase inhibitor': r'topoisomerase',
    'Proteasome inhibitor': r'proteasome',
    'ATPase inhibitor': r'ATPase',
    'Kinase inhibitor (any)': r'kinase',
    'Tubulin inhibitor': r'tubulin',
    'HSP inhibitor': r'heat shock|HSP',
    'DNA methyltransferase inhibitor': r'DNA methyltransferase|DNMT',
    'Aromatase inhibitor': r'aromatase',
}
classes = {k: df.moa.astype(str).str.contains(v, case=False, na=False) for k, v in MOA_CLASSES.items()}
classes['Cardiac glycoside (Hub moa)'] = df.moa.astype(str).str.contains(
    'cardiac glycoside', case=False, na=False)
classes['ATP1A1 target'] = df.target.astype(str).str.contains('ATP1A1', na=False)

rows = []
for name, member in classes.items():
    a = int((member & df.dual).sum())
    b = int((member & ~df.dual).sum())
    c = int((~member & df.dual).sum())
    d = int((~member & ~df.dual).sum())
    if a + b == 0:
        continue
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    rate_in = a / (a + b)
    rate_out = c / (c + d) if (c + d) else float('nan')
    rows.append(dict(mechanism_class=name, n_in_library=a + b, n_dual=a,
                     rate_in=rate_in, rate_out=rate_out, odds_ratio=orr, p=p))

res = pd.DataFrame(rows).sort_values('p')
# Benjamini-Hochberg
m = len(res)
res['padj'] = np.minimum.accumulate((res.p.values * m / (np.arange(m) + 1))[::-1])[::-1]
res['padj'] = res.padj.clip(upper=1.0)
print('\n=== mechanism-class enrichment among dual-reversal candidates (one-sided Fisher) ===')
print(res.round(4).to_string(index=False))

# same test against the annotation-filtered candidate set (annotation-based
# prioritization, i.e. dual-positive AND carrying a Hub indication)
df['principled'] = df.dual & df.annotated
rows2 = []
for name, member in classes.items():
    a = int((member & df.principled).sum())
    b = int((member & ~df.principled).sum())
    c = int((~member & df.principled).sum())
    d = int((~member & ~df.principled).sum())
    if a + b == 0:
        continue
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    rows2.append(dict(mechanism_class=name, n_in_library=a + b, n_principled=a,
                      rate_in=a / (a + b),
                      rate_out=c / (c + d) if (c + d) else float('nan'),
                      odds_ratio=orr, p=p))
res2 = pd.DataFrame(rows2).sort_values('p')
m2 = len(res2)
res2['padj'] = np.minimum.accumulate((res2.p.values * m2 / (np.arange(m2) + 1))[::-1])[::-1]
res2['padj'] = res2.padj.clip(upper=1.0)
print('\n=== mechanism-class enrichment among prioritised candidates (dual-positive with a Hub indication) ===')
print(res2.round(4).to_string(index=False))

res.to_csv(f'{OUT}/A5_enrichment_dual.csv', index=False)
res2.to_csv(f'{OUT}/A5_enrichment_principled.csv', index=False)

# launched-drug enrichment: are approved drugs over-represented?
for tag, mask in [('dual-reversal candidates', df.dual), ('prioritised candidates', df.principled)]:
    a = int((df.launched & mask).sum()); b = int((df.launched & ~mask).sum())
    c = int((~df.launched & mask).sum()); d = int((~df.launched & ~mask).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    print(f'\n[launched] {tag}: {a}/{a+b} in vs {c}/{c+d} out  OR={orr:.2f}  p={p:.4g}')

print(f'\nwrote {OUT}/A5_enrichment_dual.csv and A5_enrichment_principled.csv')
