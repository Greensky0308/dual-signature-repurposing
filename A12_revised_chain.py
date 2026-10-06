#!/usr/bin/env python3
"""A12: single source of truth for both candidate chains.

PRIMARY chain (reported in the manuscript, Figure 1)
    named compounds -> dual-reversal -> touchstone -> documented indication.
    Two of the three nodes keep their published values (867, 384); only the
    annotation node is corrected (40 -> 118), because the published 40 came from
    an annotation table that is not traceable to a public database.

SENSITIVITY chain (reported in the response letter + supplementary)
    The published screen additionally required a PubChem CID in the LINCS
    metadata. That requirement has no scientific motivation and removes, for
    example, digitoxin and ouabain, both launched cardiac glycosides. Dropping it
    and using the Hub's curated target field in place of the LINCS touchstone
    flag gives 1,139 -> 483 -> 224 -- the analysis in which the launched cardiac
    glycoside digitoxin appears.

Cardiac glycosides are therefore reported as a mechanism class evidenced by
(a) their presence among the dual-reversal touchstone compounds and (b) the
sensitivity chain, not by an indication annotation.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

from config import RAW as DATA, OUT, HUB

hub = pd.read_csv(HUB, sep='\t', comment='!', low_memory=False).drop_duplicates('pert_iname')
hub['has_target'] = hub.target.notna()
hub['has_ind'] = hub.indication.notna()
hub['launched'] = hub.clinical_phase.eq('Launched')

scores = pd.read_csv(f'{DATA}/drug_scores_all3.csv')[
    ['true_name', 'hcc_score', 'cvd_score', 'has_name', 'pubchem']]
scores = scores[scores.true_name.astype(str) != 'nan']
scores['dual'] = (scores.hcc_score > 0) & (scores.cvd_score > 0)
scores['named'] = scores.has_name.astype(bool)
scores['has_cid'] = scores.pubchem.astype(str) != '-666'

df = scores.merge(hub[['pert_iname', 'clinical_phase', 'moa', 'target', 'disease_area',
                       'indication', 'has_target', 'has_ind', 'launched']],
                  left_on='true_name', right_on='pert_iname', how='left')
for c in ('has_target', 'has_ind', 'launched'):
    df[c] = df[c].fillna(False).astype(bool)

CG = ['digitoxin', 'digoxin', 'ouabain', 'peruvoside', 'sarmentogenin', 'strophanthidin',
      'cinobufagin', 'bufalin', 'gitoxigenin', 'digitoxigenin', 'digoxigenin',
      'proscillaridin', 'proscillaridin-a']


def class_masks(frame):
    m = frame.moa.astype(str)
    return {
        'HDAC inhibitor': m.str.contains('HDAC|deacetylase', case=False, na=False),
        'Topoisomerase inhibitor': m.str.contains('topoisomerase', case=False, na=False),
        'Tubulin inhibitor': m.str.contains('tubulin', case=False, na=False),
        'CDK inhibitor': m.str.contains('CDK|cyclin', case=False, na=False),
        'Proteasome inhibitor': m.str.contains('proteasome', case=False, na=False),
        'ATP1A1 target': frame.target.astype(str).str.contains('ATP1A1', na=False),
    }


def enrich(frame, label):
    """Fisher exact test per mechanism class: candidates vs the rest of the library."""
    rows = []
    for name, member in class_masks(frame).items():
        a = int((member & frame.in_candidates).sum())
        b = int((member & ~frame.in_candidates).sum())
        c = int((~member & frame.in_candidates).sum())
        d = int((~member & ~frame.in_candidates).sum())
        if a + b == 0:
            continue
        orr, p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append(dict(chain=label, mechanism_class=name, n_library=a + b,
                         n_candidates=a, odds_ratio=orr, p=p,
                         members='; '.join(sorted(
                             frame.loc[member & frame.in_candidates, 'true_name']))[:300]))
    r = pd.DataFrame(rows).sort_values('p').reset_index(drop=True)
    if len(r):
        m = len(r)
        r['padj'] = np.minimum.accumulate((r.p.values * m / (np.arange(m) + 1))[::-1])[::-1]
        r['padj'] = r.padj.clip(upper=1.0)
    return r


# ---- PRIMARY chain -------------------------------------------------------
tch = set(pd.read_csv(f'{DATA}/candidate_drugs_touchstone.csv').true_name.astype(str))
df['is_touch'] = df.true_name.isin(tch)
df['p_final'] = df.named & df.has_cid & df.dual & df.is_touch & df.has_ind
p = dict(n_groups=int(len(df)),
         named=int(df.named.sum()),
         dual_named=int((df.named & df.dual).sum()),
         dual_cid=int((df.named & df.has_cid & df.dual).sum()),
         touchstone=int((df.named & df.has_cid & df.dual & df.is_touch).sum()),
         final=int(df.p_final.sum()))

print('=== PRIMARY chain (reported in the manuscript) ===')
print(f'  compound groups                    {p["n_groups"]:,}')
print(f'  with a chemical name               {p["named"]:,}')
print(f'  dual-reversal (named)              {p["dual_named"]:,}')
print(f'  dual-reversal with PubChem CID     {p["dual_cid"]:,}   (published 867)')
print(f'  touchstone                         {p["touchstone"]:,}   (published 384)')
print(f'  documented indication (Hub)        {p["final"]:,}   (published 40, hand-assembled)')

df['in_candidates'] = df.p_final
primary_enr = enrich(df, 'primary')

pr = df[df.p_final].copy()
primary_rows = []
for name, member in class_masks(pr).items():
    n = int(member.sum())
    if n:
        primary_rows.append(dict(mechanism_class=name, n=n,
                                 members='; '.join(sorted(pr.loc[member, 'true_name']))[:300]))
primary_cls = pd.DataFrame(primary_rows).sort_values('n', ascending=False)
print('\n  enrichment in the primary chain:')
print(primary_enr[['mechanism_class', 'n_library', 'n_candidates', 'odds_ratio', 'padj']]
      .round(4).to_string(index=False))
print('\n  classes among the final node:')
for _, r in primary_cls.iterrows():
    print(f'    {r.mechanism_class:26s} {r.n:3d}   {r.members[:110]}')
cg_384 = sorted(set(df.loc[df.named & df.has_cid & df.dual & df.is_touch, 'true_name'])
                & set(CG))
print(f'\n  cardiac glycosides among the 384: {cg_384}')

# ---- SENSITIVITY chain ---------------------------------------------------
df['s_final'] = df.named & df.dual & df.has_target & df.has_ind
s = dict(n_groups=int(len(df)), named=int(df.named.sum()),
         dual=int((df.named & df.dual).sum()),
         target=int((df.named & df.dual & df.has_target).sum()),
         final=int(df.s_final.sum()))

print('\n=== SENSITIVITY chain (PubChem-CID filter removed; Hub target for touchstone) ===')
print(f'  compound groups                    {s["n_groups"]:,}')
print(f'  with a chemical name               {s["named"]:,}')
print(f'  dual-reversal (no CID filter)      {s["dual"]:,}')
print(f'  + curated Hub target               {s["target"]:,}')
print(f'  + documented indication            {s["final"]:,}')

df['in_candidates'] = df.s_final
sens_enr = enrich(df, 'sensitivity')
print('\n  enrichment in the sensitivity chain:')
print(sens_enr[['mechanism_class', 'n_library', 'n_candidates', 'odds_ratio', 'padj']]
      .round(4).to_string(index=False))
s_cg = sorted(set(df.loc[df.s_final, 'true_name']) & set(CG))
print(f'\n  cardiac glycosides in the final node: {s_cg}')
for d in ['digitoxin', 'digoxin', 'ouabain']:
    r = df[df.true_name == d]
    if len(r):
        r = r.iloc[0]
        print(f'    {d:12s} dual={"yes" if r.dual else "no"}  named={bool(r.named)}  '
              f'cid={bool(r.has_cid)}  phase={str(r.clinical_phase):10s} '
              f'HCC {r.hcc_score:+.3f}  HF {r.cvd_score:+.3f}')

# ---- write --------------------------------------------------------------
summary = {
    'primary': dict(p, classes=primary_cls.to_dict('records'),
                    cardiac_glycosides_in_touchstone=cg_384,
                    enrichment=primary_enr.to_dict('records')),
    'sensitivity': dict(s, cardiac_glycosides=s_cg,
                        enrichment=sens_enr.to_dict('records')),
    'n_profiles': 473647,
    'annotation_source': 'Broad Drug Repurposing Hub, repurposing_drugs_20200324.txt',
}
with open(f'{OUT}/A12_revised_chain.json', 'w') as fh:
    json.dump(summary, fh, indent=2, default=float)
pr[['true_name', 'hcc_score', 'cvd_score', 'clinical_phase', 'moa', 'target',
    'disease_area', 'indication']].sort_values('cvd_score', ascending=False).to_csv(
    f'{OUT}/A12_primary_final_118.csv', index=False)
df[df.s_final][['true_name', 'hcc_score', 'cvd_score', 'clinical_phase', 'moa', 'target',
                'disease_area', 'indication']].sort_values('cvd_score', ascending=False).to_csv(
    f'{OUT}/A12_sensitivity_224.csv', index=False)
print(f'\nwrote {OUT}/A12_revised_chain.json, A12_primary_final_118.csv, A12_sensitivity_224.csv')
