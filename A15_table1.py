#!/usr/bin/env python3
"""A15: rebuild Table 1 from the Hub annotation.

Part A  representative annotated candidates (the final node of the primary
        chain), ordered by combined reversal score, with target, indication and
        mechanism class taken from the Broad Drug Repurposing Hub.
Part B  the cardiac-glycoside / ATP1A1 entries, listed separately as in the
        published table, with their evidence level stated: which layer each one
        reaches, and whether its target or indication is documented.

Outputs: tables/Table1_revised.csv and Table1_revised.tsv
"""
import os
import pandas as pd

from config import RAW as DATA, OUT, HUB, TBL

os.makedirs(TBL, exist_ok=True)

hub = pd.read_csv(HUB, sep='\t', comment='!', low_memory=False).drop_duplicates('pert_iname')
hub['has_target'] = hub.target.notna()
hub['has_ind'] = hub.indication.notna()
scores = pd.read_csv(f'{DATA}/drug_scores_all3.csv')
scores = scores[scores.true_name.astype(str) != 'nan']
scores['dual'] = (scores.hcc_score > 0) & (scores.cvd_score > 0)
scores['named'] = scores.has_name.astype(bool)
scores['has_cid'] = scores.pubchem.astype(str) != '-666'
scores['combined'] = scores.hcc_score + scores.cvd_score

tch = set(pd.read_csv(f'{DATA}/candidate_drugs_touchstone.csv').true_name.astype(str))
# The Hub spells some compounds differently from LINCS, and an exact-name join
# drops their whole annotation without a word (proscillaridin is recorded as
# proscillaridin-A: Launched, with a congestive-heart-failure indication).
HUB_ALIAS = {'proscillaridin': 'proscillaridin-A'}
scores['hub_name'] = [HUB_ALIAS.get(str(n).lower(), n) for n in scores.true_name]
df = scores.merge(hub, left_on='hub_name', right_on='pert_iname', how='left')
df['is_touch'] = df.true_name.isin(tch)
final = df[df.named & df.has_cid & df.dual & df.is_touch &
           df.has_ind.fillna(False).astype(bool)].copy()


def mech(row):
    m = str(row.get('moa', '')).lower()
    t = str(row.get('target', ''))
    for key, lab in [('hdac', 'HDAC inhibitor'), ('deacetylase', 'HDAC inhibitor'),
                     ('topoisomerase', 'Topoisomerase inhibitor'),
                     ('tubulin', 'Tubulin inhibitor'),
                     ('cdk', 'CDK inhibitor'), ('cyclin', 'CDK inhibitor'),
                     ('proteasome', 'Proteasome inhibitor')]:
        if key in m:
            return lab
    if 'ATP1A1' in t:
        return 'Cardiac glycoside / ATP1A1'
    return 'Other'


final['mechanism_class'] = final.apply(mech, axis=1)
partA = final.sort_values('combined', ascending=False)[
    ['true_name', 'hcc_score', 'cvd_score', 'combined', 'target', 'indication',
     'mechanism_class', 'clinical_phase']].head(12)
partA.columns = ['Drug', 'HCC score', 'HF score', 'Combined score', 'Target (Hub)',
                 'Indication (Hub)', 'Mechanism class', 'Clinical phase']
partA = partA.round(3)

# ---- Part B: cardiac glycosides / ATP1A1 --------------------------------
CG = ['digitoxin', 'digoxin', 'ouabain', 'peruvoside', 'sarmentogenin',
      'strophanthidin', 'cinobufagin', 'bufalin', 'proscillaridin',
      'digitoxigenin', 'gitoxigenin', 'digoxigenin']
rows = []
for d in CG:
    r = df[df.true_name == d]
    if not len(r):
        continue
    r = r.iloc[0]
    flag = lambda ok: 'yes' if ok else 'no'          # noqa: E731
    dual = bool(r.dual)
    touch = bool(r.is_touch)
    cid = bool(r.has_cid)
    ind = bool(r.has_ind) if isinstance(r.has_ind, bool) else False
    rows.append({
        'Drug': d,
        'HCC score': round(float(r.hcc_score), 3),
        'HF score': round(float(r.cvd_score), 3),
        'Combined score': round(float(r.hcc_score + r.cvd_score), 3),
        'Target (Hub)': r.target if isinstance(r.target, str) else 'not annotated',
        'Indication (Hub)': r.indication if isinstance(r.indication, str) else 'not annotated',
        'Clinical phase': r.clinical_phase if isinstance(r.clinical_phase, str) else 'not annotated',
        'Dual reversal': flag(dual),
        'Touchstone': flag(touch),
        'PubChem CID': flag(cid),
        'Documented indication': flag(ind),
        'Passes primary chain': flag(dual and touch and cid and ind),
    })
partB = pd.DataFrame(rows)

print('=== Part A: representative annotated candidates (final node of the primary chain) ===')
print(partA.to_string(index=False))
print('\n=== Part B: cardiac glycosides / ATP1A1 ===')
print(partB.to_string(index=False))

partA.to_csv(f'{TBL}/Table1_revised_partA.csv', index=False)
partB.to_csv(f'{TBL}/Table1_revised_partB.csv', index=False)
with open(f'{TBL}/Table1_revised.tsv', 'w') as fh:
    fh.write('A. Representative annotated candidates\n')
    partA.to_csv(fh, sep='\t', index=False)
    fh.write('\nB. Cardiac glycosides and ATP1A1 interactors\n')
    partB.to_csv(fh, sep='\t', index=False)
print(f'\nwrote {TBL}/Table1_revised.csv and .tsv (partA + partB)')
