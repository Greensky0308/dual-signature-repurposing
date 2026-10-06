import json
import pandas as pd
import os

import gseapy as gp
from config import OUT, RAW, HUB

GENESETS = os.path.join(os.path.dirname(HUB), 'genesets', 'GO_Biological_Process_2023.gmt')

signatures = json.load(open(f'{RAW}/disease_signatures.json'))

results = {}
for disease, genes in signatures.items():
    for direction, gene_list in [('up', genes['up']), ('down', genes['down'])]:
        if not gene_list:
            continue
        try:
            # the gene sets are frozen in data/genesets so the panel does not depend on the
            # Enrichr library as served at run time
            enr = gp.enrich(gene_list=gene_list, gene_sets=GENESETS, outdir=None)
            res = enr.results.head(10)
            results[f'{disease}_{direction}'] = res
            print(f'=== {disease} {direction} (n={len(gene_list)}) ===')
            print(res[['Term', 'Adjusted P-value', 'Overlap']].head(5).to_string(index=False))
        except Exception as e:
            print(f'{disease} {direction} enrichment failed: {repr(e)[:100]}')

# save
with open(f'{RAW}/enrichment_results.json', 'w') as f:
    json.dump({k: v[['Term', 'Adjusted P-value', 'Overlap']].to_dict('records') for k, v in results.items()}, f)
print(f'\nwrote {RAW}/enrichment_results.json')
