import pandas as pd
import numpy as np
from config import OUT, RAW

def gene_signature(deg_file, p2g_file, n=150):
    deg = pd.read_csv(deg_file)
    deg['probe'] = deg['probe'].astype(str).str.strip('"')
    p2g = pd.read_csv(p2g_file)
    p2g['ID'] = p2g['ID'].astype(str)
    p2g = p2g.dropna(subset=['symbol'])
    p2g = p2g[p2g['symbol'] != '']
    # join probes to symbols
    m = deg.merge(p2g[['ID', 'symbol']], left_on='probe', right_on='ID', how='inner')
    # one row per symbol, keeping the largest |log2FC|
    m['absfc'] = m.log2FC.abs()
    m = m.sort_values('absfc', ascending=False).drop_duplicates('symbol')
    # up / down
    up = m[(m.log2FC > 1) & (m.padj < 0.05)].sort_values('log2FC', ascending=False)
    down = m[(m.log2FC < -1) & (m.padj < 0.05)].sort_values('log2FC', ascending=True)
    up_genes = up['symbol'].tolist()[:n]
    down_genes = down['symbol'].tolist()[:n]
    return up_genes, down_genes

signatures = {}
for name, deg, p2g in [
    ('HCC', f'{RAW}/DEG_HCC.csv', f'{RAW}/GPL3921_probe2gene.csv'),
    ('CAD', f'{RAW}/DEG_CAD.csv', f'{RAW}/GPL20115_probe2gene.csv'),
    ('CVD_HF', f'{RAW}/DEG_CVD.csv', f'{RAW}/GPL11532_probe2gene.csv'),
]:
    up, down = gene_signature(deg, p2g)
    signatures[name] = {'up': up, 'down': down}
    print(f'{name}: {len(up)} up, {len(down)} down')
    print(f'  top5 up: {up[:5]}')
    print(f'  top5 down: {down[:5]}')

import json
with open(f'{RAW}/disease_signatures.json', 'w') as f:
    json.dump(signatures, f)
print(f'\nwrote {RAW}/disease_signatures.json')
