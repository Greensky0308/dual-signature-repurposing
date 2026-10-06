import gzip
import numpy as np
import pandas as pd
from scipy import stats
from config import OUT, RAW

def load_series_matrix(fn):
    titles = None
    with gzip.open(fn, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                break
        header = f.readline().rstrip('\n').split('\t')
        header = [h.strip('"') for h in header]
        rows = [line.rstrip('\n').split('\t') for line in f]
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    return df.apply(pd.to_numeric, errors='coerce'), titles

def bh_correct(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(n)
    adj[order] = pvals[order] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return np.clip(adj, 0, 1)

def dea(df, labels):
    labels = np.array(labels)
    mask = labels != 'exclude'
    g = labels[mask]
    X = df.iloc[:, mask]
    d_cols = X.columns[g == 'disease']
    c_cols = X.columns[g == 'control']
    print(f'  disease n={len(d_cols)}, control n={len(c_cols)}')
    res = []
    for probe in df.index:
        d = X.loc[probe, d_cols].values.astype(float)
        c = X.loc[probe, c_cols].values.astype(float)
        d = d[~np.isnan(d)]; c = c[~np.isnan(c)]
        if len(d) < 3 or len(c) < 3:
            continue
        log2fc = np.mean(d) - np.mean(c)
        t, p = stats.ttest_ind(d, c, equal_var=False)
        res.append({'probe': probe, 'log2FC': log2fc, 'p': p})
    deg = pd.DataFrame(res)
    deg['padj'] = bh_correct(deg['p'].values)
    return deg

# CAD
print('=== CAD (GSE113079) ===')
df_cad, titles_cad = load_series_matrix(f'{RAW}/GSE113079_series_matrix.txt.gz')
labels_cad = []
for t in titles_cad:
    if 'CAD' in t:
        labels_cad.append('disease')
    elif 'control' in t.lower() or 'healthy' in t.lower() or 'normal' in t.lower():
        labels_cad.append('control')
    else:
        labels_cad.append('exclude')
deg_cad = dea(df_cad, labels_cad)
deg_cad.to_csv(f'{RAW}/DEG_CAD.csv', index=False)
n_cad = ((deg_cad.log2FC.abs()>1)&(deg_cad.padj<0.05)).sum()
print(f'CAD significant DEGs (|log2FC|>1 & padj<0.05): {n_cad}')
