import gzip
import numpy as np
import pandas as pd
from scipy import stats
from config import OUT, RAW

def load_series_matrix(fn):
    """Read a GEO series matrix: (expression DataFrame, sample titles)."""
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
        rows = []
        for line in f:
            parts = line.rstrip('\n').split('\t')
            rows.append(parts)
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(header[0])
    # to numeric
    df = df.apply(pd.to_numeric, errors='coerce')
    return df, titles

def bh_correct(pvals):
    """Benjamini-Hochberg FDR correction."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    adj = np.empty(n)
    adj[order] = sorted_p * n / (np.arange(n) + 1)
    # enforce monotonicity
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return np.clip(adj, 0, 1)

def dea(df, titles, disease_label, control_label):
    """Two-group Welch t-test; returns the DEG table."""
    groups = []
    for t in titles:
        if disease_label(t):
            groups.append('disease')
        elif control_label(t):
            groups.append('control')
        else:
            groups.append('exclude')
    mask = np.array(groups) != 'exclude'
    g = np.array(groups)[mask]
    X = df.iloc[:, mask]
    disease_cols = X.columns[g == 'disease']
    control_cols = X.columns[g == 'control']
    n_d = len(disease_cols)
    n_c = len(control_cols)
    print(f'  disease n={n_d}, control n={n_c}')
    res = []
    for probe in df.index:
        d = X.loc[probe, disease_cols].values.astype(float)
        c = X.loc[probe, control_cols].values.astype(float)
        d = d[~np.isnan(d)]; c = c[~np.isnan(c)]
        if len(d) < 3 or len(c) < 3:
            continue
        mean_d = np.mean(d); mean_c = np.mean(c)
        log2fc = mean_d - mean_c  # already log2, so the difference is the log2FC
        t, p = stats.ttest_ind(d, c, equal_var=False)
        res.append({'probe': probe, 'log2FC': log2fc, 'p': p,
                    'mean_disease': mean_d, 'mean_control': mean_c})
    deg = pd.DataFrame(res)
    deg['padj'] = bh_correct(deg['p'].values)
    return deg

# ---- HCC: GSE14520 ----
print('=== HCC (GSE14520-GPL3921) ===')
df_hcc, titles_hcc = load_series_matrix(f'{RAW}/GSE14520-GPL3921_series_matrix.txt.gz')
deg_hcc = dea(df_hcc, titles_hcc,
              disease_label=lambda t: 'Tumor' in t and 'Non' not in t,
              control_label=lambda t: 'Non-Tumor' in t or 'Non-Tumor' in t)
deg_hcc.to_csv(f'{RAW}/DEG_HCC.csv', index=False)
print(f'HCC DEGs: {len(deg_hcc)}; |log2FC|>1 & padj<0.05: {((deg_hcc.log2FC.abs()>1)&(deg_hcc.padj<0.05)).sum()}')

# ---- CVD: GSE57345 ----
print('=== CVD (GSE57345-GPL11532) ===')
df_cvd, titles_cvd = load_series_matrix(f'{RAW}/GSE57345-GPL11532_series_matrix.txt.gz')
deg_cvd = dea(df_cvd, titles_cvd,
              disease_label=lambda t: ('Dilated' in t) or ('Ischemic' in t),
              control_label=lambda t: ('onfailing' in t) or ('on-failing' in t) or ('ontrol' in t))
deg_cvd.to_csv(f'{RAW}/DEG_CVD.csv', index=False)
print(f'CVD DEGs: {len(deg_cvd)}; |log2FC|>1 & padj<0.05: {((deg_cvd.log2FC.abs()>1)&(deg_cvd.padj<0.05)).sum()}')
