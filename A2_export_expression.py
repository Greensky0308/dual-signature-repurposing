#!/usr/bin/env python3
"""A2 step 1: export the two series matrices and group labels for limma.

Writes probe x sample expression CSVs plus a two-column group file, using the
exact same sample-label rules as the original Welch's t-test pipeline
(raw/data/DEG_HCC.csv, raw/data/DEG_CVD.csv).
"""
import gzip
import numpy as np
import pandas as pd

from config import RAW as DATA, OUT


def load_series_matrix(fn):
    titles = None
    with gzip.open(fn, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.split('\t')[1:]]
            elif line.startswith('!series_matrix_table_begin'):
                break
        header = [h.strip('"') for h in f.readline().rstrip('\n').split('\t')]
        rows = [line.rstrip('\n').split('\t') for line in f]
    rows = [r for r in rows if r and not r[0].startswith('!')]
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna(how='all')
    return df, titles


JOBS = [
    ('GSE14520-GPL3921_series_matrix.txt.gz', 'GSE14520',
     lambda t: 'Tumor' in t and 'Non' not in t,
     lambda t: 'Non-Tumor' in t),
    ('GSE57345-GPL11532_series_matrix.txt.gz', 'GSE57345',
     lambda t: ('Dilated' in t) or ('Ischemic' in t),
     lambda t: ('onfailing' in t) or ('on-failing' in t) or ('ontrol' in t)),
]

for fn, tag, dis_lab, ctl_lab in JOBS:
    df, titles = load_series_matrix(f'{DATA}/{fn}')
    grp = []
    for t in titles:
        if dis_lab(t):
            grp.append('disease')
        elif ctl_lab(t):
            grp.append('control')
        else:
            grp.append('exclude')
    grp = np.array(grp)
    mask = grp != 'exclude'
    sub = df.iloc[:, mask]
    expr = sub.copy()
    expr.columns = [c for c in sub.columns]
    expr.index.name = 'probe'
    expr.to_csv(f'{OUT}/A2_{tag}_expr.csv')
    pd.DataFrame({'sample': sub.columns, 'group': grp[mask]}).to_csv(
        f'{OUT}/A2_{tag}_groups.csv', index=False)
    n_d = int((grp == 'disease').sum())
    n_c = int((grp == 'control').sum())
    print(f'{tag}: probes {expr.shape[0]}  disease {n_d}  control {n_c}  '
          f'value range {np.nanmin(expr.values):.2f}-{np.nanmax(expr.values):.2f}  '
          f'median {np.nanmedian(expr.values):.2f}')
