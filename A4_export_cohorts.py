#!/usr/bin/env python3
"""A4 step 1: export independent validation cohorts to CSV for limma.

Cohorts
  GSE76427 (GPL10558, Illumina HumanHT-12 V4)  HCC: 115 tumour vs 52 adjacent
                                               non-tumour liver
  GSE141910 (RNA-seq)                          heart failure; expression has to
                                               come from GEO supplementary files

Group labels are taken from the GEO 'tissue:' characteristic, not from free-text
titles, so the definition is auditable.
"""
import gzip
import os
import numpy as np
import pandas as pd

from config import VAL as DATA, OUT   # validation-cohort inputs


def series_matrix(fn):
    """Return (expression DataFrame probes x samples, sample titles, characteristics)."""
    titles = None
    chars = {}
    with gzip.open(fn, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Sample_title'):
                titles = [t.strip('"') for t in line.split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                vals = [t.strip('"') for t in line.split('\t')[1:]]
                for v in vals:
                    if ':' in v:
                        k = v.split(':', 1)[0].strip()
                        chars.setdefault(k, []).append(v)
            elif line.startswith('!series_matrix_table_begin'):
                break
        header = [h.strip('"') for h in f.readline().rstrip('\n').split('\t')]
        rows = [r for r in (l.rstrip('\n').split('\t') for l in f)
                if r and not r[0].startswith('!')]
    df = pd.DataFrame(rows, columns=header).set_index(header[0])
    df.index = [str(i).strip('"') for i in df.index]
    df = df.apply(pd.to_numeric, errors='coerce').dropna(how='all')
    return df, titles, chars


def export_gse76427():
    df, titles, chars = series_matrix(f'{DATA}/GSE76427_series_matrix.txt.gz')
    tis = [v.replace('tissue:', '').strip() for v in chars.get('tissue', [])]
    assert len(tis) == df.shape[1], (len(tis), df.shape)
    grp = np.array(['tumour' if 'tumor' in t and 'non' not in t
                    else ('normal' if 'non-tumor' in t else 'exclude') for t in tis])
    keep = grp != 'exclude'
    expr = df.iloc[:, keep]
    expr.index.name = 'probe'
    # GSE76427 ships linear (background-subtracted) intensities, median ~120 with
    # negative values, not log2 like GSE14520. Convert with the standard
    # log2(pmax(x,0)+1); between-array quantile normalisation is applied in limma.
    raw_med = float(np.nanmedian(expr.values))
    expr = np.log2(expr.clip(lower=0) + 1)
    expr.to_csv(f'{OUT}/A4_GSE76427_expr.csv')
    pd.DataFrame({'sample': expr.columns, 'group': grp[keep]}).to_csv(
        f'{OUT}/A4_GSE76427_groups.csv', index=False)
    print(f'GSE76427: probes {expr.shape[0]}  tumour {int((grp=="tumour").sum())}  '
          f'normal {int((grp=="normal").sum())}')
    print(f'  raw median {raw_med:.1f} (linear) -> log2(pmax(x,0)+1): '
          f'{np.nanmin(expr.values):.2f}-{np.nanmax(expr.values):.2f}, '
          f'median {np.nanmedian(expr.values):.2f}')
    # probe -> gene symbol from the GEO annotation file
    with gzip.open(f'{DATA}/GPL10558.annot.gz', 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):
                break
        hdr = f.readline().rstrip('\n').split('\t')
        rows = [l.rstrip('\n').split('\t') for l in f
                if not l.startswith('!')]
    ann = pd.DataFrame(rows, columns=hdr)
    ann = ann[['ID', 'Gene symbol']].rename(columns={'ID': 'probe', 'Gene symbol': 'symbol'})
    ann = ann[ann.symbol.astype(str).str.strip() != '']
    ann.to_csv(f'{OUT}/A4_GPL10558_probe2gene.csv', index=False)
    print(f'GPL10558 probe->gene: {len(ann)} pairs, covering '
          f'{ann.probe.isin(expr.index.astype(str)).sum()}')


def export_gse141910():
    """Third cohort: GSE141910, per-sample RNA-seq matrices in a RAW tar.

    Disease label follows the discovery cohort definition (DCM + ICM vs
    non-failing donors); hypertrophic and peripartum aetiologies are excluded so
    the contrast stays comparable to GSE57345.
    """
    import glob
    files = sorted(glob.glob(f'{DATA}/GSE141910_raw/*.csv.gz'))
    frames, order = {}, []
    for fn in files:
        sid = os.path.basename(fn).split('_', 1)[1].replace('.csv.gz', '')
        s = pd.read_csv(fn, index_col=0, header=0)
        frames[sid] = s.iloc[:, 0].astype(float)
        order.append(sid)
    expr = pd.DataFrame(frames)
    print(f'GSE141910: merged {len(order)} per-sample files, {expr.shape[0]} ENSG ids')

    _, titles, chars = series_matrix(f'{DATA}/GSE141910_series_matrix.txt.gz')
    et = [v.replace('etiology:', '').strip() for v in chars.get('etiology', [])]
    lab = pd.DataFrame({'sample': titles, 'etiology': et})
    dcm = ['Dilated cardiomyopathy (DCM)', 'Ischemic cardiomyopathy (ICM)']
    grp = np.where(lab.etiology.isin(dcm), 'disease',
                   np.where(lab.etiology == 'Non-Failing Donor', 'normal', 'exclude'))
    lab['group'] = grp
    print('  etiology x group:', lab.groupby(['etiology', 'group']).size().to_dict())

    cols = [c for c in expr.columns if c in set(lab['sample'])]
    expr = expr[cols]
    lab = lab.set_index('sample').loc[cols].reset_index()
    keep = lab['group'] != 'exclude'
    expr = expr.loc[:, keep.values]

    # ENSG -> symbol via the GSE116250 RPKM table (same GENCODE-style ids)
    ann = pd.read_csv(f'{DATA}/GSE116250_rpkm.txt.gz', sep='\t',
                      usecols=['Gene', 'Common_name'])
    m = dict(zip(ann.Gene.astype(str), ann.Common_name.astype(str)))
    expr.index = [m.get(i, i) for i in expr.index]
    sym = pd.Series(expr.index)
    expr = expr[(sym.str.strip().ne('') & sym.ne('nan')).values]
    expr = expr.assign(_mean=expr.mean(axis=1)).sort_values('_mean', ascending=False)
    expr = expr[~expr.index.duplicated(keep='first')].drop(columns='_mean')
    expr.index.name = 'symbol'
    expr.to_csv(f'{OUT}/A4_GSE141910_expr.csv')
    lab[keep].to_csv(f'{OUT}/A4_GSE141910_groups.csv', index=False)
    print(f'GSE141910: {expr.shape[0]} genes x {expr.shape[1]} samples  '
          f'disease {int((lab["group"]=="disease").sum())}  '
          f'normal {int((lab["group"]=="normal").sum())}')


def probe_rnaseq_supplementary(acc):
    """List supplementary files for an RNA-seq GEO accession (expression lives there)."""
    import urllib.request
    import re
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief'
    with urllib.request.urlopen(url, timeout=60) as r:
        txt = r.read().decode('utf-8', 'replace')
    files = re.findall(r'!Series_supplementary_file = (.+)', txt)
    print(f'{acc} supplementary files:')
    for x in files:
        print('   ', x.strip())
    return [x.strip() for x in files]


if __name__ == '__main__':
    export_gse76427()
    try:
        export_gse141910()
    except Exception as e:                           # noqa: BLE001
        print(f'GSE141910: failed: {e}')
