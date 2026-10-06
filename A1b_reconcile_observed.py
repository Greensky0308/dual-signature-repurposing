#!/usr/bin/env python3
"""A1 step 0: reconcile the recomputed observed scores against the manuscript's
published intermediate files (raw/data/drug_scores_all3.csv, candidate_drugs.csv).

Checks (a) which signature genes fail to map to the L1000 gene axis,
(b) agreement of per-compound scores, (c) how the dual-positive count depends on
the aggregation rule (mean vs lower median) and on the compound universe.
"""
import json
import numpy as np
import pandas as pd
import h5py

from config import RAW as DATA, OUT

gene_info = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_gene_info.txt.gz', sep='\t', compression='gzip')
gene_info['pr_gene_id'] = gene_info['pr_gene_id'].astype(str)
sig_info = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_sig_info.txt.gz', sep='\t',
                       compression='gzip', low_memory=False)
sig_info['sig_id'] = sig_info['sig_id'].astype(str)
pert_info = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_pert_info.txt.gz', sep='\t',
                        compression='gzip', low_memory=False)

f = h5py.File(f'{DATA}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx', 'r')
gene_ids = [x.decode() for x in f['/0/META/ROW/id'][:]]
prof_ids = [x.decode() for x in f['/0/META/COL/id'][:]]
mat = np.asarray(f['/0/DATA/0/matrix'][:, :], dtype=np.float32)
f.close()

symbol2id = dict(zip(gene_info.pr_gene_symbol, gene_info.pr_gene_id))
id2col = {gid: i for i, gid in enumerate(gene_ids)}
syms = set(gene_info.pr_gene_symbol.dropna().astype(str))

sigs = json.load(open(f'{DATA}/disease_signatures.json'))


def map_genes(genes):
    ok, missing = [], []
    for g in genes:
        gid = symbol2id.get(g)
        if gid is not None and gid in id2col:
            ok.append(id2col[gid])
        else:
            missing.append(g)
    return np.array(sorted(set(ok)), dtype=np.int64), missing


print('=== signature genes on the L1000 gene axis ===')
mapped = {}
for dis, key in [('HCC', 'HCC'), ('CVD_HF', 'CVD_HF')]:
    for direction in ['up', 'down']:
        idx, miss = map_genes(sigs[key][direction])
        mapped[(dis, direction)] = idx
        print(f'{dis:7s} {direction:5s}: requested {len(sigs[key][direction]):3d}  mapped {len(idx):3d}  '
              f'missing {len(miss):2d}  {miss if miss else ""}')

hcc_up, hcc_dn = mapped[('HCC', 'up')], mapped[('HCC', 'down')]
hf_up, hf_dn = mapped[('CVD_HF', 'up')], mapped[('CVD_HF', 'down')]


def rev(col_up, col_dn, block=1 << 22):
    n = mat.shape[0]
    s = np.empty(n, dtype=np.float32)
    for i in range(0, n, block):
        j = min(i + block, n)
        c = mat[i:j]
        s[i:j] = c[:, col_dn].mean(1) - c[:, col_up].mean(1)
    return s


obs_hcc = rev(hcc_up, hcc_dn)
obs_hf = rev(hf_up, hf_dn)

prof = pd.DataFrame({'sig_id': prof_ids})
sg = sig_info.set_index('sig_id')
prof['pert_id'] = prof.sig_id.map(sg.pert_id).astype(str)
prof['pert_type'] = prof.sig_id.map(sg.pert_type)
pi_cp = pert_info[pert_info.pert_type == 'trt_cp'].copy()
pi_cp['has_name'] = pi_cp.pert_iname != pi_cp.pert_id
pi_cp = pi_cp.drop_duplicates('pert_id')
name_map = dict(zip(pi_cp.pert_id, pi_cp.pert_iname))
cid_map = dict(zip(pi_cp.pert_id, pi_cp.pubchem_cid))
prof['true_name'] = prof.pert_id.map(name_map)
prof['pubchem'] = prof.pert_id.map(cid_map)
prof['has_name'] = prof.pert_id.map(pi_cp.set_index('pert_id').has_name).fillna(False)

trt = prof[prof.pert_type == 'trt_cp'].copy()
trt = trt[trt.has_name & (trt.pubchem.astype(str) != '-666')].copy()
print(f'\n[universe] trt_cp profiles with a name and a PubChem CID = {len(trt)}, '
      f'compounds (by true_name) = {trt.true_name.nunique()}')

g = pd.factorize(trt.true_name)
trt['gcode'] = g[0]
names = np.array(g[1])
n_drug = len(names)

# contiguous grouping for fast per-group median
order = np.argsort(trt.gcode.values, kind='stable')
srt_codes = trt.gcode.values[order]
counts = np.bincount(srt_codes, minlength=n_drug)
offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])
lower_med_pos = offsets + (counts - 1) // 2          # lower median for even n
prof_pos = trt.index.values[order]                    # profile axis, grouped


def agg_mean(s):
    tot = np.bincount(srt_codes, weights=s[prof_pos], minlength=n_drug)
    c = np.where(counts == 0, 1, counts)
    return (tot / c).astype(np.float32)


def agg_median(s):
    ss = s[prof_pos]
    return ss[lower_med_pos].astype(np.float32)


out = {}
for agg_name, agg in [('mean', agg_mean), ('lower_median', agg_median)]:
    dh, df = agg(obs_hcc), agg(obs_hf)
    dual = int(((dh > 0) & (df > 0)).sum())
    out[agg_name] = (dh, df)
    print(f'[recompute] aggregate={agg_name:12s} compounds {n_drug}  '
          f'HCC>0 {int((dh>0).sum())}  HF>0 {int((df>0).sum())}  dual {dual} '
          f'({100*dual/n_drug:.1f}%)')

# ---- compare against the published intermediate files --------------------
prev = pd.read_csv(f'{DATA}/drug_scores_all3.csv')
prev_named = prev[prev.has_name & (prev.pubchem.astype(str) != '-666')]
print(f'\n[published] drug_scores_all3.csv named compounds with a CID = {len(prev_named)}, '
      f'dual-positive = {int(((prev_named.hcc_score>0)&(prev_named.cvd_score>0)).sum())}')
cand = pd.read_csv(f'{DATA}/candidate_drugs.csv')
print(f'[published] candidate_drugs.csv rows = {len(cand)}')

for agg_name, (dh, df) in out.items():
    m = pd.DataFrame({'compound': names, 'hcc': dh, 'hf': df}).merge(
        prev_named[['true_name', 'hcc_score', 'cvd_score']],
        left_on='compound', right_on='true_name', how='inner')
    if len(m) == 0:
        continue
    print(f'[agree] aggregate={agg_name:12s} shared compounds {len(m)}  '
          f'corr(HCC)={m.hcc.corr(m.hcc_score):.4f}  corr(HF)={m.hf.corr(m.cvd_score):.4f}  '
          f'HCC max abs difference={np.abs(m.hcc-m.hcc_score).max():.4f}')

# per-compound percentile of the manuscript's three cardiac glycosides
for agg_name, (dh, df) in out.items():
    tbl = pd.DataFrame({'compound': names, 'hcc': dh, 'hf': df})
    hcc_pct = tbl.hcc.rank(pct=True) * 100
    hf_pct = tbl.hf.rank(pct=True) * 100
    print(f'\n[cardiac glycosides] aggregate={agg_name}')
    for drug in ['peruvoside', 'sarmentogenin', 'strophanthidin', 'digoxin', 'digitoxin']:
        r = tbl[tbl.compound == drug]
        if len(r):
            i = r.index[0]
            print(f'  {drug:16s} HCC {r.hcc.iloc[0]:+.4f} ({hcc_pct[i]:5.1f}th pct)   '
                  f'HF {r.hf.iloc[0]:+.4f} ({hf_pct[i]:5.1f}th pct)   '
                  f'dual={bool(r.hcc.iloc[0]>0 and r.hf.iloc[0]>0)}')
        else:
            print(f'  {drug:16s} not in this compound set')

pd.DataFrame({'compound': names,
              'hcc_mean': out['mean'][0], 'hf_mean': out['mean'][1],
              'hcc_lower_median': out['lower_median'][0],
              'hf_lower_median': out['lower_median'][1],
              'n_sig': counts}).to_csv(f'{OUT}/A1b_recomputed_drug_scores.csv', index=False)
print(f'\nwrote {OUT}/A1b_recomputed_drug_scores.csv')
