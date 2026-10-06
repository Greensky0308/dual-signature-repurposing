#!/usr/bin/env python3
"""A3: does the 300-vs-48 gene imbalance between the HCC and HF signatures
affect reversal scores and compound rankings?

Three experiments, all on the manuscript's exact scoring/aggregation pipeline:
  1. size-matched HCC: subsample the HCC signature down to the HF effective size
     (21 up / 19 down) and ask whether the dual-reversal candidate set and the
     compound ranking survive;
  2. noise curve: score dispersion as a function of signature size;
  3. split-half reliability of each signature (full size and size-matched).
"""
import json
import numpy as np
import pandas as pd
import h5py
from scipy import stats

from config import OUT, RAW as DATA, SEED
rng = np.random.default_rng(SEED)

gi = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_gene_info.txt.gz', sep='\t', compression='gzip')
gi['pr_gene_id'] = gi.pr_gene_id.astype(str)
si = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_sig_info.txt.gz', sep='\t',
                 compression='gzip', low_memory=False)
si['sig_id'] = si.sig_id.astype(str)
pi = pd.read_csv(f'{DATA}/GSE92742_Broad_LINCS_pert_info.txt.gz', sep='\t',
                 compression='gzip', low_memory=False)

f = h5py.File(f'{DATA}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx', 'r')
gene_ids = [x.decode() for x in f['/0/META/ROW/id'][:]]
prof_ids = [x.decode() for x in f['/0/META/COL/id'][:]]
mat = np.asarray(f['/0/DATA/0/matrix'][:, :], dtype=np.float32)
f.close()

s2i = dict(zip(gi.pr_gene_symbol, gi.pr_gene_id))
i2c = {g: i for i, g in enumerate(gene_ids)}


def mg(genes):
    return np.array(sorted({i2c[s2i[g]] for g in genes if s2i.get(g) in i2c}), dtype=np.int64)


S = json.load(open(f'{DATA}/disease_signatures.json'))
hcc_up, hcc_dn = mg(S['HCC']['up']), mg(S['HCC']['down'])
hf_up, hf_dn = mg(S['CVD_HF']['up']), mg(S['CVD_HF']['down'])
print(f'HCC {len(hcc_up)}up/{len(hcc_dn)}down   HF {len(hf_up)}up/{len(hf_dn)}down')

sg = si.set_index('sig_id')
prof = pd.DataFrame({'sig_id': prof_ids})
prof['pert_id'] = prof.sig_id.map(sg.pert_id).astype(str)
prof['pert_type'] = prof.sig_id.map(sg.pert_type)
pi_cp = pi[pi.pert_type == 'trt_cp'].copy()
pi_cp['has_name'] = pi_cp.pert_iname != pi_cp.pert_id
prof['true_name'] = prof.pert_id.map(dict(zip(pi_cp.pert_id, pi_cp.pert_iname)))
trt = prof[prof.pert_type == 'trt_cp'].copy()
trt['pubchem'] = trt.pert_id.map(dict(zip(pi_cp.pert_id, pi_cp.pubchem_cid)))
trt['has_name'] = trt.pert_id.map(dict(zip(pi_cp.pert_id, pi_cp.has_name)))

codes, names = pd.factorize(trt.true_name)
n_grp = len(names)
counts = np.bincount(codes, minlength=n_grp)
offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])
lo = offsets + (counts - 1) // 2
hi = offsets + counts // 2
even = (counts % 2) == 0

grp_tbl = trt.groupby('true_name').agg(has_name=('has_name', 'first'),
                                       pubchem=('pubchem', 'first')).loc[names]
keep_grp = (grp_tbl.has_name.values.astype(bool) &
            (grp_tbl.pubchem.values.astype(str) != '-666'))

mat_t = mat[trt.index.values.astype(np.int64)]
matT = np.ascontiguousarray(mat_t.T)          # genes x profiles
del mat, mat_t


def score(up, dn):
    return matT[dn].mean(0) - matT[up].mean(0)


def drug_med(s):
    idx = np.lexsort((s, codes))
    ss = s[idx]
    med = ss[lo].astype(np.float64)
    if even.any():
        med[even] = (ss[lo][even].astype(np.float64) + ss[hi][even].astype(np.float64)) / 2.0
    return med.astype(np.float32)


o_hcc, o_hf = score(hcc_up, hcc_dn), score(hf_up, hf_dn)
g_hcc, g_hf = drug_med(o_hcc), drug_med(o_hf)
base_dual = int((keep_grp & (g_hcc > 0) & (g_hf > 0)).sum())
print(f'baseline: dual-positive {base_dual} / {keep_grp.sum()}  '
      f'sd(HCC score) {o_hcc.std():.4f}  sd(HF score) {o_hf.std():.4f}  '
      f'ratio {o_hcc.std()/o_hf.std():.2f}')

res = {'baseline_dual': base_dual, 'n_annotated': int(keep_grp.sum()),
       'sd_hcc_profile': float(o_hcc.std()), 'sd_hf_profile': float(o_hf.std())}

# ---- experiment 1: size-matched HCC --------------------------------------
k_up, k_dn = len(hf_up), len(hf_dn)          # 21 / 19
N_REP = 200
dual_counts, spear, cg_ranks = [], [], []
base_rank = pd.Series(g_hcc[keep_grp]).rank(ascending=False).values
base_names = names[keep_grp]
for r in range(N_REP):
    su = rng.choice(hcc_up, size=k_up, replace=False)
    sd = rng.choice(hcc_dn, size=k_dn, replace=False)
    g = drug_med(score(su, sd))[keep_grp]
    dual_counts.append(int(((g > 0) & (g_hf[keep_grp] > 0)).sum()))
    spear.append(stats.spearmanr(g, g_hcc[keep_grp]).statistic)
    cg_ranks.append(pd.Series(g).rank(ascending=False).values)
dual_counts = np.array(dual_counts)
spear = np.array(spear)
res['size_matched_hcc'] = {
    'k_up': int(k_up), 'k_dn': int(k_dn), 'n_rep': N_REP,
    'dual_mean': float(dual_counts.mean()), 'dual_sd': float(dual_counts.std(ddof=1)),
    'dual_q025': float(np.percentile(dual_counts, 2.5)),
    'dual_q975': float(np.percentile(dual_counts, 97.5)),
    'spearman_vs_full_hcc_mean': float(spear.mean()),
    'spearman_vs_full_hcc_q025': float(np.percentile(spear, 2.5)),
}
print(f'\n[exp1] HCC down-sampled to {k_up}up/{k_dn}down ({N_REP} replicates): dual-positive '
      f'{dual_counts.mean():.0f} ± {dual_counts.std(ddof=1):.0f} '
      f'[{np.percentile(dual_counts,2.5):.0f}, {np.percentile(dual_counts,97.5):.0f}]'
      f'  (baseline {base_dual})')
print(f'       HCC ranking Spearman vs full signature: {spear.mean():.3f} '
      f'[{np.percentile(spear,2.5):.3f}, {np.percentile(spear,97.5):.3f}]')

# cardiac glycoside ranks under size-matched HCC
for drug in ['peruvoside', 'sarmentogenin', 'strophanthidin', 'digoxin', 'digitoxin']:
    j = np.flatnonzero(base_names == drug)
    if len(j):
        j = j[0]
        rk = [cg_ranks[i][j] for i in range(N_REP)]
        res.setdefault('size_matched_cg_ranks', {})[drug] = {
            'base_rank': float(base_rank[j]), 'matched_median_rank': float(np.median(rk)),
            'matched_q025': float(np.percentile(rk, 2.5)),
            'matched_q975': float(np.percentile(rk, 97.5))}
        print(f'       {drug:16s} full-signature HCC rank {int(base_rank[j]):5d}   '
              f'down-sampled median {np.median(rk):7.0f} [{np.percentile(rk,2.5):.0f}, {np.percentile(rk,97.5):.0f}]')

# ---- experiment 2: noise curve -------------------------------------------
curve = {}
for k in [10, 20, 40, 80, 150]:
    ku, kd = k // 2, k - k // 2
    if ku > len(hcc_up):
        continue
    sc = np.empty((100, matT.shape[1]), dtype=np.float32)
    for r in range(100):
        sc[r] = score(rng.choice(hcc_up, ku, False), rng.choice(hcc_dn, kd, False))
    # inter-replicate sd relative to the observed score sd across compounds
    rep_sd = sc.std(0, ddof=1)
    curve[k] = {'mean_inter_rep_sd': float(rep_sd.mean()),
                'observed_sd': float(o_hcc.std()),
                'noise_to_signal': float(rep_sd.mean() / o_hcc.std())}
    print(f'[exp2] HCC signature k={k:3d}: inter-replicate sd {rep_sd.mean():.4f}  '
          f'observed sd {o_hcc.std():.4f}  noise/signal {rep_sd.mean()/o_hcc.std():.3f}')
res['noise_curve'] = curve

# HF reference point at its own size, same 100-replicate protocol
sc = np.empty((100, matT.shape[1]), dtype=np.float32)
for r in range(100):
    sc[r] = score(rng.choice(hf_up, len(hf_up) // 2, False),
                  rng.choice(hf_dn, len(hf_dn) // 2, False))
rep_sd = sc.std(0, ddof=1)
res['noise_hf_halfsize'] = {'k': int(len(hf_up) // 2 + len(hf_dn) // 2),
                            'mean_inter_rep_sd': float(rep_sd.mean()),
                            'observed_sd': float(o_hf.std()),
                            'noise_to_signal': float(rep_sd.mean() / o_hf.std())}
print(f'[exp2] HF half-signature k={len(hf_up)//2+len(hf_dn)//2}: inter-replicate sd {rep_sd.mean():.4f}  '
      f'observed sd {o_hf.std():.4f}  noise/signal {rep_sd.mean()/o_hf.std():.3f}')

# ---- experiment 3: split-half reliability --------------------------------
def split_half(up, dn, reps=200):
    cs = []
    for _ in range(reps):
        pu = rng.permutation(up); pd_ = rng.permutation(dn)
        a = score(pu[:len(pu) // 2], pd_[:len(pd_) // 2])
        b = score(pu[len(pu) // 2:], pd_[len(pd_) // 2:])
        cs.append(stats.spearmanr(a, b).statistic)
    cs = np.array(cs)
    # Spearman-Brown correction to full-length reliability
    sb = 2 * cs / (1 + cs)
    return float(np.median(cs)), float(np.median(sb))


for lab, up, dn in [('HCC_full', hcc_up, hcc_dn), ('HF_full', hf_up, hf_dn)]:
    r_split, r_sb = split_half(up, dn)
    res[f'split_half_{lab}'] = {'median_r': r_split, 'spearman_brown': r_sb,
                                'n_up': len(up), 'n_down': len(dn)}
    print(f'[exp3] {lab:9s} ({len(up)}up/{len(dn)}down) split-half median r={r_split:.3f}  '
          f'Spearman-Brown corrected {r_sb:.3f}')

# size-matched comparison: HCC subsampled to 40 genes, split-half
cs = []
for _ in range(200):
    su = rng.choice(hcc_up, 20, False); sd = rng.choice(hcc_dn, 20, False)
    a = score(su[:10], sd[:10]); b = score(su[10:], sd[10:])
    cs.append(stats.spearmanr(a, b).statistic)
cs = np.array(cs)
res['split_half_HCC_sizematched_40'] = {'median_r': float(np.median(cs)),
                                        'spearman_brown': float(np.median(2 * cs / (1 + cs)))}
print(f'[exp3] HCC_40  (20up/20down) split-half median r={np.median(cs):.3f}  '
      f'Spearman-Brown corrected {np.median(2*cs/(1+cs)):.3f}')

with open(f'{OUT}/A3_size_matched.json', 'w') as fh:
    json.dump(res, fh, indent=2, default=float)
print(f'\nwrote {OUT}/A3_size_matched.json')
