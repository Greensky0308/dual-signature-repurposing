#!/usr/bin/env python3
"""A1: permutation null for the dual-signature reversal screen.

Reproduces the manuscript's scoring pipeline exactly (per-profile reversal score
-> per-compound median across all trt_cp profiles -> compound-level annotation
filter), then permutes the disease-signature gene labels within the L1000
measured-gene universe to build null distributions for

  (a) the number of dual-positive compounds (the manuscript's 867 criterion),
  (b) per-compound reversal scores, giving empirical p-values for candidates.

Usage: python3 A1_permutation_null.py [N_PERM]
"""
import json, sys, time
import numpy as np
import pandas as pd
import h5py

from config import OUT, RAW as DATA, SEED
N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

t0 = time.time()
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
print(f'[load] matrix in RAM {mat.shape} {mat.nbytes/1e9:.1f} GB ({time.time()-t0:.0f}s)', flush=True)

symbol2id = dict(zip(gene_info.pr_gene_symbol, gene_info.pr_gene_id))
id2col = {gid: i for i, gid in enumerate(gene_ids)}
universe = np.array(sorted({id2col[gid] for gid in symbol2id.values() if gid in id2col}),
                    dtype=np.int64)
print(f'[pool] measured genes with a symbol: {len(universe)}', flush=True)

sigs = json.load(open(f'{DATA}/disease_signatures.json'))


def map_genes(genes):
    return np.array(sorted({id2col[symbol2id[g]] for g in genes
                            if symbol2id.get(g) in id2col}), dtype=np.int64)


hcc_up, hcc_dn = map_genes(sigs['HCC']['up']), map_genes(sigs['HCC']['down'])
hf_up, hf_dn = map_genes(sigs['CVD_HF']['up']), map_genes(sigs['CVD_HF']['down'])
print(f'[sig] HCC up/down {len(hcc_up)}/{len(hcc_dn)}  HF up/down {len(hf_up)}/{len(hf_dn)}'
      f'  (nominal {len(sigs["HCC"]["up"])}/{len(sigs["HCC"]["down"])},'
      f' {len(sigs["CVD_HF"]["up"])}/{len(sigs["CVD_HF"]["down"])})', flush=True)

# ---- compound grouping, exactly as the manuscript ------------------------
sg = sig_info.set_index('sig_id')
prof = pd.DataFrame({'sig_id': prof_ids})
prof['pert_id'] = prof.sig_id.map(sg.pert_id).astype(str)
prof['pert_type'] = prof.sig_id.map(sg.pert_type)
pi_cp = pert_info[pert_info.pert_type == 'trt_cp'].copy()
pi_cp['has_name'] = pi_cp.pert_iname != pi_cp.pert_id
nm = dict(zip(pi_cp.pert_id, pi_cp.pert_iname))
prof['true_name'] = prof.pert_id.map(nm)

trt = prof[prof.pert_type == 'trt_cp'].copy()
trt['pubchem'] = trt.pert_id.map(dict(zip(pi_cp.pert_id, pi_cp.pubchem_cid)))
trt['has_name'] = trt.pert_id.map(dict(zip(pi_cp.pert_id, pi_cp.has_name)))
prof_pos = trt.index.values.astype(np.int64)                    # positions into the profile axis

# median aggregator over the full trt_cp profile set
codes, names = pd.factorize(trt.true_name)
n_grp = len(names)
counts = np.bincount(codes, minlength=n_grp)
offsets = np.concatenate([[0], np.cumsum(counts)[:-1]])
lo = offsets + (counts - 1) // 2
hi = offsets + counts // 2
even = (counts % 2) == 0
even_lo = lo[even]
even_hi = hi[even]


def group_median(s):
    """s is indexed in trt_cp profile order (rows of mat_t).

    Values are sorted within each compound group (lexsort by group code, then
    by score) so that the middle rank of each group is its median.
    """
    idx = np.lexsort((s, codes))
    ss = s[idx]
    med = ss[lo].astype(np.float64)
    if even.any():
        med[even] = (ss[even_lo].astype(np.float64) + ss[even_hi].astype(np.float64)) / 2.0
    return med.astype(np.float32)


# compound-level annotation filter (first value per group, as in the original)
grp_tbl = trt.groupby('true_name').agg(has_name=('has_name', 'first'),
                                       pubchem=('pubchem', 'first')).loc[names]
keep_grp = (grp_tbl.has_name.values.astype(bool) &
            (grp_tbl.pubchem.values.astype(str) != '-666'))
print(f'[annot] groups {n_grp}  passing compound-level annotation filter {keep_grp.sum()}',
      flush=True)

mat_t = mat[prof_pos]                              # trt_cp profiles only
matT = np.ascontiguousarray(mat_t.T)               # genes x profiles: row select is contiguous
del mat, mat_t
print(f'[load] trt_cp transposed {matT.shape} {matT.nbytes/1e9:.1f} GB ({time.time()-t0:.0f}s)',
      flush=True)


def rev(col_up, col_dn):
    return matT[col_dn].mean(0) - matT[col_up].mean(0)


def dual_count(dh, df):
    return int((keep_grp & (dh > 0) & (df > 0)).sum())


print('[obs] computing observed scores', flush=True)
obs_hcc, obs_hf = rev(hcc_up, hcc_dn), rev(hf_up, hf_dn)
o_hcc, o_hf = group_median(obs_hcc), group_median(obs_hf)
obs_dual = dual_count(o_hcc, o_hf)
obs_p_hcc = float((obs_hcc > 0).mean())
obs_p_hf = float((obs_hf > 0).mean())
print(f'[obs] profile-level P(score>0): HCC {obs_p_hcc:.4f}  HF {obs_p_hf:.4f}  '
      f'product {obs_p_hcc*obs_p_hf:.4f}', flush=True)
print(f'[obs] dual-positive compounds (manuscript criterion) = {obs_dual}  '
      f'[manuscript reports 867]', flush=True)

# ---- permutation ---------------------------------------------------------
rng = np.random.default_rng(SEED)
sizes = (len(hcc_up), len(hcc_dn), len(hf_up), len(hf_dn))
n_pick = sum(sizes)
null_dual = np.zeros(N_PERM, dtype=np.int32)
null_p_hcc = np.zeros(N_PERM, dtype=np.float32)
null_p_hf = np.zeros(N_PERM, dtype=np.float32)
null_g_hcc = np.zeros((N_PERM, keep_grp.sum()), dtype=np.float32)
null_g_hf = np.zeros((N_PERM, keep_grp.sum()), dtype=np.float32)

tp = time.time()
for p in range(N_PERM):
    pick = rng.choice(len(universe), size=n_pick, replace=False)
    a = 0
    pu = universe[pick[a:a + sizes[0]]]; a += sizes[0]
    qd = universe[pick[a:a + sizes[1]]]; a += sizes[1]
    ru = universe[pick[a:a + sizes[2]]]; a += sizes[2]
    rd = universe[pick[a:a + sizes[3]]]
    sh, sf = rev(pu, qd), rev(ru, rd)
    gh, gf = group_median(sh), group_median(sf)
    null_g_hcc[p] = gh[keep_grp]
    null_g_hf[p] = gf[keep_grp]
    null_dual[p] = dual_count(gh, gf)
    null_p_hcc[p] = (sh > 0).mean()
    null_p_hf[p] = (sf > 0).mean()
    if N_PERM <= 5 or (p + 1) % 100 == 0:
        el = time.time() - tp
        print(f'  perm {p+1}/{N_PERM}  {el/(p+1):.2f}s/perm  eta {(el/(p+1))*(N_PERM-p-1)/60:.1f}min'
              f'  null_dual run-median {np.median(null_dual[:p+1]):.0f}', flush=True)

# ---- summary -------------------------------------------------------------
def emp_p(null, obs):
    return float((1 + int((null >= obs).sum())) / (1 + len(null)))


summary = {
    'n_perm': N_PERM, 'seed': SEED,
    'gene_universe': int(len(universe)),
    'signature_sizes_effective': {'hcc_up': sizes[0], 'hcc_down': sizes[1],
                                  'hf_up': sizes[2], 'hf_down': sizes[3]},
    'n_trt_cp_profiles': int(matT.shape[1]),
    'n_compound_groups': int(n_grp),
    'n_compounds_annotated': int(keep_grp.sum()),
    'observed': {
        'profile_p_hcc_pos': obs_p_hcc, 'profile_p_hf_pos': obs_p_hf,
        'profile_independence_product': obs_p_hcc * obs_p_hf,
        'dual_compounds': obs_dual,
        'dual_rate': obs_dual / int(keep_grp.sum()),
        'reproduced_manuscript_867': obs_dual == 867,
    },
    'null_dual': {
        'mean': float(null_dual.mean()), 'sd': float(null_dual.std(ddof=1)),
        'median': float(np.median(null_dual)),
        'q025': float(np.percentile(null_dual, 2.5)),
        'q975': float(np.percentile(null_dual, 97.5)),
        'min': int(null_dual.min()), 'max': int(null_dual.max()),
    },
    'p_dual_count': emp_p(null_dual, obs_dual),
    'enrichment_ratio_vs_null': float(obs_dual / null_dual.mean()),
    'dual_rate_null_mean': float(null_dual.mean() / int(keep_grp.sum())),
    'null_profile_p': {
        'hcc_mean': float(null_p_hcc.mean()), 'hcc_sd': float(null_p_hcc.std(ddof=1)),
        'hf_mean': float(null_p_hf.mean()), 'hf_sd': float(null_p_hf.std(ddof=1)),
    },
}

kept_names = names[keep_grp]
o_h = o_hcc[keep_grp]; o_f = o_hf[keep_grp]

# p-value uniformity test over the whole annotated compound set: under the null
# each empirical p-value is uniform, so 5% should fall below 0.05 per arm.
ph_all = (1 + (null_g_hcc >= o_h).sum(0)) / (1 + N_PERM)
pf_all = (1 + (null_g_hf >= o_f).sum(0)) / (1 + N_PERM)
n_ann = int(keep_grp.sum())
summary['arm_pvalue_calibration'] = {
    'n_annotated': n_ann,
    'hcc_p_lt_0.05': int((ph_all < 0.05).sum()),
    'hf_p_lt_0.05': int((pf_all < 0.05).sum()),
    'both_p_lt_0.05': int(((ph_all < 0.05) & (pf_all < 0.05)).sum()),
    'expected_uniform_each_arm': 0.05 * n_ann,
    'expected_uniform_both_if_independent': 0.0025 * n_ann,
}
np.save(f'{OUT}/A1_all_pvals_hcc.npy', ph_all)
np.save(f'{OUT}/A1_all_pvals_hf.npy', pf_all)

cand = np.flatnonzero((o_h > 0) & (o_f > 0))
p_hcc_c = (1 + (null_g_hcc[:, cand] >= o_h[cand]).sum(0)) / (1 + N_PERM)
p_hf_c = (1 + (null_g_hf[:, cand] >= o_f[cand]).sum(0)) / (1 + N_PERM)
pd.DataFrame({
    'compound': kept_names[cand],
    'hcc_score': o_h[cand], 'hf_score': o_f[cand],
    'hcc_p_perm': p_hcc_c, 'hf_p_perm': p_hf_c,
}).sort_values('hcc_score', ascending=False).to_csv(
    f'{OUT}/A1_candidate_perm_pvalues.csv', index=False)
pd.DataFrame({'compound': kept_names, 'hcc_score': o_h, 'hf_score': o_f}).to_csv(
    f'{OUT}/A1_observed_drug_scores.csv', index=False)
np.save(f'{OUT}/A1_null_dual.npy', null_dual)
with open(f'{OUT}/A1_summary.json', 'w') as fh:
    json.dump(summary, fh, indent=2)

print('\n=== SUMMARY ===')
print(json.dumps(summary, indent=2))
print(f'\n[total] {time.time()-t0:.0f}s')
