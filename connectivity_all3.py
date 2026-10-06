import h5py, json

import config
from config import OUT, RAW
import numpy as np
import pandas as pd

gene_info = pd.read_csv(f'{RAW}/GSE92742_Broad_LINCS_gene_info.txt.gz', sep='\t', compression='gzip')
sig_info = pd.read_csv(f'{RAW}/GSE92742_Broad_LINCS_sig_info.txt.gz', sep='\t', compression='gzip', low_memory=False)
pi = pd.read_csv(f'{RAW}/GSE92742_Broad_LINCS_pert_info.txt.gz', sep='\t', compression='gzip', low_memory=False)
gene_info['pr_gene_id'] = gene_info['pr_gene_id'].astype(str)

f = h5py.File(config.L1000, 'r')
row_ids = [x.decode() for x in f['/0/META/ROW/id'][:]]
col_ids = [x.decode() for x in f['/0/META/COL/id'][:]]

symbol2geneid = dict(zip(gene_info.pr_gene_symbol, gene_info.pr_gene_id))
geneid2row = {gid: i for i, gid in enumerate(row_ids)}
signatures = json.load(open(f'{RAW}/disease_signatures.json'))

def gene_col_indices(genes):
    idx = []
    for g in genes:
        gid = symbol2geneid.get(g)
        if gid is not None and gid in geneid2row:
            idx.append(geneid2row[gid])
    return idx

diseases = ['HCC', 'CVD_HF', 'CAD']
all_cols = set()
col_groups = {}
for d in diseases:
    up = gene_col_indices(signatures[d]['up'])
    down = gene_col_indices(signatures[d]['down'])
    col_groups[d] = (up, down)
    all_cols.update(up + down)
    print(f'{d}: up {len(up)} / down {len(down)} gene columns')

all_cols = sorted(all_cols)
print(f'total gene columns: {len(all_cols)}')

mat = np.asarray(f['/0/DATA/0/matrix'][:, all_cols], dtype=np.float32)
f.close()
col_map = {c: i for i, c in enumerate(all_cols)}

def rev_score(v, up, down):
    up_mean = v[[col_map[c] for c in up]].mean() if up else 0
    down_mean = v[[col_map[c] for c in down]].mean() if down else 0
    return down_mean - up_mean

scores = np.zeros((len(col_ids), len(diseases)), dtype=np.float32)
for i in range(len(col_ids)):
    v = mat[i]
    for j, d in enumerate(diseases):
        up, down = col_groups[d]
        scores[i, j] = rev_score(v, up, down)

res = pd.DataFrame({'sig_id': col_ids, 'hcc_score': scores[:,0], 'cvd_score': scores[:,1], 'cad_score': scores[:,2]})
sig2 = sig_info.set_index(sig_info.sig_id.astype(str))
res['pert_id'] = res.sig_id.astype(str).map(sig2.pert_id)
res['pert_type'] = res.sig_id.astype(str).map(sig2.pert_type)
res['cell_id'] = res.sig_id.astype(str).map(sig2.cell_id)
res.to_csv(f'{RAW}/connectivity_all3.csv', index=False)
print(f'wrote {RAW}/connectivity_all3.csv')

# aggregate to compound level
pi_cp = pi[pi.pert_type=='trt_cp'].copy()
pi_cp['has_name'] = pi_cp.pert_iname != pi_cp.pert_id
name_map = dict(zip(pi_cp.pert_id, pi_cp.pert_iname))
cid_map = dict(zip(pi_cp.pert_id, pi_cp.pubchem_cid))
res['true_name'] = res.pert_id.map(name_map)
res['has_name'] = res.pert_id.map(pi_cp.set_index('pert_id').has_name).fillna(False)
res['pubchem'] = res.pert_id.map(cid_map)

trt = res[res.pert_type=='trt_cp'].copy()
drug = trt.groupby('true_name').agg(
    hcc_score=('hcc_score','median'), cvd_score=('cvd_score','median'),
    cad_score=('cad_score','median'), n_sig=('sig_id','count'),
    has_name=('has_name','first'), pubchem=('pubchem','first')
).reset_index()
drug['joint_cvd'] = drug.hcc_score + drug.cvd_score
drug['joint_cad'] = drug.hcc_score + drug.cad_score
drug.to_csv(f'{RAW}/drug_scores_all3.csv', index=False)

# candidates: HCC + heart failure, and HCC + CAD
for combo, score_col, joint_col in [('HCC+HF','cvd_score','joint_cvd'), ('HCC+CAD','cad_score','joint_cad')]:
    cand = drug[(drug.hcc_score>0)&(drug[score_col]>0)&(drug.has_name)&(drug.pubchem!='-666')].sort_values(joint_col, ascending=False)
    cand.to_csv(f'{RAW}/candidate_drugs_{combo}.csv', index=False)
    print(f'\n{combo} candidates: {len(cand)}, top10:')
    print(cand[['true_name','hcc_score',score_col,joint_col]].head(10).to_string(index=False))
