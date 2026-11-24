import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

DATA = 'data'
BAD = {'32', '34', '38', '39', '40'}
GMAP = {'L2/3 IT':'L2/3 IT','L4 IT':'L4 IT','L5 IT':'L5 IT','L6 IT':'L6 IT','L6 IT Car3':'L6 IT',
        'L6 CT':'L6 CT','L6b':'L6b','L5/6 NP':'L5/6 NP','L5 ET':'L5/6 NP','Pvalb':'Pvalb','Sst':'Sst',
        'Sst Chodl':'Sst','Vip':'Vip','Sncg':'Sncg','Lamp5':'Lamp5','Lamp5 Lhx6':'Lamp5','Pax6':'Pax6',
        'Chandelier':'Chandelier','Astro':'Astro','Oligo':'Oligo','OPC':'OPC','Micro':'Myeloid',
        'Immune':'Myeloid','Endo':'Vascular','VLMC':'Vascular','SMC':'Vascular','PC':'Vascular'}

S = ad.read_h5ad(f'{DATA}/cmc_s2_clustered.h5ad', backed='r')
o = S.obs[['sample', 'celltype', 'leiden']].copy()
o.index = o.index.astype(str)
o = o[~o['leiden'].astype(str).isin(BAD)]
o['group'] = o['celltype'].astype(str).map(GMAP)
o = o.dropna(subset=['group'])
print('cells kept', o.shape[0], 'groups', o['group'].nunique(), flush=True)

A = ad.read_h5ad(f'{DATA}/cmc_all.h5ad', backed='r')
names = pd.Index(A.obs_names.astype(str))
o = o.reindex(names)
mask = o['group'].notna().to_numpy()
ikeys = pd.MultiIndex.from_frame(o.loc[mask, ['sample', 'group']])
combo = ikeys.unique()
proc = np.full(A.n_obs, -1, dtype=np.int64)
proc[mask] = combo.get_indexer(ikeys)
print('matched cells', int(mask.sum()), 'combos', len(combo), 'genes', A.n_vars, flush=True)

acc = np.zeros((len(combo), A.n_vars), dtype=np.float64)
CH = 10000
for i in range(0, A.n_obs, CH):
    j = min(i + CH, A.n_obs)
    pc = proc[i:j]
    m = pc >= 0
    if not m.any():
        continue
    Xi = A.X[i:j]
    Xi = sp.csr_matrix(Xi) if not sp.issparse(Xi) else Xi.tocsr()
    Xi = Xi[m]
    pc = pc[m]
    D = sp.csr_matrix((np.ones(pc.size, dtype=np.float32), (np.arange(pc.size), pc)),
                      shape=(pc.size, len(combo)))
    acc += np.asarray((D.T @ Xi).todense())
    print('chunk', j, '/', A.n_obs, flush=True)

df = pd.DataFrame(acc.T, index=pd.Index(A.var_names.astype(str), name='gene'))
df.columns = [f'{s}|{g}' for s, g in combo]
df.to_csv(f'{DATA}/cmc_s5a_pseudobulk.tsv', sep='\t')
meta = pd.DataFrame({'col': df.columns,
                     'sample': [s for s, _ in combo],
                     'group': [g for _, g in combo],
                     'n_cells': np.bincount(proc[proc >= 0], minlength=len(combo))})
meta.to_csv(f'{DATA}/cmc_s5a_pseudobulk_meta.tsv', sep='\t', index=False)
print('written', df.shape, 'combos >=30 cells', int((meta['n_cells'] >= 30).sum()), flush=True)
