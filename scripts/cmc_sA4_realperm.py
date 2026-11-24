import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(7)
pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
meta = pd.read_csv('data/cmc_s5a_pseudobulk_meta.tsv', sep='\t')
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
qc = pd.read_csv('data/cmc_s5c_qc_by_group.tsv', sep='\t').set_index('col')
def num(x):
    return pd.to_numeric(pd.Series(x).astype(str).str.replace('+', '', regex=False), errors='coerce')
D = []
for g in meta.group.unique():
    m = meta[(meta.group == g) & (meta.n_cells >= 30)]
    if len(m) < 30:
        continue
    sub = pb[list(m['col'])]
    cpm = np.log2(sub.div(sub.sum(axis=0), axis=1) * 1e6 + 1)
    d = pd.DataFrame({'sample': list(m['sample'])}, index=list(m['col']))
    s = d['sample'].values
    d['SCZ'] = (sheet.loc[s, 'condition'].values == 'SCZ').astype(float)
    d['male'] = (sheet.loc[s, 'Biological_Sex'].str.lower().values == 'male').astype(float)
    for v, c in [('Age', 'Age_death'), ('PMI', 'PMI'), ('pH', 'pH')]:
        d[v] = num(sheet.loc[s, c].values).values
    d['qc'] = qc.reindex(list(m['col']))['median_n_genes'].values
    d = d.dropna()
    for v in ['Age', 'PMI', 'pH', 'qc']:
        d[v] = (d[v] - d[v].mean()) / d[v].std()
    D.append((g, cpm[d.index], d))
print('cell types used', len(D))
genes = None
for _, cpm, _ in D:
    s = set(cpm.index)
    genes = s if genes is None else genes & s
genes = sorted(genes)
E = np.mean([(cpm.loc[genes].mean(axis=1) > 1).values for _, cpm, _ in D], axis=0)
genes = [g for g, x in zip(genes, E) if x >= 0.5]
print('common expressed genes', len(genes))
def betas(cpm, d, lab):
    Y = cpm.loc[genes].values.T
    X = np.column_stack([np.ones(len(d)), lab, d.male.values, d.Age.values,
                         d.PMI.values, d.pH.values, d.qc.values])
    return (np.linalg.pinv(X) @ Y)[1]
def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(m); r[o] = np.clip(q, 0, 1); return r
def total_sig(blist):
    B = np.vstack(blist)
    up = (B > 0).sum(0); dn = (B < 0).sum(0)
    neff = up + dn; k = np.minimum(up, dn)
    p = np.ones(len(k))
    for nv in np.unique(neff):
        if nv < 8:
            continue
        sel = neff == nv
        p[sel] = np.minimum(2 * stats.binom.cdf(k[sel], nv, 0.5), 1)
    return int((bh(p) < 0.05).sum())
obs = total_sig([betas(cpm, d, d.SCZ.values) for _, cpm, d in D])
print('observed #sig (real labels):', obs)
lab = (sheet['condition'] == 'SCZ').astype(float)
null = []
for r in range(200):
    pl = pd.Series(rng.permutation(lab.values), index=lab.index)
    bl = [betas(cpm, d, pl.reindex(d['sample']).values) for _, cpm, d in D]
    null.append(total_sig(bl))
    if (r + 1) % 50 == 0:
        print('perm', r + 1, flush=True)
null = np.array(null)
print('null median %.1f  95th %.1f  max %d' % (np.median(null), np.percentile(null, 95), null.max()))
print('empirical p = %.3f' % float((null >= obs).mean()))
pd.DataFrame({'n_sig': null}).to_csv('results_sim/cmc_sA4_realperm_null.csv', index=False)
