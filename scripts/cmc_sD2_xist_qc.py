import pandas as pd, numpy as np
cols = list(pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', nrows=0).columns)[1:]
meta = pd.read_csv('data/cmc_s5a_pseudobulk_meta.tsv', sep='\t')
x = pd.read_csv('/tmp/xist_row.tsv', sep='\t', header=None).values.ravel()
if len(x) == len(cols) + 1:
    x = x[1:]
print('ncol', len(cols), 'nval', len(x), 'nmeta', len(meta))
print('cols[:2]', cols[:2]); print(meta.head(2).to_string())
g = np.array([c.split('|')[-1] for c in cols])
ok = len(g) == len(meta) and (g == meta['group'].values).all()
print('align ok:', ok)
if ok:
    t = pd.DataFrame({'group': g, 'n_cells': meta['n_cells'].values, 'xist': x})
    t = t[t.n_cells >= 30]
    r = t.groupby('group')['xist'].agg(mean='mean', median='median', mx='max', n='size')
    print(r.round(2).to_string())
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
qc = pd.read_csv('data/cmc_s5c_qc_by_group.tsv', sep='\t')[['col', 'median_n_genes']]
q = qc.merge(meta[['col', 'sample', 'group', 'n_cells']], on='col')
q = q[q.n_cells >= 30]
q['sex'] = q['sample'].map(sheet['Biological_Sex'])
print('\nsex-qc 全局:')
print(q.groupby('sex')['median_n_genes'].describe().round(1).to_string())
from scipy import stats
res = []
for ct, sub in q.groupby('group'):
    a = sub.loc[sub.sex == 'male', 'median_n_genes']
    b = sub.loc[sub.sex == 'female', 'median_n_genes']
    if len(a) > 2 and len(b) > 2:
        res.append((ct, len(a), len(b), round(a.median(), 1), round(b.median(), 1),
                    stats.mannwhitneyu(a, b).pvalue))
rr = pd.DataFrame(res, columns=['group', 'nM', 'nF', 'medM', 'medF', 'p'])
rr['fdr'] = stats.false_discovery_control(rr.p)
print(rr.sort_values('p').to_string(index=False))
print('n p<0.05:', int((rr.p < 0.05).sum()))
