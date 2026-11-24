import pandas as pd, numpy as np
from scipy import stats
cols = list(pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', nrows=0).columns)[1:]
meta = pd.read_csv('data/cmc_s5a_pseudobulk_meta.tsv', sep='\t')
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
x = pd.read_csv('/tmp/xist_row.tsv', sep='\t', header=None).values.ravel()
if len(x) == len(cols) + 1:
    x = x[1:]
t = pd.DataFrame({'col': cols, 'xist': x.astype(float)})
t['group'] = [c.split('|')[-1] for c in cols]
t['sample'] = [c.split('|')[0] for c in cols]
t = t.merge(meta[['col', 'n_cells']], on='col')
t = t[t.n_cells >= 30]
t['sex'] = t['sample'].map(sheet['Biological_Sex'])
rows = []
for ct, s in t.groupby('group'):
    a = s.loc[s.sex == 'male', 'xist'].values
    b = s.loc[s.sex == 'female', 'xist'].values
    if len(a) > 2 and len(b) > 2:
        r = stats.rankdata(np.r_[a, b])[len(a):]
        auc = (r.sum() - len(b) * (len(b) + 1) / 2) / (len(a) * len(b))
        rows.append((ct, len(a), len(b), np.median(a), np.median(b), round(auc, 3),
                     stats.mannwhitneyu(a, b).pvalue))
r = pd.DataFrame(rows, columns=['group', 'nM', 'nF', 'medM', 'medF', 'aucF', 'p'])
print(r.sort_values('aucF').to_string(index=False))
print('aucF 最小 %.3f  n p<0.05 %d/%d' % (r.aucF.min(), int((r.p < 0.05).sum()), len(r)))
r.to_csv('results_sex/cmc_sD2_xist_auc.csv', index=False)
d = pd.read_csv('results_sex/cmc_s8_deseq_all.csv')
d = d[d.group != 'L6 CT']
genes = ('RPS4Y1 DDX3Y UTY KDM5D USP9Y TTTY14 TTTY10 NLGN4Y PCDH11Y PRKY ZFY '
         'XIST JPX ZFX KDM6A PUDP DDX3X KDM5C PLCXD1 LINC00685 ASMTL-AS1 PPP2R3B CD99').split()
m = d[d.gene.isin(genes)].pivot_table(index='gene', columns='group', values='log2FoldChange')
print('矩阵', m.shape, '空格(未检验)', int(m.isna().sum().sum()))
m.to_csv('results_sex/cmc_sD3_lfc_matrix.csv')
