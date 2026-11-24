import pandas as pd, numpy as np
from scipy import stats
pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t')
qc = pd.read_csv('data/cmc_s5c_qc_by_group.tsv', sep='\t').set_index('col')
ct = 'L2/3 IT'
cols = [c for c in pb.columns if c.endswith('|' + ct)]
sub = pb[cols]
cpm = np.log2(sub.div(sub.sum(axis=0), axis=1) * 1e6 + 1)
si = sheet.set_index('Individual_ID')
def num(s):
    return pd.to_numeric(s.astype(str).str.replace('+', '', regex=False), errors='coerce')
df = pd.DataFrame({'col': cols})
df['sample'] = [c.split('|')[0] for c in cols]
df['SCZ'] = (df['sample'].map(si['condition']) == 'SCZ').astype(float)
df['male'] = (df['sample'].map(si['Biological_Sex']).str.lower() == 'male').astype(float)
df['Age'] = num(df['sample'].map(si['Age_death']))
df['PMI'] = num(df['sample'].map(si['PMI']))
df['pH'] = num(df['sample'].map(si['pH']))
df['qc'] = qc.reindex(cols)['median_n_genes'].values
print('before CC', len(df), 'missing', int(df[['Age','PMI','pH','qc']].isna().any(axis=1).sum()))
df = df.dropna(subset=['Age','PMI','pH','qc']).reset_index(drop=True)
print('after CC', len(df), 'SCZ', int(df.SCZ.sum()))
for c in ['Age','PMI','pH','qc']:
    df[c] = (df[c] - df[c].mean()) / df[c].std()
X = np.column_stack([np.ones(len(df)), df.SCZ, df.male, df.Age, df.PMI, df.pH, df.qc])
Y = cpm[df.col].values.T
beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
resid = Y - X @ beta
dof = len(df) - X.shape[1]
s2 = (resid ** 2).sum(axis=0) / dof
inv = np.linalg.inv(X.T @ X)
se = np.sqrt(s2 * inv[1, 1])
tstat = beta[1] / se
p = 2 * stats.t.sf(np.abs(tstat), dof)
ols = pd.DataFrame({'gene': cpm.index, 'lfc_ols': beta[1], 'p_ols': p})
d = pd.read_csv('results_qc/cmc_s5c_deseq_all.csv')
d = d[d.group == ct][['gene', 'log2FoldChange', 'pvalue']]
m = ols.merge(d, on='gene').dropna()
m['mean'] = cpm.mean(axis=1).reindex(m.gene).values
print('n', len(m))
print('pearson  %.4f' % np.corrcoef(m.lfc_ols, m.log2FoldChange)[0, 1])
print('spearman %.4f' % stats.spearmanr(m.lfc_ols, m.log2FoldChange).statistic)
print('median|diff| %.4f' % (m.lfc_ols - m.log2FoldChange).abs().median())
print('plog10 rank corr %.4f' % stats.spearmanr(-np.log10(m.p_ols + 1e-300),
                                                 -np.log10(m.pvalue + 1e-300)).statistic)
q = m['mean'].quantile([0.25, 0.5, 0.75]).values
for lo, hi, lab in [(-1e9, q[0], 'Q1'), (q[0], q[1], 'Q2'), (q[1], q[2], 'Q3'), (q[2], 1e9, 'Q4')]:
    s = m[(m['mean'] >= lo) & (m['mean'] < hi)]
    if len(s) > 100:
        print('%-3s n=%5d pearson %.3f' % (lab, len(s),
              np.corrcoef(s.lfc_ols, s.log2FoldChange)[0, 1]))
ols = pd.DataFrame({'gene': cpm.index, 'lfc_ols': beta[1], 'p_ols': p})
d = pd.read_csv('results_qc/cmc_s5c_deseq_all.csv')
d = d[d.group == ct][['gene', 'log2FoldChange', 'pvalue']]
m = ols.merge(d, on='gene').dropna()
m['mean'] = cpm.mean(axis=1).reindex(m.gene).values
print('n', len(m))
print('pearson  %.4f' % np.corrcoef(m.lfc_ols, m.log2FoldChange)[0, 1])
print('spearman %.4f' % stats.spearmanr(m.lfc_ols, m.log2FoldChange).statistic)
print('median|diff| %.4f' % (m.lfc_ols - m.log2FoldChange).abs().median())
print('plog10 rank corr %.4f' % stats.spearmanr(
    -np.log10(m.p_ols + 1e-300), -np.log10(m.pvalue + 1e-300)).statistic)
q = m['mean'].quantile([0.25, 0.5, 0.75]).values
for lo, hi, lab in [(-1e9, q[0], 'Q1'), (q[0], q[1], 'Q2'),
                    (q[1], q[2], 'Q3'), (q[2], 1e9, 'Q4')]:
    s = m[(m['mean'] >= lo) & (m['mean'] < hi)]
    if len(s) > 100:
        print('%-3s n=%5d pearson %.3f' % (lab, len(s),
              np.corrcoef(s.lfc_ols, s.log2FoldChange)[0, 1]))
