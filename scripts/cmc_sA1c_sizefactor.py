import pandas as pd, numpy as np
pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
si = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
ct = 'L2/3 IT'; cols = [c for c in pb.columns if c.endswith('|' + ct)]
cnt = pb[cols].values.astype(float)
pos = (cnt > 0).all(axis=1)
lg = np.log(cnt[pos]); gm = lg.mean(axis=1)
sf = np.exp(np.median(lg - gm[:, None], axis=0))
norm = np.log2(cnt / sf + 1)
s = [c.split('|')[0] for c in cols]
scz = np.array([si.condition[x] == 'SCZ' for x in s], float)
lfc_sf = norm[:, scz == 1].mean(1) - norm[:, scz == 0].mean(1)
d = pd.read_csv('results_qc/cmc_s5c_deseq_all.csv'); d = d[d.group == ct]
m = pd.DataFrame({'gene': pb.index, 'lfc_sf': lfc_sf}).merge(d, on='gene')
print('genes with all-positive counts (SF basis):', int(pos.sum()))
print('size factor range: %.2f - %.2f' % (sf.min(), sf.max()))
print('LFC(size-factor) vs DESeq2 : %.3f' % np.corrcoef(m.lfc_sf, m.log2FoldChange)[0, 1])
print('spearman                   : %.3f' % pd.Series(m.lfc_sf).corr(
    pd.Series(m.log2FoldChange.values), method='spearman'))
mm = m.copy(); mm['mean'] = norm.mean(axis=1)[:len(mm)]
q = np.quantile(norm.mean(axis=1), [0.25, 0.5, 0.75])
