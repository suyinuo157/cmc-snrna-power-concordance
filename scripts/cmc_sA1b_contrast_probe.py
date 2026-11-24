import pandas as pd, numpy as np
pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
si = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
ct = 'L2/3 IT'; cols = [c for c in pb.columns if c.endswith('|' + ct)]
cpm = np.log2(pb[cols].div(pb[cols].sum(axis=0), axis=1) * 1e6 + 1)
s = [c.split('|')[0] for c in cols]
scz = np.array([si.condition[x] == 'SCZ' for x in s], float)
ph = np.array([pd.to_numeric(str(si.pH[x]).replace('+', ''), errors='coerce') for x in s])
ok = ~np.isnan(ph); X = cpm.values[:, ok]; y = ph[ok]
Xc = X - X.mean(1, keepdims=True); yc = y - y.mean()
r_ph = (Xc @ yc) / np.sqrt((Xc ** 2).sum(1) * (yc ** 2).sum())
lfc = X[:, scz[ok] == 1].mean(1) - X[:, scz[ok] == 0].mean(1)
d = pd.read_csv('results_qc/cmc_s5c_deseq_all.csv'); d = d[d.group == ct]
m = pd.DataFrame({'gene': cpm.index, 'r_ph': r_ph, 'lfc_scz': lfc}).merge(d, on='gene')
print('DESeq2 LFC vs SCZ contrast : %.3f' % np.corrcoef(m.lfc_scz, m.log2FoldChange)[0, 1])
print('DESeq2 LFC vs pH assoc     : %.3f' % np.corrcoef(m.r_ph, m.log2FoldChange)[0, 1])
print('pH assoc  vs SCZ contrast  : %.3f' % np.corrcoef(m.r_ph, m.lfc_scz)[0, 1])
