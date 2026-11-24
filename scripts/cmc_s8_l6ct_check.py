import pandas as pd
d = pd.read_csv('results_sex/cmc_s8_deseq_all.csv')
l6 = d[d.group == 'L6 CT']
print('--- L6 CT 性染色体基因 ---')
print(l6[l6.gene.isin(['XIST', 'RPS4Y1', 'UTY', 'DDX3Y', 'KDM5D'])][
    ['gene', 'log2FoldChange', 'pvalue', 'padj', 'baseMean']].to_string(index=False))
s = l6[l6.padj < 0.05]
print('\nsig', len(s), 'up', int((s.log2FoldChange > 0).sum()), 'down', int((s.log2FoldChange < 0).sum()))
print('sig baseMean 分位:', s.baseMean.quantile([.1, .5, .9]).round(1).tolist())
print('全部 baseMean 分位:', l6.baseMean.quantile([.1, .5, .9]).round(1).tolist())
meta = pd.read_csv('data/cmc_s5a_pseudobulk_meta.tsv', sep='\t')
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
for ct in ['L6 CT', 'L2/3 IT', 'L6 IT']:
    m = meta[(meta.group == ct) & (meta.n_cells >= 30)].copy()
    m['sex'] = m['sample'].map(sheet['Biological_Sex'])
    print(ct, 'n=%d' % len(m), m['sex'].value_counts().to_dict(),
          'n_cells median=%d' % m.n_cells.median())
