import pandas as pd
d = pd.read_csv('results_sex/cmc_s8_deseq_all.csv')
s = d[d.padj < 0.05].copy()
Y = set('RPS4Y1 RPS4Y2 DDX3Y UTY KDM5D USP9Y TTTY14 TTTY10 NLGN4Y PCDH11Y PRKY ZFY EIF1AY TMSB4Y TBL1Y AMELY SRY TXLNGY GOLGA2LY'.split())
s['cls'] = s.gene.map(lambda g: 'Y' if g in Y else ('XIST' if g == 'XIST' else 'other'))
print('=== 显著基因按类别 x 细胞类型 ===')
print(pd.crosstab(s['group'], s['cls']).to_string())
o = s[s.cls == 'other'].sort_values('padj')
print('\n=== 非性染色体基因：出现次数 ===')
print(o.groupby('gene').agg(n_types=('group', 'size'), med_lfc=('log2FoldChange', 'median'),
      min_padj=('padj', 'min')).sort_values('n_types', ascending=False).head(25).to_string())
print('\n=== 非性染色体 top 20（按 padj）===')
print(o[['group', 'gene', 'log2FoldChange', 'padj', 'baseMean']].head(20).to_string(index=False))
print('\nunique non-sex-chr genes:', o.gene.nunique())
