import pandas as pd
d = pd.read_csv('results_sex/cmc_s8_deseq_all.csv')
d = d[d.group != 'L6 CT']
s = d[d.padj < 0.05].copy()
Y = set('RPS4Y1 RPS4Y2 DDX3Y UTY KDM5D USP9Y TTTY14 TTTY10 TTTY12 NLGN4Y PCDH11Y PRKY ZFY EIF1AY TMSB4Y TBL1Y LINC00278 AC064829.1 AC007563.2'.split())
XE = set('KDM5C KDM6A DDX3X ZFX JPX PUDP STS PLCXD1 DHRSX ASMTL LINC00685 EIF2S3 EIF1AX XIAP'.split())
def cls(g):
    if g in Y:
        return 'Y'
    if g == 'XIST':
        return 'XIST'
    if g in XE:
        return 'X_escape'
    return 'other'
s['cls'] = s.gene.map(cls)
print(pd.crosstab(s['group'], s['cls']).to_string())
print('\n总数:'); print(s['cls'].value_counts().to_string())
print('\nX_escape 明细（应全部为负 = 女性偏高）:')
print(s[s.cls == 'X_escape'][['group', 'gene', 'log2FoldChange', 'padj']]
      .sort_values(['gene', 'group']).to_string(index=False))
print('\n仍归入 other 的基因（逐个核实）:')
print(s[s.cls == 'other'].groupby('gene').agg(
    n=('group', 'size'), med_lfc=('log2FoldChange', 'median')).sort_values('n', ascending=False).to_string())
s.to_csv('results_sex/cmc_s8_clean.csv', index=False)
