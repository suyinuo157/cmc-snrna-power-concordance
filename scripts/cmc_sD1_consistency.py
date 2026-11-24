import pandas as pd, numpy as np
d = pd.read_csv('results_sex/cmc_s8_deseq_all.csv')
d = d[d.group != 'L6 CT'].copy()
print('cols:', list(d.columns))
print('groups:', d.group.nunique())
G = {
 'Y': 'RPS4Y1 RPS4Y2 DDX3Y UTY KDM5D USP9Y TTTY14 TTTY10 NLGN4Y PCDH11Y PRKY ZFY TMSB4Y TBL1Y'.split(),
 'XIST': ['XIST'],
 'Xesc': 'JPX ZFX KDM6A PUDP DDX3X KDM5C NLGN4X'.split(),
 'PAR': 'PLCXD1 LINC00685 ASMTL-AS1 PPP2R3B CD99 DHRSX'.split(),
}
rows = []
for k, genes in G.items():
    for g in genes:
        x = d[d.gene == g]
        if len(x) == 0:
            continue
        s = x[x.padj < 0.05]
        rows.append(dict(cls=k, gene=g, n_tested=len(x), n_sig=len(s),
            n_up=int((x.log2FoldChange > 0).sum()),
            med_lfc=round(x.log2FoldChange.median(), 2),
            lo=round(x.log2FoldChange.min(), 2),
            hi=round(x.log2FoldChange.max(), 2)))
t = pd.DataFrame(rows)
pd.set_option('display.width', 220)
for k in G:
    print('\n=== %s ===' % k)
    print(t[t.cls == k].to_string(index=False))
t.to_csv('results_sex/cmc_sD1_consistency.csv', index=False)
print('\nXIST 各类型 baseMean:')
print(d[d.gene == 'XIST'][['group', 'baseMean', 'log2FoldChange', 'padj']].to_string(index=False))
