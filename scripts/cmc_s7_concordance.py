import numpy as np, pandas as pd
from scipy import stats

d = pd.read_csv('results_qc/cmc_s5c_deseq_all.csv')
g = d.groupby('gene').agg(
    n_tested=('log2FoldChange', 'size'),
    n_up=('log2FoldChange', lambda x: int((x > 0).sum())),
    n_dn=('log2FoldChange', lambda x: int((x < 0).sum())),
    med_lfc=('log2FoldChange', 'median'),
    min_p=('pvalue', 'min')).reset_index()

# 两组间符号一致性：只对每个基因统计「多少个细胞类型同向」
g['n_conc'] = g[['n_up', 'n_dn']].max(axis=1)
def signp(r):
    return stats.binomtest(int(r.n_conc), int(r.n_tested), 0.5).pvalue
g['sign_p'] = g.apply(signp, axis=1)

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(n); r[o] = np.clip(q, 0, 1); return r

g['sign_fdr'] = bh(g.sign_p.values)
g['direction'] = np.where(g.n_up >= g.n_dn, 'up', 'down')
g.to_csv('results/cmc_s7_concordance.csv', index=False)
print('genes total', len(g))
print('n_tested>=8:', int((g.n_tested >= 8).sum()))
print('sign_fdr<0.05:', int((g.sign_fdr < 0.05).sum()))
print('sign_fdr<0.05 & |med_lfc|>=0.1:',
      int(((g.sign_fdr < 0.05) & (g.med_lfc.abs() >= 0.1)).sum()))

# ---- 通路超几何检验（基因集手工内置，背景=实际检验到的基因）----
SETS = {
 'Cholesterol_biosynthesis': 'HMGCR HMGCS1 MVK PMVK MVD IDI1 FDPS FDFT1 SQLE LSS CYP51A1 TM7SF2 MSMO1 NSDHL HSD17B7 EBP SC5D DHCR7 DHCR24 LDLR INSIG1 SREBF2 SOAT1 CYP46A1 APOE'.split(),
 'Heat_shock_chaperone': 'HSPA1A HSPA1B HSPA8 HSPH1 HSP90AA1 HSP90AB1 HSPB1 DNAJA1 DNAJB1 CCT2 CCT3 CCT4 CCT5 CCT6A CCT7 CCT8 TCP1 BAG3 STIP1'.split(),
 'Myelin': 'MBP PLP1 MAG MOG MOBP CLDN11 MAL CNP OPALIN MYRF ERMN'.split(),
 'Synaptic': 'GRIN1 GRIN2A GRIN2B GRIA1 GRIA2 SYN1 SNAP25 SYT1 DLG4 SHANK3 NRXN1 NLGN1 GAD1 GAD2 SLC17A7 SLC32A1'.split(),
 'SCZ_GWAS': 'CACNA1C TCF4 ZNF804A DTNBP1 NRG1 DISC1 COMT AKT1 PRODH CHRNA7 CACNB2 HCN1 RIMS1 MEF2C SLC39A8 GRM3'.split(),
}
bg = set(g.loc[g.n_tested >= 8, 'gene'])
sig = set(g.loc[(g.sign_fdr < 0.05) & (g.med_lfc.abs() >= 0.1), 'gene'])
rows = []
for name, genes in SETS.items():
    S = set(genes) & bg
    hit = S & sig
    # 双尾超几何
    p = stats.hypergeom.sf(len(hit) - 1, len(bg), len(S), len(sig))
    rows.append((name, len(S), len(hit), len(sig), round(p, 4),
                 ','.join(sorted(hit)[:12])))
r = pd.DataFrame(rows, columns=['set', 'n_in_bg', 'n_hit', 'n_sig', 'p', 'hits'])
r['fdr'] = bh(r.p.values)
r = r.sort_values('p')
r.to_csv('results/cmc_s7_pathway.csv', index=False)
print('\n=== pathway hypergeometric ===')
print(r.round(4).to_string(index=False))
print('\n=== top 20 concordant (|med_lfc|>=0.1) ===')
print(g[(g.sign_fdr < 0.05) & (g.med_lfc.abs() >= 0.1)]
      .reindex(g.med_lfc.abs().sort_values(ascending=False).index)
      .dropna(subset=['gene']).head(20)
      [['gene', 'n_tested', 'n_conc', 'direction', 'med_lfc', 'sign_fdr']]
      .to_string(index=False))
