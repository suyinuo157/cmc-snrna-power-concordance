"""置换稳健性（论文 A 第 3 条）：方案与判据的敏感性。

  --mode plain : 全部 donor 一起打乱标签
  --mode strat : 在 性别 x 年龄三分位 分层内打乱
每次抽样同时计算 bh（已发表判据）与 c80（>=80% 同向）。
"""
import argparse
import time

import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument('--mode', choices=['plain', 'strat'], default='plain')
ap.add_argument('--nperm', type=int, default=20)
ap.add_argument('--seed', type=int, default=20261124)
a = ap.parse_args()

rng = np.random.default_rng(a.seed)

pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
meta = pd.read_csv('data/cmc_s5a_pseudobulk_meta.tsv', sep='\t')
sheet = pd.read_csv('data/cmc_s0_sheet.tsv', sep='\t').set_index('Individual_ID')
qc = pd.read_csv('data/cmc_s5c_qc_by_group.tsv', sep='\t').set_index('col')


def num(x):
    return pd.to_numeric(pd.Series(x).astype(str).str.replace('+', '', regex=False),
                         errors='coerce')


D = []
for g in meta.group.unique():
    m = meta[(meta.group == g) & (meta.n_cells >= 30)]
    if len(m) < 30:
        continue
    sub = pb[list(m['col'])]
    cpm = np.log2(sub.div(sub.sum(axis=0), axis=1) * 1e6 + 1)
    d = pd.DataFrame({'sample': list(m['sample'])}, index=list(m['col']))
    s = d['sample'].values
    d['SCZ'] = (sheet.loc[s, 'condition'].values == 'SCZ').astype(float)
    d['male'] = (sheet.loc[s, 'Biological_Sex'].str.lower().values == 'male').astype(float)
    for v, c in [('Age', 'Age_death'), ('PMI', 'PMI'), ('pH', 'pH')]:
        d[v] = num(sheet.loc[s, c].values).values
    d['qc'] = qc.reindex(list(m['col']))['median_n_genes'].values
    d = d.dropna()
    for v in ['Age', 'PMI', 'pH', 'qc']:
        d[v] = (d[v] - d[v].mean()) / d[v].std()
    D.append((g, cpm[d.index], d))
print('cell types used', len(D))

genes = None
for _, cpm, _ in D:
    sg = set(cpm.index)
    genes = sg if genes is None else genes & sg
genes = sorted(genes)
E = np.mean([(cpm.loc[genes].mean(axis=1) > 1).values for _, cpm, _ in D], axis=0)
genes = [g for g, x in zip(genes, E) if x >= 0.5]
print('common expressed genes', len(genes))

CELLS = []
for g, cpm, d in D:
    CELLS.append({
        'group': g,
        'sample': d['sample'].values,
        'Y': cpm.loc[genes].values.T,
        'Z': np.column_stack([d.male.values, d.Age.values, d.PMI.values,
                              d.pH.values, d.qc.values]),
    })


def betas(labv, Y, Z):
    X = np.column_stack([np.ones(len(labv)), labv, Z])
    return (np.linalg.pinv(X) @ Y)[1]


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(m)
    r[o] = np.clip(q, 0, 1)
    return r


def crit_bh(B):
    up = (B > 0).sum(0)
    dn = (B < 0).sum(0)
    neff = up + dn
    k = np.minimum(up, dn)
    p = np.ones(len(k))
    for nv in np.unique(neff):
        if nv < 8:
            continue
        sel = neff == nv
        p[sel] = np.minimum(2 * stats.binom.cdf(k[sel], nv, 0.5), 1)
    return int((bh(p) < 0.05).sum())


def crit_c80(B):
    up = (B > 0).sum(0)
    dn = (B < 0).sum(0)
    neff = up + dn
    frac = np.maximum(up, dn) / np.maximum(neff, 1)
    return int(((neff >= 8) & (frac >= 0.8)).sum())


lab_all = (sheet['condition'] == 'SCZ').astype(float)


def collect(lab):
    return np.vstack([betas(lab.reindex(c['sample']).values, c['Y'], c['Z'])
                      for c in CELLS])


B0 = collect(lab_all)
obs_bh, obs_c80 = crit_bh(B0), crit_c80(B0)
print('observed: bh %d   c80 %d' % (obs_bh, obs_c80))

sex = sheet['Biological_Sex'].astype(str).str.lower()
age = num(sheet['Age_death'])
tert = pd.Series(pd.qcut(age, 3, labels=False, duplicates='drop'), index=age.index)
strata = sex + '_' + tert.fillna(-1).astype(int).astype(str)
print('strata sizes:', strata.value_counts().to_dict())


def permute(lab):
    if a.mode == 'plain':
        return pd.Series(rng.permutation(lab.values), index=lab.index)
    out = lab.copy()
    st = strata.reindex(lab.index).fillna('NA')
    for s in st.unique():
        idx = st.index[st == s]
        out.loc[idx] = rng.permutation(lab.loc[idx].values)
    return out


rows = []
t0 = time.time()
for r in range(a.nperm):
    pl = permute(lab_all)
    B = collect(pl)
    rows.append({'perm': r + 1, 'n_bh': crit_bh(B), 'n_c80': crit_c80(B)})
    if (r + 1) % 50 == 0:
        el = time.time() - t0
        print('perm %d  %.0fs elapsed  est. total %.1f min'
              % (r + 1, el, el / (r + 1) * a.nperm / 60.0), flush=True)

df = pd.DataFrame(rows)
out = 'results_sim/cmc_sA4b_%s_null.csv' % a.mode
df.to_csv(out, index=False)
print('wrote', out, len(df), 'permutations')
for c, ob in [('n_bh', obs_bh), ('n_c80', obs_c80)]:
    v = df[c].values
    P = (1.0 + (v >= ob).sum()) / (1.0 + len(v))
    print('%s: obs %d  null median %.0f  95th %.0f  max %d  emp p = %.4f'
          % (c, ob, np.median(v), np.percentile(v, 95), v.max(), P))
