import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(42)
pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
CTS = ['L2/3 IT','L4 IT','L5 IT','L6 IT','L6 CT','L6b','Pvalb','Sst','Vip','Astro','Oligo','OPC']
samp = {c: set(x.split('|')[0] for x in pb.columns if x.endswith('|' + c)) for c in CTS}
common = sorted(set.intersection(*samp.values()))
print('cell types', len(CTS), 'common samples', len(common))
mats = {}
for c in CTS:
    col = [s + '|' + c for s in common]
    sub = pb[col]
    mats[c] = np.log2(sub.div(sub.sum(axis=0), axis=1) * 1e6 + 1).values
keep = mats['L2/3 IT'].mean(axis=1) > 1
G = min(3000, int(keep.sum()))
idx = rng.choice(np.where(keep)[0], G, replace=False)
Y = {c: mats[c][idx] for c in CTS}
mu = {c: Y[c].mean(axis=1) for c in CTS}
sd = np.median(np.vstack([Y[c].std(axis=1) for c in CTS]), axis=0)
print('genes', G, 'median sd', round(float(np.median(sd)), 3))
n = len(common); n1 = 41
lab = np.zeros(n); lab[:n1] = 1
a = (rng.random(G) < 0.3).astype(float)
def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(m); r[o] = np.clip(q, 0, 1); return r
def sign_sig(Yd, lb):
    k = np.zeros(G)
    for c in CTS:
        k += (Yd[c][:, lb == 1].mean(1) - Yd[c][:, lb == 0].mean(1)) > 0
    k = np.minimum(k, len(CTS) - k)
    return int((bh(np.minimum(2 * stats.binom.cdf(k, len(CTS), 0.5), 1)) < 0.05).sum())
res = []
for lam in [0, 0.25, 0.5, 0.75]:
    cnt = []
    for r in range(30):
        t = rng.normal(size=n)
        Yd = {c: mu[c][:, None] + lam * a[:, None] * t[None, :]
                 + rng.normal(scale=sd[:, None], size=(G, n)) for c in CTS}
        cnt.append(sign_sig(Yd, lab))
    res.append((lam, float(np.mean(cnt)), float(np.median(cnt)), int(np.max(cnt))))
r1 = pd.DataFrame(res, columns=['lambda','mean_sig','median_sig','max_sig'])
print('=== 无真实效应，仅技术轴 ===')
print(r1.round(2).to_string(index=False))
lam = 0.5
t = rng.normal(size=n)
Yd = {c: mu[c][:, None] + lam * a[:, None] * t[None, :]
         + rng.normal(scale=sd[:, None], size=(G, n)) for c in CTS}
obs = sign_sig(Yd, lab)
null_perm = np.array([sign_sig(Yd, rng.permutation(lab)) for _ in range(200)])
null_wrong = np.array([sign_sig({c: mu[c][:, None] + rng.normal(scale=sd[:, None], size=(G, n))
                                 for c in CTS}, rng.permutation(lab)) for _ in range(200)])
print('observed #sig', obs)
print('correct null (sample permutation, keeps axis): median %.0f  95th %.0f'
      % (np.median(null_perm), np.percentile(null_perm, 95)))
print('wrong null   (no technical axis):             median %.0f  95th %.0f'
      % (np.median(null_wrong), np.percentile(null_wrong, 95)))
print('empirical p  correct null %.3f   wrong null %.4f'
      % ((null_perm >= obs).mean(), (null_wrong >= obs).mean()))
r1.to_csv('results_sim/cmc_sA3_inflation.csv', index=False)
res2 = []
for lam in [0.75, 1.0, 1.25, 1.5]:
    cnt = []
    for r in range(50):
        t = rng.normal(size=n)
        Yd = {c: mu[c][:, None] + lam * a[:, None] * t[None, :]
                 + rng.normal(scale=sd[:, None], size=(G, n)) for c in CTS}
        cnt.append(sign_sig(Yd, lab))
    res2.append((lam, float(np.mean(cnt)), float(np.percentile(cnt, 90))))
r2 = pd.DataFrame(res2, columns=['lambda', 'mean_sig', 'p90_sig'])
r2['frac_pct'] = (r2.mean_sig / G * 100).round(2)
print('=== 校准：匹配真实 4.77% ===')
print(r2.round(2).to_string(index=False))
print('real data: 393 / 8245 = %.2f%%' % (393 / 8245 * 100))
r2.to_csv('results_sim/cmc_sA3_calib.csv', index=False)
