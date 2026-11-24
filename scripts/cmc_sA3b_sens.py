"""第 4 条：技术轴模拟的敏感性分析（论文 A）。

基线（cmc_sA3_inflation.py）：12 种细胞类型、30% 基因携带等权载荷、单一技术轴。
本脚本逐个放宽这三处假设，各自扫描 lambda 并插值出发观测一致率所需的 lambda*。

  base    : 30% 基因，等权载荷，单轴（原设定，作为复现基准）
  frac10  : 10% 基因携带载荷
  frac50  : 50% 基因携带载荷
  hetero  : 30% 基因，载荷按细胞类型抖动（log-normal, sd = 0.5）
  twoax   : 两个独立轴，各 30% 基因，合并时除以 sqrt(2) 保持总方差

每个场景在 lambda = 0.75 / 1.0 / 1.25 / 1.5 各跑 50 次，报告一致基因百分比，
并线性插值出达到观测率 4.26%（630/14,804，线性引擎）所需的 lambda。
输出 results_sim/cmc_sA3b_sens.csv
"""
import numpy as np
import pandas as pd
from scipy import stats

rng = np.random.default_rng(42)

pb = pd.read_csv('data/cmc_s5a_pseudobulk.tsv', sep='\t', index_col=0)
CTS = ['L2/3 IT', 'L4 IT', 'L5 IT', 'L6 IT', 'L6 CT', 'L6b',
       'Pvalb', 'Sst', 'Vip', 'Astro', 'Oligo', 'OPC']
samp = {c: set(x.split('|')[0] for x in pb.columns if x.endswith('|' + c)) for c in CTS}
common = sorted(set.intersection(*samp.values()))

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

n = len(common)
n1 = 41
lab = np.zeros(n)
lab[:n1] = 1
TARGET = 630 / 14804 * 100.0
print('cell types %d  samples %d  genes %d  target %.2f%%'
      % (len(CTS), n, G, TARGET))


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty(m)
    r[o] = np.clip(q, 0, 1)
    return r


def sign_sig(Yd, lb):
    k = np.zeros(G)
    for c in CTS:
        k += (Yd[c][:, lb == 1].mean(1) - Yd[c][:, lb == 0].mean(1)) > 0
    k = np.minimum(k, len(CTS) - k)
    return int((bh(np.minimum(2 * stats.binom.cdf(k, len(CTS), 0.5), 1)) < 0.05).sum())


def make(scen, lam, a1, a2, het):
    noise = dict((c, rng.normal(scale=sd[:, None], size=(G, n))) for c in CTS)
    if scen == 'twoax':
        t1 = rng.normal(size=n)
        t2 = rng.normal(size=n)
        return {c: mu[c][:, None]
                   + lam / np.sqrt(2) * (a1[:, None] * t1[None, :]
                                         + a2[:, None] * t2[None, :])
                   + noise[c] for c in CTS}
    t = rng.normal(size=n)
    out = {}
    for c in CTS:
        load = a1 if het is None else a1 * het[c]
        out[c] = mu[c][:, None] + lam * load[:, None] * t[None, :] + noise[c]
    return out


FRACS = {'base': 0.30, 'frac10': 0.10, 'frac50': 0.50, 'hetero': 0.30, 'twoax': 0.30}
LAMS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
R = 200

rows = []
for scen, f in FRACS.items():
    a1 = (rng.random(G) < f).astype(float)
    a2 = (rng.random(G) < f).astype(float)
    het = ({c: float(np.exp(rng.normal(0, 0.5))) for c in CTS}
           if scen == 'hetero' else None)
    for lam in LAMS:
        cnt = [sign_sig(make(scen, lam, a1, a2, het), lab) for _ in range(R)]
        m = float(np.mean(cnt))
        sdv = float(np.std(cnt, ddof=1))
        rows.append({'scenario': scen, 'lambda': lam, 'mean_sig': m,
                     'mean_pct': m / G * 100,
                     'sd_sig': sdv, 'se_sig': sdv / np.sqrt(R), 'n_runs': R,
                     'p90_sig': float(np.percentile(cnt, 90))})
        print('%-7s lam %.2f  mean %.1f  %.2f%%  se %.1f'
              % (scen, lam, m, m / G * 100, sdv / np.sqrt(R)), flush=True)

d = pd.DataFrame(rows)
d.to_csv('results_sim/cmc_sA3b_sens.csv', index=False)
print('wrote results_sim/cmc_sA3b_sens.csv')

print('=== 达到 %.2f%% 所需的 lambda（分段线性插值） ===' % TARGET)
for scen in FRACS:
    s = d[d.scenario == scen]
    x = s['lambda'].values
    y = s.mean_pct.values
    if y.min() > TARGET:
        print('%-7s 最小 lambda %.2f 已超过目标' % (scen, x.min()))
    elif y.max() < TARGET:
        print('%-7s 最大 lambda %.2f 仍未达到目标' % (scen, x.max()))
    else:
        lam_star = float(np.interp(TARGET, y, x))
        print('%-7s lambda* = %.2f' % (scen, lam_star))
