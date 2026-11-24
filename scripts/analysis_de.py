"""L5.IT pseudobulk SCZ vs Control 差异表达（本地替代集群 DESeq2 流程）。

官方 pseudobulk 已做归一化（非原始 counts），因此用带协变量的线性模型 + BH 校正，
等价于 limma 的思路；如后续拿到 raw counts，可直接换 pydeseq2。
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT = r"D:\brainscope\analysis\de"
os.makedirs(OUT, exist_ok=True)

pb = pd.read_csv(r"D:\brainscope\expression\L5IT_pb.bed.gz", sep="\t")
meta_cols = ["#chr", "start", "end", "gene", "length", "strand"]
samples = [c for c in pb.columns if c not in meta_cols]
print("genes:", pb.shape[0], "samples:", len(samples))

expr = pb[samples].astype(np.float32)
expr.index = pb["gene"].astype(str)
expr = expr.groupby(level=0).mean()  # 官方 pseudobulk 已是 log 归一化值，不再取对数

meta = pd.read_csv(r"D:\brainscope\metadata\PEC2_sample_metadata.txt", sep="\t")
meta = meta[meta.Individual_ID.isin(samples)].copy()
print("metadata matched:", meta.shape)
print(meta.Disorder.value_counts().to_string())

keep_samples = meta.Individual_ID.tolist()
expr = expr[keep_samples]
meta = meta.set_index("Individual_ID").loc[keep_samples]

meta["is_scz"] = (meta.Disorder.str.lower() == "schizophrenia").astype(int)
meta["is_ctl"] = (meta.Disorder.str.lower() == "control").astype(int)
meta = meta[(meta.is_scz == 1) | (meta.is_ctl == 1)].copy()
expr = expr[meta.index]
print("SCZ:", int(meta.is_scz.sum()), "Control:", int(meta.is_ctl.sum()))

# 过滤低表达基因
expressed = (expr > 0).mean(axis=1) >= 0.2
expr = expr[expressed]
print("genes after expression filter:", expr.shape[0])

# 设计矩阵：diagnosis + cohort + sex + age
X_parts = [np.ones((len(meta), 1))]
names = ["intercept"]
X_parts.append(meta[["is_scz"]].values.astype(float))
names.append("SCZ_vs_Control")
for cov in ["Cohort", "Biological_Sex"]:
    dummies = pd.get_dummies(meta[cov], drop_first=True, dtype=float)
    if dummies.shape[1]:
        X_parts.append(dummies.values)
        names.extend([f"{cov}_{c}" for c in dummies.columns])
meta["age_num"] = pd.to_numeric(meta["Age_death"], errors="coerce")
X_parts.append(meta[["age_num"]].fillna(meta.age_num.mean()).values.astype(float))
names.append("Age_death")
X = np.hstack(X_parts)
print("design:", X.shape, names)

Y = expr.values  # genes x samples
beta, *_ = np.linalg.lstsq(X, Y.T, rcond=None)  # p x genes
resid = Y.T - X @ beta
dof = X.shape[0] - X.shape[1]
sigma2 = (resid ** 2).sum(axis=0) / dof
XtX_inv = np.linalg.pinv(X.T @ X)
se = np.sqrt(sigma2[None, :] * np.diag(XtX_inv)[:, None])
tvals = beta / np.maximum(se, 1e-12)
pvals = 2 * stats.t.sf(np.abs(tvals), dof)

idx = names.index("SCZ_vs_Control")
res = pd.DataFrame({
    "gene": expr.index,
    "log2FC": beta[idx],
    "stat": tvals[idx],
    "pvalue": pvals[idx],
})
res["padj"] = multipletests(res.pvalue, method="fdr_bh")[1]
res = res.sort_values("pvalue")
res.to_csv(os.path.join(OUT, "L5.IT_SCZ_vs_Control_DE.csv"), index=False)

sig = res[(res.padj < 0.05) & (res.log2FC.abs() > 0.5)]
print("significant:", len(sig), "up:", int((sig.log2FC > 0).sum()), "down:", int((sig.log2FC < 0).sum()))
sig.to_csv(os.path.join(OUT, "L5.IT_significant_DE.csv"), index=False)

# 与官方 DEG 对照
off = pd.read_csv(r"D:\brainscope\DEG\Schizophrenia_DEGcombined.csv")
off_l5 = off[(off.cell_type == "L5.IT") & off.padj.notna()]
off_l5 = off_l5[(off_l5.padj < 0.05) & (off_l5.log2FoldChange.abs() > 0.5)]
print("official L5.IT significant:", len(off_l5))

mine = set(sig.gene)
theirs = set(off_l5.gene)
shared_genes = mine & theirs
merged = sig[sig.gene.isin(shared_genes)].merge(
    off_l5[["gene", "log2FoldChange", "padj"]], on="gene", how="inner")

# 全基因层面的方向一致性
common = sorted(set(res.gene) & set(off_l5.gene))
mine_all = res.set_index("gene").loc[common]
off_all = off_l5.set_index("gene").loc[common]
r_all = np.corrcoef(mine_all.log2FC, off_all.log2FoldChange)[0, 1]
sign_agree = (np.sign(mine_all.log2FC) == np.sign(off_all.log2FoldChange)).mean()
sig_mask = (off_all.padj < 0.05) & (off_all.log2FoldChange.abs() > 0.5)
sign_agree_sig = (np.sign(mine_all.log2FC[sig_mask]) ==
                  np.sign(off_all.log2FoldChange[sig_mask])).mean()

from scipy.stats import fisher_exact
universe = len(set(expr.index) | set(off.gene))
a = len(shared_genes)
b = len(mine - theirs)
c = len(theirs - mine)
d = universe - a - b - c
odds, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")

lines = [
    f"local significant: {len(sig)}",
    f"official L5.IT significant: {len(off_l5)}",
    f"overlap: {a}",
    f"jaccard: {a / max(1, len(mine | theirs)):.3f}",
    f"pearson r of log2FC (all {len(common)} shared genes): {r_all:.3f}",
    f"sign agreement (all shared genes): {sign_agree:.3f}",
    f"sign agreement on official significant genes: {sign_agree_sig:.3f}",
    f"fisher one-sided p: {p_fisher:.3e}",
]
with open(os.path.join(OUT, "DE_vs_official_summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))

# 火山图
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 5))
y = -np.log10(res.padj.clip(lower=1e-300))
up = (res.padj < 0.05) & (res.log2FC > 0.5)
down = (res.padj < 0.05) & (res.log2FC < -0.5)
ax.scatter(res.log2FC, y, s=3, c="#bdbdbd", label="ns")
ax.scatter(res.log2FC[up], y[up], s=6, c="#d62728", label=f"up ({up.sum()})")
ax.scatter(res.log2FC[down], y[down], s=6, c="#1f77b4", label=f"down ({down.sum()})")
ax.axvline(0.5, ls="--", lw=0.8, c="k")
ax.axvline(-0.5, ls="--", lw=0.8, c="k")
ax.set_xlabel("log2 fold change (SCZ vs Control)")
ax.set_ylabel("-log10(padj)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "volcano_L5.IT.png"), dpi=150)
print("done ->", OUT)
