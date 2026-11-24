"""对全部细胞类型做 SCZ vs Control 供体级差异表达（协变量线性模型 + BH-FDR）。

输入：
  D:\\brainscope\\expression\\pseudobulk\\<celltype>.expr.bed.gz
  D:\\brainscope\\metadata\\cell_type_fraction_count_with_meta.csv
输出：
  D:\\brainscope\\analysis\\de\\all_celltypes\\<celltype>_DE.csv
  D:\\brainscope\\analysis\\de\\all_celltypes_summary.csv
"""
import os
import glob
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

PB_DIR = r"D:\brainscope\expression\pseudobulk"
META = r"D:\brainscope\metadata\cell_type_fraction_count_with_meta.csv"
OUT = r"D:\brainscope\analysis\de\all_celltypes"
os.makedirs(OUT, exist_ok=True)

# 样本级 metadata：individualID 与 pseudobulk 列名一致
meta = pd.read_csv(META).drop_duplicates("individualID").set_index("individualID")
meta["Age"] = pd.to_numeric(meta["Age"], errors="coerce")
meta["Sex_male"] = (meta["Sex"].astype(str).str.lower().str[0] == "m").astype(float)

rows = []
for path in sorted(glob.glob(os.path.join(PB_DIR, "*.expr.bed.gz"))):
    ct = os.path.basename(path).replace(".expr.bed.gz", "")
    pb = pd.read_csv(path, sep="\t")
    samples = [c for c in pb.columns[6:] if c in meta.index]
    if len(samples) < 10:
        print(f"[skip] {ct}: only {len(samples)} matched samples")
        continue

    expr = pb[samples].astype(np.float32)
    expr.index = pb["gene"].astype(str)
    expr = expr.groupby(level=0).mean()

    m = meta.loc[samples].copy()
    diag = m["primaryDiagnosis"].astype(str).str.lower()
    keep = diag.isin(["schizophrenia", "control"])          # 只比 SCZ 与 Control
    m = m[keep]
    # 只保留同时含 SCZ 和 Control 的队列：只有单一诊断的队列作哑变量会与分组变量共线，
    # 导致分组系数标准误爆炸、结果全不显著
    ctd = pd.crosstab(m["Cohort"], m["primaryDiagnosis"].astype(str).str.lower())
    n_scz = ctd["schizophrenia"] if "schizophrenia" in ctd else pd.Series(0, index=ctd.index)
    n_ctl = ctd["control"] if "control" in ctd else pd.Series(0, index=ctd.index)
    ok_cohorts = ctd.index[(n_scz >= 2) & (n_ctl >= 2)]
    m = m[m["Cohort"].isin(ok_cohorts)]
    expr = expr[m.index]
    if len(m) < 20 or m["primaryDiagnosis"].nunique() < 2:
        print(f"[skip] {ct}: too few balanced samples ({len(m)})")
        continue

    age = m["Age"].fillna(m["Age"].mean())
    cohort = pd.get_dummies(m["Cohort"], drop_first=True, dtype=float)
    X = np.column_stack([
        np.ones(len(m)),
        (m["primaryDiagnosis"].astype(str).str.lower() == "schizophrenia").astype(float),
        m["Sex_male"].values,
        age.values,
        cohort.values,
    ])
    Y = expr.values.T.astype(np.float64)                     # 样本 × 基因

    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    dof = X.shape[0] - X.shape[1]
    sigma2 = (resid ** 2).sum(axis=0) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(sigma2 * XtX_inv[1, 1])
    tstat = beta[1] / se
    from scipy import stats
    pval = 2 * stats.t.sf(np.abs(tstat), dof)

    padj = multipletests(pval, method="fdr_bh")[1]
    is_scz = X[:, 1] == 1
    res = pd.DataFrame({
        "gene": expr.index,
        "log2FC": beta[1],
        "t": tstat,
        "pval": pval,
        "padj": padj,
        "mean_SCZ": Y[is_scz].mean(axis=0),
        "mean_Control": Y[~is_scz].mean(axis=0),
    }).sort_values("padj")
    res.to_csv(os.path.join(OUT, f"{ct}_DE.csv"), index=False)

    sig = res[(res["padj"] < 0.05) & (res["log2FC"].abs() > 0.5)]
    rows.append({
        "celltype": ct,
        "n_samples": len(m),
        "n_SCZ": int(is_scz.sum()),
        "n_Control": int((~is_scz).sum()),
        "n_genes_tested": len(res),
        "n_sig": len(sig),
        "n_up": int((sig["log2FC"] > 0).sum()),
        "n_down": int((sig["log2FC"] < 0).sum()),
        "top_gene": res.iloc[0]["gene"] if len(res) else "",
        "top_padj": float(res.iloc[0]["padj"]) if len(res) else np.nan,
    })
    print(f"{ct:22s} samples={len(m):3d} (SCZ={int(is_scz.sum()):2d}) sig={len(sig):4d}", flush=True)

summary = pd.DataFrame(rows).sort_values("n_sig", ascending=False)
summary.to_csv(r"D:\brainscope\analysis\de\all_celltypes_summary.csv", index=False)
print("\n" + summary.to_string(index=False))
