import os
import pandas as pd

SRC = r"D:\brainscope\DEG\Schizophrenia_DEGcombined.csv"
OUT = r"C:\Users\Lenovo\AppData\Roaming\Trae CN\ModularData\ai-agent\work-mode-projects\6ac845f8860e188e4c3bd6c5"

d = pd.read_csv(SRC)
d = d.drop(columns=[c for c in d.columns if c.startswith("Unnamed")], errors="ignore")
for c in ["log2FoldChange", "pvalue", "padj", "baseMean"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")

summ = []
for ct, g in d.groupby("cell_type"):
    padj = g["padj"]
    lfc = g["log2FoldChange"]
    sig = g[(padj < 0.05)]
    summ.append({
        "celltype": ct,
        "n_tested": len(g),
        "n_padj05": len(sig),
        "n_padj05_fc025": int(((padj < 0.05) & (lfc.abs() > 0.25)).sum()),
        "n_up": int((sig["log2FoldChange"] > 0).sum()),
        "n_down": int((sig["log2FoldChange"] < 0).sum()),
        "min_padj": float(padj.min()),
        "top_gene": sig.loc[sig["padj"].idxmin(), "gene"] if len(sig) else "",
    })
summ = pd.DataFrame(summ).sort_values("n_padj05", ascending=False)

sig = d[d["padj"] < 0.05].copy()
sig["pass_fc"] = sig["log2FoldChange"].abs() > 0.25
sig = sig.sort_values(["cell_type", "padj"])

summ.to_csv(os.path.join(OUT, "official_DEG_summary.csv"), index=False)
sig[["cell_type", "gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "pass_fc"]] \
    .to_csv(os.path.join(OUT, "official_DEG_significant.csv"), index=False)

print(summ.to_string(index=False))
print("\n官方 DEG 合计:", len(sig), "| 其中 |log2FC|>0.25:", int(sig["pass_fc"].sum()),
      "| 唯一基因:", sig["gene"].nunique(), "| 细胞类型:", sig["cell_type"].nunique())
