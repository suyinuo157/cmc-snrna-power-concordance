import pandas as pd

p = r"D:\brainscope\DEG\Schizophrenia_DEGcombined.csv"
d = pd.read_csv(p)
print("shape:", d.shape)
print("cell_type 数:", d["cell_type"].nunique(), "| 基因数:", d["gene"].nunique())
print("padj 非空:", int(d["padj"].notna().sum()))

rows = []
for ct, g in d.groupby("cell_type"):
    padj = g["padj"].astype(float)
    lfc = g["log2FoldChange"].astype(float)
    pv = g["pvalue"].astype(float)
    rows.append({
        "celltype": ct,
        "n": len(g),
        "padj<0.05": int((padj < 0.05).sum()),
        "padj<0.05&|lfc|>0.25": int(((padj < 0.05) & (lfc.abs() > 0.25)).sum()),
        "pval<0.05&|lfc|>0.25": int(((pv < 0.05) & (lfc.abs() > 0.25)).sum()),
        "min_pval": float(pv.min()),
        "min_padj": float(padj.min()) if padj.notna().any() else float("nan"),
    })
s = pd.DataFrame(rows).sort_values("padj<0.05&|lfc|>0.25", ascending=False)
print(s.to_string(index=False))
print("\n合计 padj<0.05&|lfc|>0.25:", int(s["padj<0.05&|lfc|>0.25"].sum()),
      "| pval<0.05&|lfc|>0.25:", int(s["pval<0.05&|lfc|>0.25"].sum()))
