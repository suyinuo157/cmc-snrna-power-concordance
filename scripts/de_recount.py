import glob, os
import pandas as pd

DE = r"D:\brainscope\analysis\de\all_celltypes"
rows = []
all_sig = []
for f in sorted(glob.glob(os.path.join(DE, "*_DE.csv"))):
    ct = os.path.basename(f)[:-len("_DE.csv")]
    if ct.endswith("_limma") or ct.endswith("_meta"):
        continue
    d = pd.read_csv(f)
    if "log2FC" not in d.columns:
        continue
    d = d.dropna(subset=["log2FC", "padj"])
    # 官方口径：qval(padj)<0.05 且 |log2FC|>0.25
    s = d[(d["padj"] < 0.05) & (d["log2FC"].abs() > 0.25)].copy()
    s["celltype"] = ct
    all_sig.append(s)
    rows.append({
        "celltype": ct,
        "n_tested": len(d),
        "n_padj05": int((d["padj"] < 0.05).sum()),
        "n_official_cutoff": len(s),
        "n_up": int((s["log2FC"] > 0).sum()),
        "n_down": int((s["log2FC"] < 0).sum()),
        "min_padj": d["padj"].min(),
    })

summ = pd.DataFrame(rows).sort_values("n_official_cutoff", ascending=False)
out = pd.concat(all_sig, ignore_index=True)[
    ["celltype", "gene", "log2FC", "pval", "padj", "mean_SCZ", "mean_Control"]
].sort_values(["celltype", "padj"])
OUT = r"C:\Users\Lenovo\AppData\Roaming\Trae CN\ModularData\ai-agent\work-mode-projects\6ac845f8860e188e4c3bd6c5"
summ.to_csv(os.path.join(OUT, "official_cutoff_summary.csv"), index=False)
out.to_csv(os.path.join(OUT, "official_cutoff_genes.csv"), index=False)
print(summ.to_string(index=False))
print("\n总计官方口径 DEG:", len(out), "| 唯一基因:", out["gene"].nunique(),
      "| 有 DEG 的细胞类型:", out["celltype"].nunique(), "/", len(summ))
