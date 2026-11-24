import glob, os
import pandas as pd

DE = r"D:\brainscope\analysis\de\all_celltypes"
rows = []
for f in sorted(glob.glob(os.path.join(DE, "*_DE.csv"))):
    ct = os.path.basename(f)[:-len("_DE.csv")]
    if ct.endswith("_limma") or ct.endswith("_meta"):
        continue
    d = pd.read_csv(f)
    if "padj" not in d.columns or len(d) == 0:
        continue
    rows.append({
        "celltype": ct,
        "n": len(d),
        "padj<0.05": int((d["padj"] < 0.05).sum()),
        "padj<0.10": int((d["padj"] < 0.10).sum()),
        "padj<0.25": int((d["padj"] < 0.25).sum()),
        "p<0.05": int((d["pval"] < 0.05).sum()),
        "p<0.05&|FC|>0.25": int(((d["pval"] < 0.05) & (d["log2FC"].abs() > 0.25)).sum()),
        "padj<0.25&|FC|>0.25": int(((d["padj"] < 0.25) & (d["log2FC"].abs() > 0.25)).sum()),
        "min_padj": round(float(d["padj"].min()), 4),
    })
s = pd.DataFrame(rows).sort_values("padj<0.25", ascending=False)
print(s.to_string(index=False))
print("\n各列合计:", {c: int(s[c].sum()) for c in s.columns if c not in ("celltype", "min_padj")})
