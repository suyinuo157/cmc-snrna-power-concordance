import gzip, os, glob, csv
from collections import Counter

BASE = r"D:\brainscope"
MAT = os.path.join(BASE, "expression", "snrna_expr_matrices", "snrna_expr_matrices")
OUT = os.path.join(BASE, "analysis", "inspection")
os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(os.path.join(MAT, "**", "*-annotated_matrix.txt.gz"), recursive=True))
print("matrix files:", len(files))

# 1) 文件清单（对应文档 命令1 matrix_files.txt）
with open(os.path.join(OUT, "matrix_files.txt"), "w", encoding="utf-8") as f:
    for i, p in enumerate(files, 1):
        f.write("%d\t%s\t%d\n" % (i, os.path.relpath(p, BASE), os.path.getsize(p)))

# 2) 前 12 个矩阵的前 5 行（对应文档 命令2 matrix_headers.txt）
with open(os.path.join(OUT, "matrix_headers.txt"), "w", encoding="utf-8") as f:
    for p in files[:12]:
        f.write("=" * 70 + "\n" + os.path.relpath(p, BASE) + "\n")
        with gzip.open(p, "rt", errors="replace") as g:
            for _ in range(5):
                line = g.readline()
                if not line:
                    break
                f.write(line.rstrip("\n")[:600] + "\n")

# 3) 样本 metadata 前 10 行（对应文档 命令3）
meta = os.path.join(BASE, "metadata", "PEC2_sample_metadata.txt")
with open(meta, encoding="utf-8", errors="replace") as f:
    head10 = [next(f) for _ in range(10)]
with open(os.path.join(OUT, "PEC2_sample_metadata_head.txt"), "w", encoding="utf-8") as f:
    f.writelines(head10)

# 4) 逐样本细胞数 + 细胞类型构成（QC 指标：cells per sample / cells per cell type）
rows = []
ct_global = Counter()
for p in files:
    cohort = os.path.basename(os.path.dirname(p))
    sample = os.path.basename(p).replace("-annotated_matrix.txt.gz", "")
    with gzip.open(p, "rt", errors="replace") as g:
        header = g.readline().rstrip("\n").split("\t")
    cells = header[1:]
    cnt = Counter(cells)
    ct_global.update(cnt)
    rows.append({"cohort": cohort, "sample": sample, "n_cells": len(cells),
                 "n_celltypes": len(cnt),
                 "top_celltypes": ";".join("%s:%d" % kv for kv in cnt.most_common(3))})

rows.sort(key=lambda r: (r["cohort"], r["sample"]))
with open(os.path.join(OUT, "sample_qc_cells.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["cohort", "sample", "n_cells", "n_celltypes", "top_celltypes"])
    w.writeheader()
    w.writerows(rows)

with open(os.path.join(OUT, "celltype_counts_all_samples.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["celltype", "n_cells"])
    for k, v in ct_global.most_common():
        w.writerow([k, v])

# 5) 抽查 3 个文件的基因数（行数）
for p in files[:3]:
    n = 0
    with gzip.open(p, "rt", errors="replace") as g:
        for _ in g:
            n += 1
    print("genes", os.path.relpath(p, BASE), n - 1)

print("n_cells total:", sum(r["n_cells"] for r in rows))
print("n_celltype labels:", len(ct_global))
print("done ->", OUT)
