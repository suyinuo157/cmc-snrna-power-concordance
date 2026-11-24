import csv, collections, math

genes = {}
for line in open('data/deca2026_autosomal.txt'):
    if line.strip():
        g, b = line.split()
        genes[g] = 1 if b == 'male' else -1

tested = collections.defaultdict(list)
with open('results_sex/cmc_s8_deseq_all.csv') as f:
    for r in csv.DictReader(f):
        if r['gene'] in genes and r['group'].replace('_', ' ') != 'L6 CT' \
           and r['padj'] not in ('', 'NA'):
            tested[r['gene']].append((r['group'], float(r['log2FoldChange']), float(r['padj'])))

rows, n_same, n_tot = [], 0, 0
for g, d in genes.items():
    hits = tested.get(g, [])
    same = sum(1 for _, l, _ in hits if l * d > 0)
    sig = sum(1 for _, l, p in hits if l * d > 0 and p < 0.05)
    n_same += same; n_tot += len(hits)
    rows.append((g, len(hits), same, sig))

def wilson(k, n):
    z = 1.96; p = k / n; den = 1 + z*z/n
    c = (p + z*z/(2*n)) / den
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return c - h, c + h

lo, hi = wilson(n_same, n_tot)
print('genes tested at least once:', sum(1 for r in rows if r[1] > 0), '/', len(genes))
print('gene-celltype tests:', n_tot)
print('same direction:', n_same, 'rate:', round(n_same/n_tot, 3),
      '95%CI:', round(lo, 3), '-', round(hi, 3))
print('genes with >=1 significant same-direction:', sum(1 for r in rows if r[3] > 0))
tg = [r for r in rows if r[1] > 0]
print('mean per-gene same-direction rate:', round(sum(r[2]/r[1] for r in tg)/len(tg), 3))
print('gene tested same sig')
for r in sorted(rows, key=lambda x: (-x[3], -x[2]))[:15]:
    print(*r)
