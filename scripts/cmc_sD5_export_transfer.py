import csv, collections, math

genes = {}
for line in open('data/deca2026_autosomal.txt'):
    if line.strip():
        g, b = line.split()
        genes[g] = (1 if b == 'male' else -1, 'autosomal')
for line in """CHM female
DDX3X female
JPX female
KDM5C female
KDM6A female
LINC00894 female
RPS4Y1 male
TTTY14 male
UBA1 female
USP9X female
USP9Y male
UTY male
XIST female
ZFX female""".strip().splitlines():
    g, b = line.split()
    genes[g] = (1 if b == 'male' else -1, 'sexchr')

tested = collections.defaultdict(list)
with open('results_sex/cmc_s8_deseq_all.csv') as f:
    for r in csv.DictReader(f):
        if r['gene'] in genes and r['group'].replace('_', ' ') != 'L6 CT' \
           and r['padj'] not in ('', 'NA'):
            tested[r['gene']].append((float(r['log2FoldChange']), float(r['padj'])))

def wilson(k, n):
    z = 1.96; p = k / n; den = 1 + z*z/n
    c = (p + z*z/(2*n)) / den
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return c - h, c + h

rows = []
for g, (d, mod) in genes.items():
    h = tested.get(g, [])
    same = sum(1 for l, _ in h if l * d > 0)
    sig = sum(1 for l, p in h if l * d > 0 and p < 0.05)
    rows.append((g, mod, d, len(h), same, sig))

with open('results_sex/cmc_sD4_transfer_genes.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene', 'module', 'bias', 'tested', 'same', 'sig', 'rate'])
    for g, mod, d, n, k, s in sorted(rows, key=lambda x: (x[1], -(x[4]/x[3] if x[3] else 0))):
        w.writerow([g, mod, 'male' if d > 0 else 'female', n, k, s,
                    round(k/n, 4) if n else ''])

with open('results_sex/cmc_sD4_transfer_summary.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['module', 'comparisons', 'same', 'rate', 'ci_lo', 'ci_hi'])
    for mod in ('autosomal', 'sexchr'):
        n = sum(r[3] for r in rows if r[1] == mod)
        k = sum(r[4] for r in rows if r[1] == mod)
        lo, hi = wilson(k, n)
        w.writerow([mod, n, k, round(k/n, 4), round(lo, 4), round(hi, 4)])
        print(mod, n, k, round(k/n, 4), round(lo, 4), round(hi, 4))
