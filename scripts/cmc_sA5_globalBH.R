d <- read.csv("results_sex/cmc_s8_deseq_all.csv", stringsAsFactors = FALSE)
cat("rows:", nrow(d), "\n")
d <- d[!is.na(d$pvalue), ]
p <- d$pvalue
q <- p.adjust(p, method = "BH")
cat("tested:", length(p), "\n")
cat("min p:", signif(min(p), 3), "\n")
cat("min q:", signif(min(q), 3), "\n")
cat("n raw p<0.05:", sum(p < 0.05),
    "| expected if uniform:", round(0.05 * length(p)), "\n")
cat("n raw p<0.001:", sum(p < 0.001),
    "| expected if uniform:", round(0.001 * length(p)), "\n")
cat("n q<0.05:", sum(q < 0.05), "\n")
o <- head(order(p), 5)
print(data.frame(gene = d$gene[o], group = d$group[o],
                 p = signif(p[o], 3), q = signif(q[o], 3)), row.names = FALSE)
