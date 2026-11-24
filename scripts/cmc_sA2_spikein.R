suppressMessages(library(DESeq2))
DATA <- 'data'; OUT <- 'results_sim'
dir.create(OUT, showWarnings = FALSE)
cnt <- round(as.matrix(read.delim(file.path(DATA, 'cmc_s5a_pseudobulk.tsv'), row.names = 1, check.names = FALSE)))
meta <- read.delim(file.path(DATA, 'cmc_s5a_pseudobulk_meta.tsv'), stringsAsFactors = FALSE)
sheet <- read.delim(file.path(DATA, 'cmc_s0_sheet.tsv'), stringsAsFactors = FALSE)
rownames(sheet) <- sheet$Individual_ID
qctab <- read.delim(file.path(DATA, 'cmc_s5c_qc_by_group.tsv'), stringsAsFactors = FALSE)
grp <- 'L2/3 IT'
m <- meta[meta$group == grp & meta$n_cells >= 30, ]
cnt <- cnt[, m$col, drop = FALSE]
colnames(cnt) <- m$sample
num <- function(x) suppressWarnings(as.numeric(gsub('+', '', as.character(x), fixed = TRUE)))
d <- data.frame(
  SCZ = as.numeric(sheet[m$sample, 'condition'] == 'SCZ'),
  male = as.numeric(tolower(sheet[m$sample, 'Biological_Sex']) == 'male'),
  Age = num(sheet[m$sample, 'Age_death']),
  PMI = num(sheet[m$sample, 'PMI']),
  pH  = num(sheet[m$sample, 'pH']),
  qc  = qctab$median_n_genes[match(m$col, qctab$col)],
  row.names = m$sample)
d <- d[complete.cases(d), , drop = FALSE]
cnt <- cnt[, rownames(d), drop = FALSE]
for (v in c('Age', 'PMI', 'pH', 'qc')) d[[v]] <- as.numeric(scale(d[[v]]))
sczCols <- rownames(d)[d$SCZ == 1]
keep <- rowSums(cnt >= 10) >= max(5, ceiling(0.5 * ncol(cnt)))
cnt <- cnt[keep, ]
set.seed(1)
cnt <- cnt[sample(rownames(cnt), min(3000, nrow(cnt))), ]
cat('samples', ncol(cnt), 'SCZ', length(sczCols), 'genes', nrow(cnt), '\n')
res <- list(); i <- 0
for (L in c(0, 0.2, 0.3, 0.5, 1.0)) {
  for (b in 1:10) {
    i <- i + 1; cc <- cnt; inj <- character(0)
    if (L > 0) {
      inj <- sample(rownames(cc), 300)
      cc[inj[1:150], sczCols] <- round(cc[inj[1:150], sczCols] * 2^L)
      cc[inj[151:300], sczCols] <- round(cc[inj[151:300], sczCols] / 2^L)
    }
    dds <- DESeqDataSetFromMatrix(cc, colData = d, design = ~ qc + SCZ + male + Age + PMI + pH)
    dds <- DESeq(dds, quiet = TRUE)
    r <- results(dds, name = 'SCZ', alpha = 0.05)
    sig <- rownames(r)[!is.na(r$padj) & r$padj < 0.05]
    tp <- length(intersect(sig, inj)); fp <- length(setdiff(sig, inj))
    cat(sprintf('L=%.1f rep=%d sig=%d tp=%d fp=%d\n', L, b, length(sig), tp, fp))
    res[[i]] <- data.frame(L = L, rep = b, n_sig = length(sig), tp = tp, fp = fp,
      power = if (L > 0) tp / length(inj) else NA,
      fdr = if (length(sig) > 0) fp / length(sig) else NA)
  }
}
write.csv(do.call(rbind, res), file.path(OUT, 'cmc_sA2_sim.csv'), row.names = FALSE)
cat('DONE\n')
