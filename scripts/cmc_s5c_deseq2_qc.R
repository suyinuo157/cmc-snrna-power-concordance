suppressPackageStartupMessages(library(DESeq2))
setwd('/data2/penghuiminRL/scz')
DATA <- 'data'; OUT <- 'results_qc'
dir.create(OUT, showWarnings = FALSE)

pb <- read.delim(file.path(DATA, 'cmc_s5a_pseudobulk.tsv'), row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(DATA, 'cmc_s5a_pseudobulk_meta.tsv'), stringsAsFactors = FALSE)
sheet <- read.delim(file.path(DATA, 'cmc_s0_sheet.tsv'), stringsAsFactors = FALSE, check.names = FALSE)
rownames(sheet) <- sheet[['Individual_ID']]
qctab <- read.delim(file.path(DATA, 'cmc_s5c_qc_by_group.tsv'), stringsAsFactors = FALSE)
cat('pseudobulk', dim(pb), '| meta', dim(meta), '\n')
print(table(meta$group[meta$n_cells >= 30]))

groups <- sort(unique(meta$group))
MIN_CELLS <- 30; MIN_SAMP <- 12
allres <- list()
for (grp in groups) {
  m <- meta[meta$group == grp & meta$n_cells >= MIN_CELLS, ]
  if (nrow(m) < MIN_SAMP) { cat('SKIP', grp, 'n =', nrow(m), '\n'); next }
  cnt <- round(as.matrix(pb[, m$col, drop = FALSE]))
  colnames(cnt) <- m$sample
  d <- data.frame(row.names = m$sample,
                  SCZ  = as.integer(sheet[m$sample, 'condition'] == 'SCZ'),
                  male = as.integer(sheet[m$sample, 'Biological_Sex'] == 'male'),
                  Age  = as.numeric(sheet[m$sample, 'Age_death']),
                  PMI  = as.numeric(sheet[m$sample, 'PMI']),
                  pH   = as.numeric(sheet[m$sample, 'pH']),
                  qc   = qctab$median_n_genes[match(m$col, qctab$col)])
  d$qc[is.na(d$qc)] <- median(d$qc, na.rm = TRUE)
  d <- d[complete.cases(d), , drop = FALSE]
  cnt <- cnt[, rownames(d), drop = FALSE]
  d$Age <- as.numeric(scale(d$Age)); d$PMI <- as.numeric(scale(d$PMI)); d$pH <- as.numeric(scale(d$pH)); d$qc <- as.numeric(scale(d$qc))
  if (length(unique(d$SCZ)) < 2) { cat('SKIP', grp, 'no contrast\n'); next }
  keep <- rowSums(cnt >= 10) >= max(5, ceiling(0.5 * ncol(cnt)))
  cnt <- cnt[keep, , drop = FALSE]
  cat(grp, ': samples', ncol(cnt), 'genes', nrow(cnt), 'SCZ', sum(d$SCZ), '\n')
  r <- try({
    dds <- DESeqDataSetFromMatrix(cnt, colData = d, design = ~ qc + SCZ + male + Age + PMI + pH)
    dds <- DESeq(dds, quiet = TRUE)
    results(dds, name = 'SCZ', alpha = 0.05)
  }, silent = TRUE)
  if (inherits(r, 'try-error')) { cat('FAIL', grp, as.character(r), '\n'); next }
  df <- as.data.frame(r); df$gene <- rownames(df); df$group <- grp
  allres[[grp]] <- df
  safe <- gsub('[^A-Za-z0-9._-]+', '_', grp)
  write.csv(df, file.path(OUT, paste0('cmc_s5c_deseq_', safe, '.csv')), row.names = FALSE)
  cat('done', grp, 'padj<0.05 =', sum(df$padj < 0.05, na.rm = TRUE), '\n')
}

if (length(allres)) {
  res <- do.call(rbind, allres)
  write.csv(res, file.path(OUT, 'cmc_s5c_deseq_all.csv'), row.names = FALSE)
  tab <- do.call(rbind, lapply(allres, function(x) data.frame(
    group = x$group[1],
    n_tested = sum(!is.na(x$padj)),
    n_sig = sum(x$padj < 0.05, na.rm = TRUE),
    n_up = sum(x$padj < 0.05 & x$log2FoldChange > 0, na.rm = TRUE),
    n_down = sum(x$padj < 0.05 & x$log2FoldChange < 0, na.rm = TRUE))))
  write.csv(tab, file.path(OUT, 'cmc_s5c_summary.csv'), row.names = FALSE)
  cat('TOTAL tested', nrow(res), 'sig', sum(res$padj < 0.05, na.rm = TRUE), '\n')
  print(tab)
}
