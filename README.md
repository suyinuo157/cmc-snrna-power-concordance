# Power boundaries and false cross-cell-type concordance in 100-donor snRNA-seq

Analysis code accompanying the manuscript *"Statistical power and technical covariation limit
cross-cell-type differential expression in 100-donor single-nucleus transcriptomics"*.

The study asks two questions about a case-control single-nucleus cohort: what effect sizes the design
can resolve, and how much cross-cell-type concordance a purely technical, donor-level factor can create.
The code below reproduces every number and figure in the manuscript.

## Data

All input data are publicly available from the brainSCOPE resource (PsychENCODE Phase II,
https://brainscope.gersteinlab.org): sample metadata, per-sample annotated expression matrices,
per-cell-type pseudobulk expression, cell-type fractions and case-control differential expression results.
The analysed cohort is the 100-donor CommonMind Consortium (CMC) component (Mount Sinai brain bank).

Raw sequencing reads and genotypes are **not** redistributed here. They are available under controlled
access through Synapse (CMC snRNA-seq `syn22032544`; CMC genotypes `syn22000412`; PsychENCODE repository
`syn51111084`) and the NIMH Data Archive, subject to a data access request.

## Layout

```
scripts/
  download_brainscope.py      download the public brainSCOPE files used here
  inspect_matrices.py         inventory of per-sample expression matrices
  official_deg.py             reformat the published case-control DE table
  deg_official_stats.py       summary statistics for the published DE table
  de_diag.py                  diagnostics of the linear-engine DE fits
  de_recount.py               recount of significant genes across cell types
  analysis_de.py              donor-level DE for one cell type (linear model + BH)
  analysis_de_all.py          donor-level DE across all cell types
  cmc_s5a_pseudobulk.py       donor x cell-type pseudobulk construction
  cmc_s5b_deseq2.R            per-cell-type DESeq2 model without the quality term
  cmc_s5b2_deseq2_rin.R       DESeq2 sensitivity analysis using RIN
  cmc_s5c_deseq2_qc.R         main model (~ qc + SCZ + male + Age + PMI + pH)
  cmc_s7_concordance.py       cross-cell-type sign-concordance statistic
  cmc_sA1_engine_check.py     comparison of the linear engine with DESeq2
  cmc_sA1b_contrast_probe.py  contrast-specification probe
  cmc_sA1c_sizefactor.py      size-factor diagnostic
  cmc_sA2_spikein.R           spike-in power simulation
  cmc_sA2_summary.py          power and empirical FDR summary of the spike-in runs
  cmc_sA3_inflation.py        single technical axis, first version
  cmc_sA3b_sens.py            technical-axis grid and lambda* calibration
  cmc_sA4_realperm.py         donor-label permutation on the observed data
  cmc_sA4b_perm_robust.py     1,000 permutations, plain and stratified
  cmc_sA5_globalBH.R          global BH correction across all 142,722 tests
  cmc_s8_deseq2_sex.R         sex-difference DESeq2 models
  cmc_s8_annotate.py          annotation of the sex-difference results
  cmc_s8_clean.py             filtering of the sex-difference tables
  cmc_s8_l6ct_check.py        size-factor check of the L6 CT cell type
  cmc_sD1_consistency.py      four-module concordance table
  cmc_sD2_xist_qc.py          XIST separation positive control
  cmc_sD3_figdata.py          figure input tables
  cmc_sD4_replicate_deca.py   replication of the external sex-bias list
  cmc_sD5_export_transfer.py  export of the transfer-test tables
figures/
  make_figures.py             generates Fig. 1-4 (and the sex-difference panels)
  data/                       small derived tables read by make_figures.py
```

Scripts were run in the original analysis environment. The absolute data paths at the top of each script
(`/data2/penghuiminRL/scz`) must be adjusted to a local copy of the brainSCOPE files before rerunning.

## Software

- R 4.2.1 with DESeq2 1.42.1 (models and spike-in simulation; `apeglm` for shrinkage where used).
- Python 3.10 with numpy, scipy, pandas, matplotlib (simulations, permutation null, bootstrap, figures).
- `requirements.txt` lists the Python packages.

## Reproducing the analyses

1. Download the public input files with `scripts/download_brainscope.py` (about 120 MB for the core set;
   the per-sample matrices are larger).
2. Fit the cell-type models and extract effect sizes (DESeq2 scripts).
3. Run `cmc_sA3b_sens.py` for the technical-axis calibration and `cmc_sA4b_perm_robust.py` for the
   donor-label permutation null.
4. Run `figures/make_figures.py` to regenerate the figures from the tables in `figures/data/`.

Absolute paths inside the scripts point to the original analysis environment and need to be adjusted to
the local data directory.

## Citation

If you use this code, please cite the manuscript (citation to be added upon publication) and the
brainSCOPE resource:

> Emani PS, Liu JJ, Clarke D, et al. Single-cell genomics and regulatory networks for 388 human brains.
> *Science*. 2024;384(6698):eadi5199. doi:10.1126/science.adi5199

## License

MIT (see `LICENSE`).
