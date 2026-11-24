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
  cmc_sA3b_sens.py            sample-level technical-axis simulation and lambda* calibration
  cmc_sA4b_perm_robust.py     donor-label permutation null, plain and stratified
figures/
  make_figures.py             generates Fig. 1-4 (and the sex-difference panels)
  data/                       small derived tables read by make_figures.py
```

Server-side scripts that fit the DESeq2 models, run the spike-in simulation, compute the global-BH
sensitivity analysis, and export the external-list transfer test are deposited in the same `scripts/`
directory of the archived release.

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
