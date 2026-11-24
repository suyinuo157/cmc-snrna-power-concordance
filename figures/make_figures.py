"""Route A (methodological) and Route D (sex differences) figures.

Inputs : figures/data/*.csv  (exported from /data2/penghuiminRL/scz)
Outputs: figures/*.png, figures/*.pdf
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import FuncNorm

FIG = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(FIG, "data")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
    "pdf.fonttype": 42,
})

C_MALE = "#2c3e50"
C_FEMALE = "#c0392b"
C_GREY = "#95a5a6"


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- D: XIST
def fig_d1():
    d = pd.read_csv(os.path.join(DATA, "cmc_sD2_xist_auc.csv")).sort_values("aucF")
    y = np.arange(len(d))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.3),
                             gridspec_kw={"width_ratios": [1.0, 1.25]})

    ax = axes[0]
    ax.barh(y, d.aucF - 0.95, left=0.95, height=0.7,
            color="#4a6fa5", edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels(d.group, fontsize=6)
    ax.set_xlim(0.95, 1.025)
    ax.set_xticks([0.95, 0.96, 0.97, 0.98, 0.99, 1.0])
    ax.set_xlabel("AUC (female > male, XIST counts)\naxis truncated at 0.95")
    for i, v in enumerate(d.aucF):
        ax.text(v + 0.002, i, "%.3f" % v, va="center", ha="left", fontsize=5.5)
    ax.set_title("XIST separates sex in every cell type", loc="left")

    ax2 = axes[1]
    for i, r in enumerate(d.itertuples()):
        ax2.plot([max(r.medM, 0.3), max(r.medF, 0.3)], [i, i],
                 color="#dddddd", lw=1.0, zorder=1)
    ax2.scatter(np.maximum(d.medM, 0.3), y, s=13, color=C_MALE, label="male", zorder=2)
    ax2.scatter(np.maximum(d.medF, 0.3), y, s=13, color=C_FEMALE, label="female", zorder=2)
    ax2.set_xscale("log")
    ax2.set_xlim(0.2, 8000)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlabel("median XIST counts per pseudobulk sample\n(0 counts plotted at 0.3)")
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
               fontsize=6, markerscale=1.2, handletextpad=0.3, columnspacing=1.2)

    save(fig, "FigD1_XIST_sex_separation")


# ------------------------------------------------------- D: LFC heat map
def fig_d2():
    m = pd.read_csv(os.path.join(DATA, "cmc_sD3_lfc_matrix.csv"), index_col=0)
    disp = {"L2_3_IT": "L2/3 IT", "L4_IT": "L4 IT", "L5_IT": "L5 IT",
            "L5_6_NP": "L5/6 NP", "L6_IT": "L6 IT"}
    order = ["L2_3_IT", "L4_IT", "L5_IT", "L5_6_NP", "L6_IT", "L6b",
             "Pvalb", "Sst", "Vip", "Lamp5", "Sncg", "Pax6", "Chandelier",
             "Astro", "OPC", "Oligo", "Vascular", "Myeloid"]
    mods = [
        ("Y-linked", "#c0392b",
         ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "USP9Y", "TTTY14", "TTTY10",
          "NLGN4Y", "PCDH11Y", "PRKY", "ZFY"]),
        ("XIST", "#8e44ad", ["XIST"]),
        ("XCI escape", "#2c3e50",
         ["JPX", "ZFX", "KDM6A", "PUDP", "DDX3X", "KDM5C"]),
        ("PAR1", "#e67e22", ["PLCXD1", "LINC00685", "ASMTL-AS1", "PPP2R3B", "CD99"]),
    ]
    genes = [g for _, _, gs in mods for g in gs]
    arr = m.loc[genes, order].values.astype(float)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eeeeee")
    fwd = lambda x: np.sign(x) * np.abs(x) ** 0.5
    inv = lambda y: np.sign(y) * np.abs(y) ** 2.0
    norm = FuncNorm((fwd, inv), vmin=-9, vmax=9)

    fig = plt.figure(figsize=(7.2, 6.0))
    ax = fig.add_axes([0.30, 0.30, 0.56, 0.62])
    cax = fig.add_axes([0.30, 0.115, 0.20, 0.016])

    im = ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([disp.get(c, c) for c in order], rotation=90, fontsize=6)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontstyle="italic", fontsize=6)

    ax.set_xticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.5)
    ax.tick_params(which="minor", length=0)

    for x in (5.5, 12.5, 15.5):
        ax.axvline(x, color="#555555", lw=0.7)

    start = 0
    for label, colour, gs in mods:
        stop = start + len(gs)
        mid = (start + stop - 1) / 2
        ax.text(len(order) - 0.2, mid, label, va="center", ha="left",
                color=colour, fontweight="bold", fontsize=6.5, clip_on=False)
        if stop < len(genes):
            ax.axhline(stop - 0.5, color="#555555", lw=0.7)
        start = stop

    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    ticks = [-8, -2, -0.5, 0.5, 2, 8]
    cb.set_ticks(ticks)
    cb.set_ticklabels(["%g" % t for t in ticks])
    cb.ax.tick_params(labelsize=5.5, width=0.5, length=2)
    cb.set_label("log2 fold change (male / female); grey = not tested",
                 fontsize=6, labelpad=3)

    save(fig, "FigD2_sexDE_heatmap")


# --------------------------------------------------------------- A: power
def fig_a1():
    L = np.array([0.2, 0.3, 0.5, 1.0])
    pw = np.array([0.121, 0.343, 0.670, 0.952])
    fd = np.array([0.055, 0.078, 0.075, 0.065])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))
    ax = axes[0]
    ax.axvspan(0.05, 0.323, color="#f1c40f", alpha=0.22, lw=0)
    ax.axhline(0.8, color="#7f8c8d", lw=0.6, ls=":")
    ax.plot(L, pw, "o-", color=C_FEMALE, ms=3.6, lw=1.1)
    ax.scatter([0.2], [0.0], marker="x", s=18, color="#7f8c8d", zorder=3)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("injected log2 fold change")
    ax.set_ylabel("power (alpha = 0.05)")
    ax.text(0.06, 0.90, "observed effect size\nmedian 0.099, p90 0.323", fontsize=5.6)
    ax.text(0.40, 0.815, "80% power", fontsize=5.8, color="#7f8c8d")
    ax.text(0.46, 0.09, "LFC = 0: 0/10 detected", fontsize=5.6)
    ax.set_title("a  spike-in power", loc="left")

    ax = axes[1]
    ax.plot(L, fd * 100, "o-", color="#2c3e50", ms=3.6, lw=1.1)
    ax.axhline(5, color="#7f8c8d", lw=0.6, ls=":")
    ax.set_ylim(0, 10)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("injected log2 fold change")
    ax.set_ylabel("empirical FDR (%)")
    ax.text(0.55, 5.4, "nominal 5%", fontsize=5.8, color="#7f8c8d")
    ax.set_title("b  false discovery rate", loc="left")

    save(fig, "FigA1_power_curve")


# -------------------------------------------------- A: technical-axis scan
def fig_a2():
    # baseline scenario of scripts/cmc_sA3b_sens.py (200 replicates per point,
    # G = 3000 genes); results_sim/cmc_sA3b_sens.csv
    lam = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
    fp = np.array([20.5, 82.2, 159.2, 213.8, 292.9, 425.8])
    se = np.array([4.6, 12.2, 18.1, 20.7, 25.1, 27.8])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))
    ax = axes[0]
    ax.errorbar(lam, fp, yerr=se, fmt="o-", color="#2c3e50", ms=3.6, lw=1.1,
                elinewidth=0.7, capsize=2)
    ax.axvline(0.90, color=C_FEMALE, ls="--", lw=0.9)
    ax.text(0.56, 350, "lambda* = 0.90\n(observed 4.26%\nconcordance)", fontsize=5.6,
            color=C_FEMALE)
    ax.set_xlabel("strength of sample-level technical axis (lambda)")
    ax.set_ylabel("false-positive genes per run")
    ax.set_title("a  simulated inflation", loc="left")

    ax = axes[1]
    ax.errorbar(lam, fp / 30.0, yerr=se / 30.0, fmt="o-", color="#2c3e50", ms=3.6,
                lw=1.1, elinewidth=0.7, capsize=2)
    ax.axhline(4.26, color=C_FEMALE, ls=":", lw=0.9)
    ax.axvline(0.90, color=C_FEMALE, ls="--", lw=0.9)
    ax.scatter([0.90], [4.26], s=22, facecolor="white", edgecolor=C_FEMALE, zorder=4)
    ax.set_xlim(0.3, 2.1)
    ax.set_ylim(0, 17)
    ax.set_xlabel("lambda")
    ax.set_ylabel("empirical FDR (%)")
    ax.text(1.12, 4.9, "observed 4.26%", fontsize=5.6, color=C_FEMALE)
    ax.set_title("b  calibration", loc="left")

    save(fig, "FigA2_lambda_inflation")


# ------------------------------------------------ A: permutation null
def fig_a3():
    d = pd.read_csv(os.path.join(DATA, "cmc_sA4b_plain_null.csv"))
    v = d.n_bh.values

    fig, ax = plt.subplots(figsize=(3.7, 2.6))
    ax.hist(v, bins=np.arange(0, 3801, 190), color=C_GREY, edgecolor="none")
    ax.set_yscale("log")
    ax.set_ylim(0.7, 1200)
    ax.axvline(127, color=C_MALE, ls="--", lw=0.8)
    ax.axvline(630, color=C_FEMALE, lw=1.3)
    handles = [Line2D([0], [0], color=C_MALE, ls="--", lw=0.9, label="null median 127"),
               Line2D([0], [0], color=C_FEMALE, lw=1.4, label="observed 630")]
    ax.legend(handles=handles, loc="upper right", fontsize=5.8, handlelength=1.8,
              labelspacing=0.25)
    ax.text(1180, 1.5, "empirical p = 0.226\n(1,000 label permutations)", fontsize=5.8)
    ax.set_xlabel("concordant genes passing FDR < 0.05")
    ax.set_ylabel("permutations")
    ax.set_xlim(0, 3800)

    save(fig, "FigA3_permutation_null")


# --------------------------------------- A: external-list transfer test
def fig_a4():
    s = pd.read_csv(os.path.join(DATA, "cmc_sD4_transfer_summary.csv")).set_index("module")
    g = pd.read_csv(os.path.join(DATA, "cmc_sD4_transfer_genes.csv"))
    g = g[(g.tested > 0) & g.rate.notna()]

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5),
                             gridspec_kw={"width_ratios": [1.0, 1.15]})

    ax = axes[0]
    x = np.arange(2)
    vals = s.loc[["autosomal", "sexchr"]]
    # 基因聚类 bootstrap 95% 区间（20,000 次，seed 20261124）；Wilson 区间仍在 summary CSV 中
    boot = {"autosomal": (0.4812, 0.5654), "sexchr": (0.8194, 0.9810)}
    lo = np.array([boot[m][0] for m in vals.index])
    hi = np.array([boot[m][1] for m in vals.index])
    ax.bar(x, vals.rate, width=0.55, color=[C_GREY, C_FEMALE], edgecolor="none", zorder=2)
    ax.errorbar(x, vals.rate,
                yerr=[vals.rate - lo, hi - vals.rate],
                fmt="none", ecolor=C_MALE, elinewidth=0.9, capsize=3, zorder=3)
    ax.axhline(0.5, color="#7f8c8d", ls=":", lw=0.8, zorder=1)
    for i, (r, k, n, l_, h_) in enumerate(zip(vals.rate, vals.same, vals.comparisons, lo, hi)):
        ax.text(i, max(r, h_) + 0.03,
                "%.1f%% (%.1f–%.1f)\n%d/%d" % (r * 100, l_ * 100, h_ * 100, k, n),
                ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["autosomal\n(%d genes)" % (g.module == "autosomal").sum(),
                        "sex-chromosomal\n(%d genes)" % (g.module == "sexchr").sum()],
                       fontsize=6)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("fraction same direction")
    ax.text(1.58, 0.515, "chance", fontsize=5.6, color="#7f8c8d", ha="right")
    ax.set_title("a  transfer of an external list", loc="left")

    ax = axes[1]
    for mod, colour in [("autosomal", C_GREY), ("sexchr", C_FEMALE)]:
        d = g[g.module == mod]
        v = np.sort(d.rate.values)
        name = "autosomal" if mod == "autosomal" else "sex-chromosomal"
        ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post",
                color=colour, lw=1.2, label="%s (%d genes)" % (name, len(v)))
    ax.axvline(0.5, color="#7f8c8d", ls=":", lw=0.8)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("per-gene agreement with reported direction")
    ax.set_ylabel("cumulative fraction of genes")
    ax.legend(loc="upper left", fontsize=5.6, handlelength=1.6, labelspacing=0.25)
    ax.set_title("b  per-gene agreement", loc="left")

    save(fig, "FigA4_list_transfer")


if __name__ == "__main__":
    fig_d1()
    fig_d2()
    fig_a1()
    fig_a2()
    fig_a3()
    fig_a4()
