"""
Pore-size distribution (PSD) analysis and plotting.

Reads PSD exports (e.g. from a BET / gas sorption instrument), calculates
summary statistics for each sample and plots them side by side for comparison.

Each input file needs three columns:
    Pore Width (Å)
    Cumulative Pore Volume (cm³/g)
    Incremental Pore Volume (cm³/g)

Usage:
    python psd_analysis.py examples/sample_A.csv examples/sample_B.csv
    python psd_analysis.py data/*.xlsx --labels Raw Treated --out results
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# Column names expected in the input files, mapped to short internal names
COLUMNS = {
    "Pore Width (Å)": "width_A",
    "Cumulative Pore Volume (cm³/g)": "cumV",
    "Incremental Pore Volume (cm³/g)": "dV",
}

MICRO_LIMIT_NM = 1.0   # micropores: w < 1 nm
MESO_LIMIT_NM = 5.0    # small mesopores: 1 ≤ w < 5 nm


def load_psd(path):
    """Read a PSD file (.xlsx or .csv) and return width (nm), cumV and dV."""
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} is missing column(s): {missing}")

    df = df.rename(columns=COLUMNS)[list(COLUMNS.values())]
    df["width_nm"] = df["width_A"] / 10.0
    return df.sort_values("width_nm").reset_index(drop=True)


def half_max_crossing(w, dV, peak_idx, half_max, step):
    """Find where the PSD drops below half its peak height, moving left (step=-1) or right (step=+1)."""
    idx = peak_idx
    while 0 < idx < len(dV) - 1 and dV[idx] > half_max:
        idx += step
    if dV[idx] > half_max:          # never dropped below half max: use the edge of the data
        return w[idx]
    prev = idx - step               # linear interpolation between the last two points
    w1, w2, y1, y2 = w[prev], w[idx], dV[prev], dV[idx]
    return w1 if y2 == y1 else w1 + (half_max - y1) * (w2 - w1) / (y2 - y1)


def psd_stats(df, label):
    """Calculate volume-weighted pore statistics for one sample."""
    w = df["width_nm"].to_numpy()
    dV = df["dV"].to_numpy()
    V_total = df["cumV"].iloc[-1]

    # Treat dV as a distribution: normalise so the weights sum to 1
    weights = dV / dV.sum() if dV.sum() > 0 else np.zeros_like(dV)

    micro = w < MICRO_LIMIT_NM
    meso = (w >= MICRO_LIMIT_NM) & (w < MESO_LIMIT_NM)
    frac_micro, frac_meso = weights[micro].sum(), weights[meso].sum()

    w_avg = (weights * w).sum()
    w_std = np.sqrt((weights * (w - w_avg) ** 2).sum())

    def region_avg(mask, frac):
        return (weights[mask] * w[mask]).sum() / frac if frac > 0 else np.nan

    peak_idx = int(np.argmax(dV))
    half_max = 0.5 * dV[peak_idx]
    fwhm = (half_max_crossing(w, dV, peak_idx, half_max, +1)
            - half_max_crossing(w, dV, peak_idx, half_max, -1))

    return {
        "Sample": label,
        "V_total (cm³/g)": V_total,
        "V_<1nm (cm³/g)": frac_micro * V_total,
        "Frac_<1nm": frac_micro,
        "V_1-5nm (cm³/g)": frac_meso * V_total,
        "w_avg (nm)": w_avg,
        "w_std (nm)": w_std,
        "w_micro_avg (nm)": region_avg(micro, frac_micro),
        "w_meso_avg (nm)": region_avg(meso, frac_meso),
        "w_peak (nm)": w[peak_idx],
        "FWHM (nm)": fwhm,     # smaller = sharper peak
    }


def plot_psds(datasets, labels, out_path):
    """Plot each sample's incremental and cumulative PSD on a shared-scale grid."""
    n = len(datasets)
    cols = min(2, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)
    axes = axes.flatten()

    # Shared axis limits so the samples can be compared directly
    dV_max = max(d["dV"].max() for d in datasets) * 1.35
    cumV_max = max(d["cumV"].max() for d in datasets) * 1.35

    for i, (df, label) in enumerate(zip(datasets, labels)):
        ax1 = axes[i]
        ax2 = ax1.twinx()

        l1, = ax1.plot(df["width_nm"], df["dV"], color="black", linewidth=1.2,
                       label="Incremental")
        ax1.set_xscale("log")
        ax1.set_xticks([0.3, 0.5, 1, 2, 5])
        ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax1.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax1.set_ylim(0, dV_max)
        ax1.set_xlabel("Pore width (nm)")
        ax1.set_ylabel("Incremental pore volume (cm³ g⁻¹)")
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
        ax1.grid(True, which="both", ls="--", lw=0.4)

        l2, = ax2.plot(df["width_nm"], df["cumV"], color="blue", linestyle="--",
                       linewidth=1.0, label="Cumulative")
        ax2.set_ylim(0, cumV_max)
        ax2.set_ylabel("Cumulative pore volume (cm³ g⁻¹)", color="blue")
        ax1.legend(handles=[l1, l2], loc="upper right", fontsize=8, frameon=False)

        ax1.set_title(f"({chr(97 + i)}) {label}", loc="left", fontsize=11)

    for ax in axes[n:]:          # hide any unused panels
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Compare pore-size distributions across samples.")
    parser.add_argument("files", nargs="+", help="PSD files (.xlsx or .csv)")
    parser.add_argument("--labels", nargs="+", help="Sample labels (default: file names)")
    parser.add_argument("--out", default="results", help="Output folder (default: results)")
    args = parser.parse_args()

    labels = args.labels or [Path(f).stem for f in args.files]
    if len(labels) != len(args.files):
        parser.error("Give one label per file.")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = [load_psd(f) for f in args.files]
    summary = pd.DataFrame([psd_stats(df, lab) for df, lab in zip(datasets, labels)])

    fig_path = plot_psds(datasets, labels, out_dir / "psd_comparison.png")
    csv_path = out_dir / "psd_summary.csv"
    summary.to_csv(csv_path, index=False, encoding="utf-8-sig")

    print("\nPore structure summary:\n")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:0.4f}"))
    print(f"\nSaved: {fig_path}\nSaved: {csv_path}")


if __name__ == "__main__":
    main()
