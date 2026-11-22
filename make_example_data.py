"""
Generate synthetic PSD files so the tool can be tried without real data.
These are made-up distributions, not measurements.
"""

from pathlib import Path

import numpy as np
import pandas as pd

out = Path("examples")
out.mkdir(exist_ok=True)

width_nm = np.logspace(np.log10(0.3), np.log10(5), 80)


def lognormal(w, peak, spread, height):
    return height * np.exp(-((np.log(w) - np.log(peak)) ** 2) / (2 * spread ** 2))


samples = {
    "sample_A": lognormal(width_nm, 0.55, 0.25, 0.012),
    "sample_B": lognormal(width_nm, 0.60, 0.30, 0.009) + lognormal(width_nm, 2.0, 0.35, 0.002),
    "sample_C": lognormal(width_nm, 0.75, 0.40, 0.007),
    "sample_D": lognormal(width_nm, 0.50, 0.20, 0.010) + lognormal(width_nm, 1.3, 0.25, 0.003),
}

for name, dV in samples.items():
    pd.DataFrame({
        "Pore Width (Å)": width_nm * 10,
        "Cumulative Pore Volume (cm³/g)": np.cumsum(dV),
        "Incremental Pore Volume (cm³/g)": dV,
    }).to_csv(out / f"{name}.csv", index=False)
    print(f"Wrote {out / name}.csv")
