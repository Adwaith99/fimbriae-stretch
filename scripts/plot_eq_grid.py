#!/usr/bin/env python3
"""
plot_eq_grid.py
---------------
Rows:
  - (optional) EM: Potential  — spans all 3 columns
  - For each (T,FC) stage: Temperature (if NVT), Pressure (if NPT), Density (if NPT)
File name pattern (unchanged):
  nvt_<TEMP>K_<FC>_T.xvg
  npt_<TEMP>K_<FC>_P.xvg
  npt_<TEMP>K_<FC>_rho.xvg
Final unrestrained NPT is typically npt_303K_0_{P,rho}.xvg.
"""

from pathlib import Path
import re
import matplotlib.pyplot as plt

xvg = Path("xvg")
assert xvg.is_dir(), "Need ./xvg with XVG files from the extractor"

def read_xy(fname):
    xs, ys = [], []
    with open(fname) as fh:
        for ln in fh:
            if ln.startswith(("#", "@")):
                continue
            parts = ln.split()
            if len(parts) < 2:
                continue
            t, v = map(float, parts[:2])
            xs.append(t); ys.append(v)
    return xs, ys

# Gather NVT/NPT files by (T, FC)
pattern = re.compile(r"(nvt|npt)_(\d+)K_(\d+).*_(T|P|rho)\.xvg$")
series  = {}   # {(T,kbb): {"T":Path, "P":Path, "rho":Path}}
for f in xvg.glob("*.xvg"):
    m = pattern.match(f.name)
    if not m:
        continue
    step, T, kbb, kind = m.groups()
    key = (int(T), int(kbb))
    series.setdefault(key, {})
    if step == "nvt" and kind == "T":
        series[key]["T"]   = f
    elif step == "npt" and kind == "P":
        series[key]["P"]   = f
    elif step == "npt" and kind == "rho":
        series[key]["rho"] = f

# Optional EM potential panel (spans 3 columns)
em_file = (xvg / "em_potential.xvg")
has_em = em_file.exists()

# Sort rows by temperature then restraint
rows = sorted(series.keys())           # e.g. (100,1000) … (303,0)

# Figure layout
n_stage_rows = max(len(rows), 1) if rows else 0
nrows = n_stage_rows + (1 if has_em else 0)
ncols = 3

fig = plt.figure(figsize=(12, 2.5*nrows if nrows else 3))
gs = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.8, wspace=0.4)

row_offset = 0
# EM row at top, if present
if has_em:
    ax_em = fig.add_subplot(gs[0, :])  # span all columns
    tx, ty = read_xy(em_file)
    ax_em.plot(tx, ty, lw=1)
    ax_em.set_title("EM: Potential")
    ax_em.set_xlabel("Step")
    ax_em.set_ylabel("Potential (kJ/mol)")
    ax_em.grid(True, alpha=0.3)
    row_offset = 1

# Stage rows
axes_row = []
for i in range(n_stage_rows):
    axes_row.append([
        fig.add_subplot(gs[row_offset + i, 0]),
        fig.add_subplot(gs[row_offset + i, 1]),
        fig.add_subplot(gs[row_offset + i, 2]),
    ])

# Populate stage rows
for r, (T, kbb) in enumerate(rows):
    axT, axP, axR = axes_row[r]
    lab = f"{T} K / {kbb} kJ"

    # Temperature (NVT only)
    if "T" in series[(T,kbb)]:
        x, y = read_xy(series[(T,kbb)]["T"])
        axT.plot(x, y, lw=1); axT.set_ylabel("T (K)")
    else:
        axT.set_visible(False)

    # Pressure (NPT)
    if "P" in series[(T,kbb)]:
        x, y = read_xy(series[(T,kbb)]["P"])
        axP.plot(x, y, lw=1); axP.set_ylabel("P (bar)")
    else:
        axP.set_visible(False)

    # Density (NPT)
    if "rho" in series[(T,kbb)]:
        x, y = read_xy(series[(T,kbb)]["rho"])
        axR.plot(x, y, lw=1); axR.set_ylabel("ρ (kg m$^{-3}$)")
    else:
        axR.set_visible(False)

    # Titles / labels
    for c, ax in enumerate((axT, axP, axR)):
        if ax.get_visible():
            ax.set_title(lab if c == 0 else "")
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.3)

# Bottom x labels
if n_stage_rows:
    # last visible P/R axes get x-labels
    last_row_axes = axes_row[-1]
    for ax in last_row_axes[1:]:
        if ax.get_visible():
            ax.set_xlabel("Time (ps)")
else:
    # no stages, only EM — x label already set
    pass

# Save
fig.savefig("eq_grid.png", dpi=300)
fig.savefig("eq_grid.pdf")
print("Wrote eq_grid.png / pdf")
