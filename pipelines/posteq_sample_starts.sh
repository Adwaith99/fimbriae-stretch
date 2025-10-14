#!/usr/bin/env bash
set -euo pipefail

# Writes ALL candidate start times (after warmup, spaced by sample_every_ps)
# to manifests/starts/<system>__<variant>.csv:
#   system,variant,start_time_ps,final_tpr,final_xtc,start_id
# (start_id = s001, s002, ...)

python3 - <<'PY'
import sys, os, glob, yaml

cfg = yaml.safe_load(open("config.yaml"))
eq   = cfg["globals"]["equilibrium_md"]
warmup_ps = float(eq["warmup_ns"]) * 1000.0
length_ps = float(eq["length_ns"]) * 1000.0
step_ps   = float(eq["sample_every_ps"])

os.makedirs("manifests/starts", exist_ok=True)

def find_final_npt(system_dir):
    # Prefer dirs with *npt* and final/prod/production, fallback to npt*final*, then any .tpr/.xtc
    cands = []
    for path in glob.glob(os.path.join(system_dir, "**"), recursive=True):
        lp = path.lower()
        if os.path.isdir(path) and "npt" in lp and any(k in lp for k in ("final","prod","production")):
            cands.append(path)
    for root in cands:
        tprs = sorted(glob.glob(os.path.join(root, "*.tpr")))
        xtcs = sorted(glob.glob(os.path.join(root, "*.xtc")))
        if tprs and xtcs:
            return tprs[-1], xtcs[-1]
    tprs = sorted(glob.glob(os.path.join(system_dir, "**/*npt*final*.tpr"), recursive=True))
    xtcs = sorted(glob.glob(os.path.join(system_dir, "**/*npt*final*.xtc"), recursive=True))
    if tprs and xtcs:
        return tprs[-1], xtcs[-1]
    tprs = sorted(glob.glob(os.path.join(system_dir, "**/*.tpr"), recursive=True))
    xtcs = sorted(glob.glob(os.path.join(system_dir, "**/*.xtc"), recursive=True))
    if tprs and xtcs:
        return tprs[-1], xtcs[-1]
    return None, None

for s in cfg["systems"]:
    system = s["name"]
    sysdir = os.path.join("systems", system)
    for var in s.get("variants", []):
        variant = var["id"]
        tpr, xtc = find_final_npt(sysdir)
        if not tpr or not xtc:
            raise SystemExit(f"[posteq_sample_starts] ERROR: final NPT .tpr/.xtc not found for {system}/{variant}")

        # Build ALL candidate timestamps
        times = []
        t = warmup_ps
        while t <= length_ps:
            times.append(round(t, 3))
            t += step_ps
        if not times:
            raise SystemExit(f"[posteq_sample_starts] ERROR: no sampleable times for {system}/{variant}")

        out = os.path.join("manifests/starts", f"{system}__{variant}.csv")
        with open(out, "w") as f:
            f.write("system,variant,start_time_ps,final_tpr,final_xtc,start_id\n")
            width = len(str(len(times)))
            for i, st in enumerate(times, start=1):
                sid = f"s{str(i).zfill(max(3,width))}"
                f.write(f"{system},{variant},{st},{tpr},{xtc},{sid}\n")
        print(f"[posteq_sample_starts] Wrote {out} with {len(times)} candidate starts")
PY
