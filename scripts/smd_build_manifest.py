#!/usr/bin/env python3
import os, csv, yaml, sys, random

CFG = "config.yaml"
STARTS_DIR = "manifests/starts"
OUT = "manifests/smd_manifest.csv"

os.makedirs("manifests", exist_ok=True)
random.seed()  # non-deterministic, as requested

cfg = yaml.safe_load(open(CFG))
G   = cfg["globals"]
dt_ps       = float(G["dt_ps"])
k_kj        = float(G["k_kj_mol_nm2"])
target_ext  = float(G["target_extension_nm"])
axis        = cfg["globals"]["smd"]["axis"]
array_cap   = int(G["slurm"]["array_cap"])
default_speeds = G.get("pull_speeds_nm_per_ns", [])
n_reps_default = int(G.get("n_reps_default", 1))

FIELDNAMES = [
    "system","variant","replicate","speed_nm_per_ns",
    "k_kj_mol_nm2","dt_ps","target_extension_nm","axis",
    "perf_ns_per_day","start_time_ps","final_tpr","final_xtc",
    "start_id","anchor_chain","array_cap"
]

import argparse

parser = argparse.ArgumentParser(description="Build SMD manifest from sampled starts")
parser.add_argument("--force", action="store_true", help="regenerate manifest for all systems (ignore existing manifest)")
args = parser.parse_args()

rows = []
existing = {}
if os.path.isfile(OUT) and not args.force:
    # load existing manifest to preserve start assignments
    with open(OUT) as f:
        r = csv.DictReader(f)
        for rec in r:
            key = (rec['system'], rec['variant'], rec['replicate'], rec['speed_nm_per_ns'])
            existing[key] = rec
for sysent in cfg["systems"]:
    system = sysent["name"]
    perf   = float(sysent.get("perf_ns_per_day", 50.0))
    for var in sysent.get("variants", []):
        variant = var["id"]
        anchor_chain = var["anchor"]["chain"]
        speeds = var.get("speeds", default_speeds)
        n_reps = int(var.get("n_reps", n_reps_default))
        # variant-specific override for target extension (nm)
        v_te = var.get("target_extension_nm", None)
        if v_te is None:
            v_te = var.get("target_extension", None)
        if v_te is None:
            target_ext_variant = float(target_ext)
        else:
            target_ext_variant = float(v_te)

        starts_csv = os.path.join(STARTS_DIR, f"{system}__{variant}.csv")
        if not os.path.isfile(starts_csv):
            print(f"[smd-manifest] ERROR: missing starts CSV {starts_csv}. Run `make posteq-sample` first.", file=sys.stderr)
            sys.exit(2)

        all_starts = []
        with open(starts_csv) as f:
            r = csv.DictReader(f)
            all_starts = list(r)
        if not all_starts:
            print(f"[smd-manifest] ERROR: no starts in {starts_csv}", file=sys.stderr)
            sys.exit(2)

        # We need one random start PER TASK (speed x replicate).
        # Strategy: shuffle ALL starts once; if enough, take unique without replacement.
        # If not enough unique starts to cover all tasks, cycle through (rare).
        total_tasks = len(speeds) * n_reps
        starts_shuffled = all_starts[:]
        random.shuffle(starts_shuffled)

        if len(starts_shuffled) < total_tasks:
            # Not enough unique starts — reuse (cycle) after exhausting uniques
            pool = starts_shuffled[:]
            # Extend cycle to required length
            while len(starts_shuffled) < total_tasks:
                starts_shuffled.extend(pool)
        chosen = starts_shuffled[:total_tasks]

        # Assign sequentially across speeds × reps
        idx = 0
        for speed in speeds:
            for rep in range(1, n_reps+1):
                key = (system, variant, str(rep), str(speed))
                if key in existing:
                    # preserve previous assignment if source files still exist
                    er = existing[key]
                    # verify the referenced start still exists in manifests/starts CSV (best-effort)
                    rows.append({
                        "system": er["system"],
                        "variant": er["variant"],
                        "replicate": int(er["replicate"]),
                        "speed_nm_per_ns": float(er["speed_nm_per_ns"]),
                        "k_kj_mol_nm2": float(er.get("k_kj_mol_nm2", k_kj)),
                        "dt_ps": float(er.get("dt_ps", dt_ps)),
                        "target_extension_nm": float(er.get("target_extension_nm", target_ext)),
                        "axis": er.get("axis", axis),
                        "perf_ns_per_day": float(er.get("perf_ns_per_day", perf)),
                        "start_time_ps": er.get("start_time_ps"),
                        "final_tpr": er.get("final_tpr"),
                        "final_xtc": er.get("final_xtc"),
                        "start_id": er.get("start_id"),
                        "anchor_chain": er.get("anchor_chain", anchor_chain),
                        "array_cap": er.get("array_cap", array_cap)
                    })
                else:
                    st = chosen[idx]
                    idx += 1
                    rows.append({
                        "system": system,
                        "variant": variant,
                        "replicate": rep,
                        "speed_nm_per_ns": speed,
                        "k_kj_mol_nm2": k_kj,
                        "dt_ps": dt_ps,
                        "target_extension_nm": target_ext_variant,
                        "axis": axis,
                        "perf_ns_per_day": perf,
                        "start_time_ps": st["start_time_ps"],
                        "final_tpr": st["final_tpr"],
                        "final_xtc": st["final_xtc"],
                        "start_id": st["start_id"],
                        "anchor_chain": anchor_chain,
                        "array_cap": array_cap
                    })

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES)
    w.writeheader()
    w.writerows(rows)

print(f"[smd-manifest] Wrote {OUT} ({len(rows)} rows)")
