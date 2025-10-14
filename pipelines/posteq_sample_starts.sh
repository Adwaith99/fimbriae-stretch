#!/usr/bin/env bash
set -euo pipefail

# Read config with yq or python; we’ll use python (no external deps) to keep portable.
cfg="config.yaml"

py() {
python3 - "$@" <<'PY'
import sys, yaml, json, os, glob, random
random.seed()  # non-deterministic (user requested random)
cfg = yaml.safe_load(open("config.yaml"))
glob_equil = cfg["globals"]["equilibrium_md"]
warmup_ps = float(glob_equil["warmup_ns"]) * 1000.0
length_ps = float(glob_equil["length_ns"]) * 1000.0
step_ps = float(glob_equil["sample_every_ps"])
n_reps_default = int(cfg["globals"].get("n_reps_default", 1))

def find_final_npt(system_dir):
    # Heuristic: look for a dir/file containing 'npt' and 'final' and having .tpr/.xtc
    cands = []
    for path in glob.glob(os.path.join(system_dir, "**"), recursive=True):
        lp = path.lower()
        if ("npt" in lp and ("final" in lp or "prod" in lp or "production" in lp)):
            cands.append(path)
    tpr = None; xtc = None
    for root in cands:
        if os.path.isdir(root):
            tprs = glob.glob(os.path.join(root, "*.tpr"))
            xtcs = glob.glob(os.path.join(root, "*.xtc"))
            if tprs and xtcs:
                # pick the latest xtc in that folder
                xtc = sorted(xtcs)[-1]
                tpr = sorted(tprs)[-1]
                return tpr, xtc
        else:
            # single files in mixed dirs
            pass
    # fallback: search entire system_dir for any .tpr/.xtc containing 'npt'
    tprs = sorted(glob.glob(os.path.join(system_dir, "**/*npt*final*.tpr"), recursive=True))
    xtcs = sorted(glob.glob(os.path.join(system_dir, "**/*npt*final*.xtc"), recursive=True))
    if not tprs or not xtcs:
        tprs = sorted(glob.glob(os.path.join(system_dir, "**/*.tpr"), recursive=True))
        xtcs = sorted(glob.glob(os.path.join(system_dir, "**/*.xtc"), recursive=True))
    if tprs and xtcs:
        return tprs[-1], xtcs[-1]
    return None, None

os.makedirs("manifests/starts", exist_ok=True)

for sysent in cfg["systems"]:
    system = sysent["name"]
    n_reps = int( sysent.get("variants", [{}])[0].get("n_reps", cfg["globals"].get("n_reps_default",1)) )  # default; variant loop below will override
    sysdir = os.path.join("systems", system)
    for var in sysent.get("variants", []):
        variant = var["id"]
        n_reps_v = int( var.get("n_reps", n_reps_default) )
        tpr, xtc = find_final_npt(sysdir)
        if not tpr or not xtc:
            print(f"[posteq_sample_starts] ERROR: final NPT files not found for {system}/{variant}", file=sys.stderr)
            sys.exit(2)
        # build candidate timestamps
        times = []
        t = warmup_ps
        while t <= length_ps:
            times.append(round(t, 3))
            t += step_ps
        if len(times) == 0:
            print(f"[posteq_sample_starts] ERROR: No sampleable times for {system}/{variant}", file=sys.stderr)
            sys.exit(2)
        if n_reps_v <= len(times):
            selected = random.sample(times, n_reps_v)
        else:
            # sample with replacement if needed
            selected = [random.choice(times) for _ in range(n_reps_v)]
        selected.sort()
        # write CSV
        out = os.path.join("manifests/starts", f"{system}__{variant}.csv")
        with open(out, "w") as f:
            f.write("system,variant,start_time_ps,final_tpr,final_xtc,start_id\n")
            for i, st in enumerate(selected, start=1):
                f.write(f"{system},{variant},{st},{tpr},{xtc},r{i}\n")
        print(f"[posteq_sample_starts] Wrote {out}")
PY
}

py
