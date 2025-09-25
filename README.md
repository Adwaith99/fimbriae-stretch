# fimA-project (v1)

Minimal, reproducible pipeline for GROMACS SMD pulls on HPC:
- Build/equilibration via Slurm.
- Two-step index: selections from `clean.pdb`, rename, import onto solvated system, append `NonProtein`.
- Local-numbered posre inside the molecule `.itp`.
- Self-requeue Slurm array with `%` cap; single-file `-append` continuation.

## Quick start
1) Edit `config.yaml` (systems, variants, speeds).
2) Commit & push to your remote (or use rsync deploy).
3) On cluster:
   - `bash scripts/submit_builds.sh`  (builds any new systems)
   - once builds finish → `bash scripts/submit_arrays.sh`  (submits pulls with array % cap)
4) Monitor with `squeue` or `make status`.
5) Pull small artifacts back: `make pull-manifests`.

See inline comments in each script. This is v1 — we’ll iterate.
