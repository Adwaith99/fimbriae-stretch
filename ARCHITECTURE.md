# Project Architecture & Pipeline Overview

This document captures a deep-dive understanding of the current `fimbriae-stretch` repository based on code inspection (Oct 2025). It summarizes workflow stages, data / file flows, responsibilities of each script, and proposes concrete improvements.

---
## 1. High-Level Purpose
Automated preparation, staged equilibration, and steered molecular dynamics (SMD) pulling of fimbrial protein systems in GROMACS on a SLURM GPU cluster. A single `config.yaml` + generated CSV manifests drive both build and pulling phases.

---
## 2. Core Concepts & Entities
| Concept | Description | Primary Source |
|---------|-------------|----------------|
| System (`systems/<SYS>`) | A protein construct with its PDB and build outputs. | `config.yaml` (systems list) |
| Variant | Logical pulling geometry (anchor/pulled selection). | `config.yaml` under a system's `variants` |
| Manifest row | A pulling job replicate (speed × k × rep × variant). | `manifests/manifest.csv` |
| Build chain | EM → (NVT,NPT)* at descending restraint force constants → final long NPT. | `pipelines/fimA_build_submit_staged.sh`, `pipelines/equil_stage.sh` |
| Equilibrium MD | Final unrestrained NPT (length configurable). | `config.yaml: globals.equilibrium_md` |
| Pull start frames | Subsampled frames from final NPT for replicate starts. | `pipelines/posteq_sample_starts.sh` |
| Pull template | Parameterized SMD `.mdp` file. | `templates/pull.tpl.mdp` |
| Index + local posre | Generated index containing Anchor/Pulled and per-anchor local posres. | `scripts/build_indices_and_posres.py` |

---
## 3. Configuration Model (`config.yaml`)
Key Blocks:
* `globals.dt_ps` — Base integration timestep (used for final equilibrium NPT and pulling). Staged equilibration uses 0.001 ps unless final.
* `globals.k_kj_mol_nm2` — Default pulling spring constant (k).
* `globals.pull_speeds_nm_per_ns` — Default velocity set for manifest generation.
* `globals.n_reps_default` — Default replicate count per (variant, speed) unless overridden per system.
* `globals.target_extension_nm` — Used by runner to decide DONE vs requeue for SMD.
* `globals.slurm.*` — Cluster resource template + GROMACS module string.
* `globals.posre.force` — Base position restraint triplet (scaled through FC schedule externally — implicit, not directly enumerated yet for build). 
* `globals.equilibrium_md` — Parameters for final NPT: length, warmup discard, sampling cadence, write frequencies, walltime budgets.
* System entry:
  - `name`, `pdb`, `box` (vector in nm), `anchor_molecule_itp` (for future per-chain posre inclusion), `variants` (list). 
  - Variant: `id`, and anchor/pulled residue range definitions currently structured but not yet consumed by index builder (which instead expects `anchor_pattern` / `pulled_pattern`— see Gap #1 below).
  - Optional overrides: `speeds`, `n_reps`.

Observed Gap (#1): `build_indices_and_posres.py` looks for `anchor_pattern` / `pulled_pattern` keys, but config uses nested objects `{ chain: "D", res: "135-150" }`. A translation or adaptation layer is missing; current script would fail without those pattern keys.

---
## 4. Build & Equilibration Workflow
Initiated by `make build` → `scripts/submit_builds.sh` → `pipelines/fimA_build_submit_staged.sh` which:
1. Seeds `systems/<SYS>/00_build/input.pdb`.
2. Submits EM job.
3. Iteratively submits ordered dependency chain of paired NVT then NPT stages over force constants: 1000, 500, 200, 100, 50 kJ/mol/nm² (naming `nvt_fcXXX`, `npt_fcXXX`). Temperatures determined via `temp_for_fc()` mapping (override capable via `globals.build.temp_schedule_by_fc`).
4. Submits final long NPT (`npt_final`), determining total steps from length_ns / dt.

`equil_stage.sh` handles each stage:
* Chooses previous structure (fallback precedence list).
* Dynamically computes thread / domain decomposition based on environment and config fallback.
* Creates `stage.mdp` with stride logic: default coarse outputs; final NPT strides recalculated from `equilibrium_md` (xtc_write_ps etc.).
* Copies final unrestrained NPT artifacts to canonical `npt_final.*` and marks `BUILD_DONE`.

Energy extraction and QC:
* `make analyze-eq SYS=<SYS>` → `pipelines/analyze_eq.sh`: loads GROMACS module, runs `scripts/extract_eq.sh` to pick fixed indices (Temperature 16, Pressure 17, Density 23 in staged, Pressure 16 & Density 22 for final). Then `scripts/plot_eq_grid.py` composes grid figure.

Sampling start frames for SMD:
* `pipelines/posteq_sample_starts.sh` slices `npt_final.xtc` after warmup, every `sample_every_ps` ps, linking first N replicate frames as `repNNN.gro` under `_starts/`.

---
## 5. Manifest Generation & Status Tracking
`scripts/generate_manifest.py` composes two CSVs:
* `manifests/systems.csv` — Build status (presence of `BUILD_DONE`).
* `manifests/manifest.csv` — Cartesian product of (system variants × speed choices × replicates) with persistent status (PENDING|DONE|RETIRED). UID structure: `system|variant|vSPEED|kK|repNNN`.

`scripts/submit_arrays.sh` filters pending rows and launches a single SLURM array (`pipelines/pulls_array_submit.sh`). The array size equals number of PENDING rows (capped externally only by `%1` concurrency currently— effectively serializing pulls; can be tuned).

---
## 6. Pulling (SMD) Workflow
Per array task, `scripts/runner.sh`:
1. Resolves nth PENDING manifest row (independent of original manifest order – it re-filters). Robust trimming of whitespace & CRLF.
2. Ensures build artifacts and (variant) index exist; lazily generates index & local anchor posres via `build_indices_and_posres.py`.
3. Determines run directory layout: `systems/<SYS>/20_pulls/<variant>/vSPEED/kK/repNNN`.
4. Chooses start structure: sampled frame or fallback final NPT .gro.
5. Computes COM distance between Anchor & Pulled groups with `gmx distance` and sets pull direction sign on x-axis.
6. Renders `pull.mdp` from `templates/pull.tpl.mdp` performing simple placeholder substitution.
7. Appends local anchor posre include to a copied `topol.top`.
8. Runs `gmx grompp` + `gmx mdrun` with dynamic CPU/GPU layout.
9. Monitors extension progress via `pullx.xvg` last data record; updates manifest status to DONE once ≥ target extension or requeues otherwise.

Idempotency & resilience:
* Requeue loop continues until target extension reached.
* Fallback behavior if COM distance extraction fails (init=0, positive x) avoids hard failure.

---
## 7. Index & Posre Generation Details
`build_indices_and_posres.py` expects variant selection patterns (string patterns referencing groups) and performs fuzzy normalization to assign `Anchor` and `Pulled`. It writes:
* `<variant>_pull.ndx` (with standardized group names) under `indices/<system>/`.
* `<variant>_posre_anchor_local.itp` with per-atom restraint entries at configurable force constant.

Potential mismatch: Current config variant definition uses structured `{ chain: "D", res: "135-150" }` objects rather than a textual pattern; a translation layer is needed (see Improvements #2).

---
## 8. Energy & Plotting Utilities
* `extract_eq.sh` hardcodes energy menu indices (could be brittle across GROMACS versions). A supplemental `detect_energy_index.py` exists but is not plumbed into the extractor yet.
* `plot_eq_grid.py` auto-grids by (Temperature, FC), optional EM row spanning columns, seasoned for multi-stage visual QA.

---
## 9. Resource / Performance Strategy
* Threading scheme heuristics in `equil_stage.sh` and `runner.sh`: derive GPU count from SLURM env or config fallback; assign `ntmpi = gpus * 2` then compute `ntomp = cpus / ntmpi`. This is a reasonable starting heuristic but may oversubscribe for small systems.
* GPU directives use full GPU offload: `-nb gpu -bonded gpu -pme gpu -update gpu` + `GMX_ENABLE_DIRECT_GPU_COMM=1` (modern 2024.x best practice).

---
## 10. File & Directory Lifecycle
```
systems/<SYS>/00_build/
  input.pdb → clean.pdb → aligned.pdb → boxed.gro → solv.gro → ions.gro
  em.* → nvt_fcXXXX.* → npt_fcXXXX.* → npt_final.*
  BUILD_DONE (flag)
systems/<SYS>/20_pulls/
  _starts/repNNN.gro (links)
  <variant>/vSPEED/kK/repNNN/(pull.* outputs, pullx.xvg, etc.)
indices/<SYS>/<variant>_pull.ndx, <variant>_posre_anchor_local.itp
manifests/manifest.csv (status evolves)
logs/*.out (SLURM output)
```

---
## 11. Error Handling & Robustness Observations
| Area | Current Behavior | Risk |
|------|------------------|------|
Config validation | Minimal (only presence of `globals` & `systems`) | Silent schema drift / mismatch |
Energy index extraction | Hardcoded indices | GROMACS index ordering shifts break analysis |
Variant pattern resolution | Assumes text patterns; config provides structured objects | Indices not generated → pulling blocked |
Manifest concurrency | `%1` array throttle | Under-utilization of cluster resources |
Requeue termination | Based only on extension; no max wall-clock / retry cap | Infinite loops if extension stagnates |
Module load | Single module string; fallback logic attempts alternate load | Site module changes could break pipeline |

---
## 12. Security / Reproducibility Notes
* No pinning of Python dependencies (no `requirements.txt`).
* Relies on externally available GROMACS module (not containerized) → environment drift risk.
* No provenance capture (e.g., commit hash) embedded into output logs.

---
## 13. Proposed Improvements (Prioritized)
Priority legend: [H]=High impact / low-medium effort; [M]=Medium; [L]=Nice-to-have.

1. [H] Align variant selection config: Add auto synthesis of pattern strings from `{chain,res}` objects or update `config.yaml` spec & README. Provide a helper function used by both build and pulling scripts.
2. [H] Expand config validation: Schema validation (e.g. `cerberus` or `jsonschema`) ensuring required keys, numeric ranges, and variant definitions.
3. [H] Integrate dynamic energy index detection: Optionally call `detect_energy_index.py` inside `extract_eq.sh` when a fixed index extraction fails; fallback to hardcoded indices for speed.
4. [H] Add retry/backoff cap for pulling requeue: e.g., annotate a counter file; abort + mark row FAILED after N cycles without progress.
5. [M] Parallel pulling: Increase array concurrency (e.g., `%${ARRAY_CAP}` from config) and implement per-row resource shaping if needed (maybe 1 GPU per task via `--gpus-per-node=1`).
6. [M] Introduce Python dependency manifest + optional virtualenv bootstrap script.
7. [M] Automatically embed Git metadata (commit SHA, dirty flag) into each SLURM job log header for reproducibility.
8. [M] Containerization option (Singularity/Apptainer) for GROMACS + Python to decouple from site modules.
9. [M] Central logging utility: replace ad-hoc `echo` with a small Python logger that timestamps + context-tags lines.
10. [M] Unit tests for pure-Python utilities (index normalization, posre patching) using `pytest` + synthetic inputs.
11. [M] Add `ARCHITECTURE.md` (this file) to README table-of-contents.
12. [L] Performance heuristic refinement: Query system size (`gmx check -s`) to scale `ntmpi` / `ntomp` adaptively.
13. [L] Add automatic frame selection heuristic (e.g., evenly spaced after RMSD stabilization) instead of fixed time spacing.
14. [L] Generate lightweight HTML report summarizing build & early QC (embed PNG grid + key stats).
15. [L] Optional Slack / email notifications when final NPT or all pulls finish.
16. [L] Pre-flight environment probe script: verifies module availability, Python libs, disk quota, GPU visibility.
17. [L] Failure classification: Tag manifest rows FAILED vs DONE (distinct from RETIRED) enabling targeted reruns.

---
## 14. Immediate Low-Risk Quick Wins
1. Add schema validation + config pattern translation layer.
2. Introduce dynamic energy term fallback.
3. Add `ARRAY_CAP` usage for pulling arrays (reuse `globals.slurm.array_cap`).
4. Generate `requirements.txt` with pinned minimal versions (`PyYAML`, `matplotlib`, `numpy`).
5. Inject git commit hash into `submit_builds.sh` and `runner.sh` logs (e.g., `echo "[GIT] $(git rev-parse --short HEAD)"`).

---
## 15. Suggested Future Refactor (Medium Scope)
Create a unified Python CLI (`fimpipe.py`) exposing subcommands: `manifest`, `prep`, `submit-build`, `analyze-eq`, `sample-starts`, `submit-pulls`, `resume-pulls`, centralizing YAML parsing and eliminating duplication in bash. Gradually phase bash scripts into thin wrappers calling the CLI.

---
## 16. Risks & Mitigations Summary
| Risk | Mitigation |
|------|-----------|
Hardcoded energy indices drift | Implement auto-detect fallback + version stamp in logs |
Variant config mismatch | Add translation or enforce textual pattern keys | 
Indefinite pull requeues | Add max retry & extension stagnation threshold |
Environment drift (modules) | Optional container / recorded module snapshot |
Silent partial failures | Add explicit DONE/FAILED markers & manifest status transitions |

---
## 17. Glossary
* FC — Position restraint force constant applied during staged equilibration.
* SMD — Steered Molecular Dynamics; constant-velocity pulling between two groups.
* Anchor / Pulled — Defined atom selections whose COM separation is biased.
* Manifest — Tabular declaration of (variant, speed, rep) pulling tasks and their lifecycle status.

---
## 18. Observed Inconsistencies / Action Items
| Issue | File(s) | Action |
|-------|---------|--------|
Variant pattern vs structured definitions | `config.yaml`, `build_indices_and_posres.py` | Implement transformation or update config schema |
Unused helper script(s) (e.g. `rename_ndx_by_match.py`) | `scripts/rename_ndx_by_match.py` | Either integrate or remove to reduce surface area |
Duplicate shebang in heredoc (double `#!/bin/bash`) | `pipelines/fimA_build_submit_staged.sh` heredoc | Remove duplicate (cosmetic) |

---
## 19. Conclusion
The repository provides a solid, modular shell-driven orchestration around GROMACS with good idempotency and clear artifact naming. Addressing config schema alignment, dynamic energy index fallback, and improved status/error semantics will materially enhance robustness and maintainability. Transitioning repetitive shell/YAML parsing logic into a consolidated Python CLI will reduce drift and speed feature iteration.

*Document generated via repository introspection (no external runtime execution).* 

---
## 20. Addendum (Recent Enhancements)
Date: Post-initial audit (Oct 2025)

Implemented changes:
1. Structured variant selector support: `build_indices_and_posres.py` now auto-synthesizes textual patterns from objects like:
   ```yaml
   variants:
     - id: AtoD
       anchor: { chain: "D", res: "135-150" }
       pulled: { chain: "A", res: "285-300" }
   ```
   No need to add `anchor_pattern` / `pulled_pattern` manually (still accepted if present).
2. Per-system pull arrays: `scripts/submit_arrays.sh` now splits pending manifest rows by system and submits one SLURM array per system (job name `pulls_<SYS>`). This enables distinct walltime estimation and clearer monitoring/cancellation.
3. Dynamic walltime estimation: For each system array submission we approximate required walltime from:
   - `target_extension_nm`
   - slowest speed among that system's pending rows
   - optional performance mapping: `pulling.performance_ns_per_day: { <system>: value }`
   - integration timestep `dt_ps`
   Formula (ns needed / ns_per_day * 1.25 safety) → rounded up hours (capped 168h).
4. Array concurrency: Uses `globals.slurm.array_cap` to set the `%CAP` throttle per system instead of hard-coded `%1`.

Additional simplification (latest):
5. Removed SMD self-requeue logic. Each replicate runs once; manifest status becomes `DONE` if target extension reached during that single run, else `PARTIAL`. This avoids indefinite recycling and keeps bookkeeping minimal. Re-runs (if desired) can be triggered by resetting the row status back to `PENDING` manually or via a future helper script.

Global timing keys now supported (preferred over older pulling.performance map):
| Key | Meaning |
|-----|---------|
| `globals.perf_ns_per_day` | Estimated pulling throughput (ns/day) used for walltime prediction. |
| `globals.walltime_safety` | Multiplicative safety factor on predicted walltime. |
| `globals.max_wall_hours`  | Upper bound cap on requested walltime hours per system array. |

Per-system override precedence for throughput: `system.perf_ns_per_day` > `system.pulling.perf_ns_per_day` > `globals.perf_ns_per_day` > default 10.0.

Speeds & replicate precedence (manifest generation):
`variant.speeds` / `variant.n_reps` > `system.speeds` / `system.n_reps` > `globals.pull_speeds_nm_per_ns` / `globals.n_reps_default`.

Optional new config block (illustrative):
```yaml
pulling:
  performance_ns_per_day:
    fimA_WT: 15.0   # Empirical pulling throughput for walltime estimation
```

No changes were made to the equilibration/build pipeline per optimization request.
