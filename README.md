# fimbriae-stretch: Automated GROMACS Pipeline for Fimbrial Protein Simulations

This repository contains an automated, modular pipeline for running **protein stretching (pulling) molecular dynamics simulations** in **GROMACS**.  
It supports complete system preparation, staged equilibration, and initial quality-control analysis.  
All steps are configured by a single `config.yaml` file and executed using a `Makefile` interface.

---

## Project Overview

Each system (e.g., `fimA_WT`) has its own subdirectory under `systems/`, containing the PDB and build artifacts.  
All jobs are run on a Compute Canada–style SLURM cluster with GROMACS 2024.4 and CUDA 12.2.

The pipeline handles:

| Stage | Description | Job type |
|:------|:-------------|:---------|
| EM | Energy minimization of the solvated, ionized system | separate SLURM job |
| NVT / NPT (staged) | Five equilibration pairs with decreasing position restraints and increasing temperature | sequential SLURM jobs |
| Final NPT | Unrestrained equilibrium MD at 303 K and 1 bar (production-like) | single SLURM job |
| Analysis | Extracts temperature, pressure, and density traces and generates a grid plot (PDF + PNG) | local |

After equilibration, snapshots from the final NPT run will be used as starting configurations for **steered MD pulling simulations** (next development stage).

---

## Directory Layout

```
fimbriae-stretch/
│
├── config.yaml               # master configuration (temperatures, restraints, SLURM settings)
├── Makefile                  # main task interface
│
├── pipelines/                # main workflow scripts
│   ├── fimA_build_submit_staged.sh   # submits staged equilibration jobs
│   ├── equil_stage.sh                # generates MDPs and runs grompp/mdrun
│   ├── analyze_eq.sh                 # energy extraction + plotting wrapper
│
├── scripts/                  # helper utilities
│   ├── extract_eq.sh         # extract Temperature / Pressure / Density / Potential from .edr
│   ├── detect_energy_index.py # optional helper for automatic index detection
│   ├── plot_eq_grid.py       # produces eq_grid.png/pdf
│
├── manifests/                # auto-generated run manifests
│
└── systems/
    └── fimA_WT/
        ├── nfimA_5chains_beta_deleted.pdb
        └── 00_build/
            ├── em.edr, nvt_fc*.edr, npt_fc*.edr, npt.edr
            ├── topol.top, posre_*.itp, ...
            └── xvg/               # analysis outputs
```

---

## Setup and Usage

### 1. Clone and sync
```bash
git clone git@github.com:<USER>/fimbriae-stretch.git
cd fimbriae-stretch
```

Edit `config.yaml` to define:
- Global simulation parameters (`dt_ps`, restraint constants, number of replicates)
- SLURM resource parameters (`partition`, `gpus_per_node`, `cpus_per_task`, `time_limit`)
- Equilibration length and walltimes under `globals.equilibrium_md`

Commit and push to synchronize between your laptop and the cluster.

---

### 2. Build and equilibrate systems

#### a) Prepare input files
Copy each PDB into `systems/<SYS>/` and set its path in `config.yaml`.

#### b) Submit equilibration
On the **cluster login node**:
```bash
make build
```

- Each stage (EM, NVT/NPT pairs, final NPT) is submitted automatically.  
- If a job completes successfully, it will be skipped in subsequent runs.  
- If interrupted, rerun `make build` — the pipeline resumes from the last finished step.

Logs are written to `logs/`.

---

### 3. Analyze equilibration quality

Run locally or on the cluster:
```bash
make analyze-eq SYS=fimA_WT
```

This runs:
1. `scripts/extract_eq.sh` → extracts `T`, `P`, `ρ`, and `Potential` from each `.edr`.
   - Staged NVT/NPT use fixed menu indices:  
     - Temperature = 16  
     - Pressure = 17  
     - Density = 23  
   - Final NPT (`npt.edr` or `npt_final.edr`) uses:  
     - Pressure = 16  
     - Density = 22  
   - EM potential uses index = 11
2. `scripts/plot_eq_grid.py` → generates `eq_grid.png` and `eq_grid.pdf` in `00_build/`.

Each row corresponds to one temperature/restraint stage, with Temperature, Pressure, and Density plotted side-by-side.  
The EM potential appears as a single top-panel trace.

---

## Job Management Notes

- **Re-runs:** Each stage checks for its output `.gro`; if found, it is skipped.  
- **Final NPT:** Currently runs as one job; walltime and length are defined in `config.yaml`.  
- **Checkpointing:** GROMACS writes periodic `.cpt` files; you can resume manually by re-submitting the same stage.

---

## Default Output Frequencies

(Coarse for lightweight equilibration)

| Quantity | MDP key | Default interval |
|-----------|----------|------------------|
| Energies | `nstenergy` | 50 000 steps (100 ps) |
| Trajectory (XTC) | `nstxout-compressed` | 10 000 steps (20 ps) |
| Logs | `nstlog` | 50 000 steps (100 ps) |

These can be tightened in `equil_stage.sh` if smoother plots or higher-resolution trajectories are desired.

---

## Next Development Stage: Pulling MD (SMD)

The next step is implementing **steered MD pulling** from equilibrated configurations.  
Upcoming features:
- Automatic selection of starting frames from the final NPT trajectory  
- Generation of pulling `.mdp` files with specified spring constants and velocities  
- Batch submission of replicate SMD runs via SLURM arrays  
- Integrated pulling-force and extension analysis

---

## Troubleshooting and Tips

- To rerun a single system’s equilibration from scratch, remove its `00_build/BUILD_DONE` and stage `.gro` files, then run `make build`.  
- To reset all manifests:
  ```bash
  make clean
  ```
- If `make analyze-eq` only produces one curve, check `.edr` file sizes — empty `.edr` means that stage didn’t finish properly.

---

## Dependencies

- **GROMACS 2024.4** (CUDA build recommended)  
- **Python ≥ 3.8** with:  
  - `PyYAML`  
  - `matplotlib`  
  - `numpy`  
- SLURM job scheduler

---

## Citation

If you use or adapt this pipeline for publications, please cite GROMACS and note this workflow as  
“Automated staged equilibration and pulling simulation pipeline (fimbriae-stretch)”.
