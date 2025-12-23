# Fimbriae Steered Molecular Dynamics (SMD) Pipeline

Automated workflow for building, equilibrating, and running steered molecular dynamics (SMD) simulations on fimbrial systems.

## Prerequisites

### Software
- **GROMACS** 2024.4 or newer (CPU and/or GPU builds)
- **Python 3.x** with:
  - PyYAML
  - NumPy (for `make analyze-eq`)
  - Matplotlib (for `make analyze-eq`)
- **VMD** (for movie generation; must be in PATH)
- **GNU Parallel** (optional, for batch trajectory preprocessing)
- **ffmpeg** (for movie encoding)
- **Slurm** workload manager

### Cluster Environment
This pipeline is designed for HPC clusters with Slurm. Adjust partition names, module loads, and resource requests in `config.yaml` to match your site.

---

## Quick Start

### 1. Clone the Repository
```bash
git clone <repo-url>
cd fimbriae-stretch
```

### 2. Configure Your Systems
Edit `config.yaml`:
- **Systems block**: define your PDB path, box size, anchor/pulled residue ranges, speeds, replicates, and performance estimates.
- **GPU vs CPU**: 
  - For **GPU jobs** (recommended for equilibration and fast SMD):
    - Set a GPU-capable partition in `globals.slurm.partition`
    - Uncomment `gpus_per_node` (e.g., `h100:1`)
    - Use the **CUDA-enabled** `gromacs_module` line (comment out the non-CUDA one)
  - For **CPU jobs**:
    - Use a CPU partition
    - Keep `gpus_per_node` commented
    - Use the **non-CUDA** `gromacs_module` line (comment out the CUDA one)

**Template guidance**: The `systems` block is a template; adjust system names, PDB paths, box dimensions, and variant definitions per your needs.

---

## Typical Workflow

### 3. Generate Manifests
```bash
make manifest
```

**What it does**: Scans your `config.yaml` systems block and generates `manifests/systems.csv`, a CSV listing all system/variant combinations that need to be built.

**Why**: This manifest drives the build stage, ensuring each system is prepared consistently before SMD runs.

**Output**: `manifests/systems.csv`

---

### 4. Submit Build Jobs (Equilibration)
**Important**: Equilibration runs best on GPU. Ensure GPU settings are active in `config.yaml` before submitting.

```bash
make build
```

**What it does**: For each system in the manifest, submits a staged Slurm job chain that:
1. **System preparation** (on login node): Adds solvent and ions to your PDB, creating a solvated box.
2. **Energy minimization (EM)**: Removes steric clashes and bad contacts.
3. **NVT equilibration**: Heats the system to target temperature (303 K) with position restraints.
4. **NPT equilibration (short)**: Adds pressure coupling to equilibrate box density.
5. **NPT final**: Runs a longer production-like equilibration (default 10 ns) with no restraints; this trajectory is used for sampling SMD starting frames.

**Why**: Proper equilibration ensures the system reaches a stable state before pulling simulations. Poorly equilibrated systems produce unreliable force profiles.

**Output**: 
- `systems/<system>/00_build/ions.gro`: solvated structure
- `systems/<system>/00_build/npt_final.tpr` and `npt_final.xtc`: final equilibrium trajectory
- `systems/<system>/00_build/topol.top`: system topology
- `systems/<system>/00_build/index.ndx`: GROMACS index groups

**Check queue status**:
```bash
make status
# or directly:
squeue -u $USER
```

---

### 5. Analyze Equilibration
Once equilibration jobs complete, plot energy/RMSD/temperature to verify stability:
```bash
make analyze-eq SYS=fimA_WT
```

**What it does**: Extracts thermodynamic properties (temperature, pressure, density, potential energy) from the equilibration `.edr` files and generates diagnostic plots showing their evolution over time.

**Why**: Visual inspection confirms the system has reached equilibrium. Look for:
- Stable temperature around 303 K
- Converged potential energy (no drift)
- RMSD plateau (structural stability)

**Output**: `systems/<SYS>/00_build/eq_analysis.png`

Inspect the plots to ensure systems are equilibrated before proceeding to SMD.

---

### 6. Sample Starting Frames
Extract starting configurations from equilibrated trajectories:
```bash
make posteq-sample
```

**What it does**: For each system/variant, samples frames from the `npt_final.xtc` trajectory:
- Discards an initial warmup period (default: 2 ns from `equilibrium_md.warmup_ns`)
- Picks frames at regular intervals (default: every 100 ps from `equilibrium_md.sample_every_ps`)
- Writes frame metadata (time, start_id) to `manifests/starts/<system>_<variant>__*.csv`

**Why**: Using multiple diverse starting frames from equilibrium improves statistical sampling of pulling behavior. Each replicate can start from a different equilibrated snapshot, capturing different initial conformations.

**Output**: `manifests/starts/<system>_<variant>__*.csv` with sampled frame times and IDs

---

### 7. Build SMD Manifest
Generate the complete SMD job manifest:
```bash
make smd-manifest
```

**What it does**: 
1. **Reads sampled frames** from `manifests/starts/` (created in step 6)
2. **Generates `manifests/smd_manifest.csv`** with one row per SMD run, combining:
   - Each system/variant
   - Each pull speed
   - Each replicate (multiple independent pulls)
   - Each sampled start frame (diversifies initial conditions)

**Why**: Using multiple starting frames from equilibrium improves statistical sampling of force-extension behavior. Each replicate starts from a different equilibrated snapshot.

**Output**: 
- `manifests/smd_manifest.csv`: complete SMD job list
- `manifests/starts/<system>_<variant>__*.csv`: metadata for sampled frames

---

### 8. Performance Testing (Optional but Recommended)
Run short test jobs (~6 min with `-maxh 0.1`) to gauge performance for different box sizes and tune `perf_ns_per_day` in `config.yaml`.

**GPU test** (lines 2 and 5 from manifest):
```bash
make smd-test-gpu LINES=2,5 CLEAN=1
```

**CPU test** (line 12):
```bash
make smd-test-cpu LINES=12 CLEAN=1
```

**What it does**: Submits short Slurm test jobs that run `mdrun` for ~6 minutes on specific manifest rows. Prints performance metrics (ns/day) to the log.

**Why**: Real-world throughput varies by box size, GPU/CPU hardware, and system composition. Accurate `perf_ns_per_day` estimates prevent walltime underruns (job killed before completion) or overruns (wasted queue time).

`CLEAN=1` removes the test run directory after completion to avoid leftover data.

**Check logs**:
```bash
tail -f logs/smdt_<jobid>.out
# Extract performance:
grep -E "Performance|ns/day" logs/smdt_*.out
```

---

### 8. Adjust Performance Estimates
Based on test results, update `perf_ns_per_day` for each system in `config.yaml`. This ensures accurate walltime requests.

**What it does**: You manually edit `config.yaml` to set realistic throughput values per system based on observed performance from step 8.

**Why**: The submission script calculates Slurm walltime as: `(target_extension_nm / speed_nm_per_ns) / perf_ns_per_day × 24 hours × safety_factor`. Accurate estimates ensure jobs complete successfully without requesting excessive time.

**Example**: If test logs show 85 ns/day for `fimA_WT`, set `perf_ns_per_day: 85` under that system.

---

### 10. Dry-Run SMD Submission
Preview which jobs will be submitted without actually queuing them:

**GPU dry-run**:
```bash
SMD_DEBUG=1 SLURM_TEST_ONLY=1 make smd-submit-new
```

**CPU dry-run**:
```bash
SMD_DEBUG=1 SLURM_TEST_ONLY=1 make smd-submit-new-cpu
```

**What it does**: Parses `manifests/smd_manifest.csv`, groups runs by system and speed, computes walltimes and resource requests, and prints what would be submitted—but does not actually call `sbatch`.

**Why**: Lets you verify job parameters (walltime, partition, GPU/CPU allocation) before flooding the queue. Catch configuration errors early.

Review the printed job details (system, speed, replicate counts, walltime, resources).

---

### 11. Submit SMD Arrays
Once satisfied with the dry-run output:

**GPU mode** (default):
```bash
make smd-submit-new
```

**CPU mode**:
```bash
make smd-submit-new-cpu
```

**What it does**: Submits Slurm array jobs for all NEW (not yet completed) SMD runs. Each array task invokes `scripts/smd_runner.sh`, which:
1. Extracts a starting frame from `npt_final.xtc`
2. Builds run-local topology with anchor position restraints
3. Generates `pull.mdp` with the configured pull rate and spring constant
4. Runs `grompp` to create `pull.tpr`
5. Executes `mdrun` to perform the steered MD pull
6. Writes pull distance/force traces (`pull_pullx.xvg`, `pull_pullf.xvg`)

The submission script:
- Groups runs by system and speed
- Caps concurrent array tasks per `globals.slurm.array_cap`
- Skips already-completed runs (checks for `pull.xtc` and expected `nsteps`)
- Automatically computes walltimes from `perf_ns_per_day` and target extension

**Why**: SMD simulations generate force-extension curves by gradually stretching the system at constant velocity. Each run produces one pulling trajectory; aggregating many replicates across speeds reveals mechanical unfolding pathways and rupture forces.

**Output** (per run): `smd/<system>/<variant>/v<speed>/rep<rep>/<start_id>/`
- `pull.xtc`: SMD trajectory (all atoms, 10 ps stride)
- `pull_pullx.xvg`: pull distance vs. time
- `pull_pullf.xvg`: pull force vs. time
- `run.json`: metadata (speed, rate, nsteps, etc.)

**Monitor progress**:
```bash
make status
# Tail a specific array task log:
tail -f logs/pulls_<array_id>_<task_id>.out
```

---

### 12. Post-Processing: Trajectories and Movies

#### Preprocess Trajectories
Extract protein-only, PBC-fixed, centered trajectories with time-based striding:

**All systems**:
```bash
make preproc-traj J=8
```

**Specific system/variant**:
```bash
make preproc-traj-system SYS=fimA_WT VAR=AtoD J=8
```

**What it does**: For each completed SMD run, uses `gmx trjconv` to:
1. Extract only protein atoms (removes water, ions)
2. Fix periodic boundary conditions (`-pbc nojump`) to avoid molecules jumping across box edges
3. Center the protein in the box
4. Apply time-based striding to reduce file size (e.g., one frame every 1000 ps for a 0.02 nm/ns pull)

`J=<jobs>` controls parallel GNU Parallel job count (default: 4).

**Why**: Raw SMD trajectories include all atoms and may have PBC artifacts. Preprocessing creates clean, compact protein-only trajectories suitable for visualization and further analysis (RMSD, distance measurements, etc.).

**Output**: `analysis/preproc_traj/<system>/<variant>/v<speed>/rep<rep>/<start_id>/`
- `traj_protein_dt<dt>.xtc`: strided protein trajectory
- `start_protein.gro`: initial frame

**Skip behavior**: Already processed runs are skipped (non-destructive).

---

#### Generate Movies
Render MP4 movies from preprocessed trajectories using VMD:

**All movies**:
```bash
make movies
```

**Specific system/variant**:
```bash
make movies-system SYS=fimA_WT VAR=AtoD
```

**What it does**: For each preprocessed trajectory:
1. Loads `start_protein.gro` and `traj_protein_dt*.xtc` into VMD
2. Auto-detects protein chains by residue number resets (or uses a manual chain ranges file)
3. Renders each chain in a different color using NewCartoon representation
4. Saves frames as PPM images (via VMD's `render snapshot`)
5. Encodes frames into MP4 video using ffmpeg (H.264, YUV420p)

**Why**: Movies provide intuitive visual confirmation of pulling dynamics: watch the system stretch, domains unfold, and chains separate. Useful for presentations and spotting anomalies.

**Customize rendering**:
```bash
STRIDE=2 FPS=24 WIDTH=1280 HEIGHT=720 make movies
```
- `STRIDE`: skip frames in the already-strided trajectory (default: 1 = render every frame)
- `FPS`: playback speed
- `WIDTH` / `HEIGHT`: resolution

**Output**: `analysis/movies_out/<system>__<variant>__v<speed>__rep<rep>__<start_id>__traj_protein_dt<dt>.mp4`

**Skip behavior**: Existing non-empty MP4s are skipped.

**Note**: VMD must be in PATH. If not found, the script will prompt you to load the appropriate module.

---

## Adding More Systems/Runs

The pipeline is designed to be **incremental and non-destructive**. You can expand your study at any time by modifying `config.yaml` and re-running the appropriate steps.

### Adding New Systems
1. Add a new system entry to the `systems[]` block in `config.yaml`
2. Run `make manifest` — new system appears in the build manifest
3. Run `make build` — only the new system is built; existing systems are unchanged

### Adding New Variants
1. Add a new variant to an existing system's `variants[]` list in `config.yaml`
2. Run `make smd-manifest` — new variant rows are appended to `manifests/smd_manifest.csv`
3. Run `make smd-submit-new` — only NEW runs are submitted; completed runs are skipped

### Adding More Speeds
1. Add speeds to a variant's `speeds[]` list in `config.yaml`
2. Run `make smd-manifest` — new speed combinations are added to the manifest
3. Run `make smd-submit-new` — submits only the new speed runs

### Adding More Replicates
1. Increase `n_reps` for a variant in `config.yaml`
2. Run `make smd-manifest` — new replicate rows are inserted (preserves existing replicate assignments)
3. Run `make smd-submit-new` — submits the additional replicates

**Key behavior**: 
- `make smd-manifest` preserves existing row assignments and appends/inserts new rows without changing completed runs.
- `make smd-submit-new` checks each run directory for completion (`pull.xtc` + expected steps) and skips already-finished runs.
- You can safely re-run these commands after config changes; they only affect new or incomplete work.

---

## Folder Structure

```
fimbriae-stretch/
├── config.yaml                     # Main configuration (systems, Slurm, mdrun settings)
├── Makefile                        # Convenience targets for the workflow
├── README.md                       # This file
│
├── systems/                        # Per-system data
│   └── <system_name>/
│       ├── <input>.pdb             # Source structure
│       └── 00_build/               # Build/equilibration outputs
│           ├── topol.top           # Processed topology
│           ├── ions.gro            # Solvated, ionized structure
│           ├── npt_final.tpr       # Final NPT run file
│           ├── npt_final.xtc       # Final NPT trajectory (for start sampling)
│           ├── eq_analysis.png     # Equilibration plots (from make analyze-eq)
│           └── index.ndx           # GROMACS index with [Anchor] and [Pulled] groups
│
├── smd/                            # SMD run outputs
│   └── <system>/<variant>/v<speed>/rep<rep>/<start_id>/
│       ├── start.gro               # Starting structure (extracted from NPT)
│       ├── topol.top               # Run-local topology with posre includes
│       ├── pull.mdp                # SMD MDP (generated)
│       ├── pull.tpr                # SMD run file
│       ├── pull.xtc                # SMD trajectory (compressed, all atoms, 10 ps stride)
│       ├── pull_pullx.xvg          # Pull distance trace
│       ├── pull_pullf.xvg          # Pull force trace
│       ├── run.json                # Run metadata (speed, rate, nsteps, etc.)
│       └── expected_nsteps.txt     # Expected step count (for completion checks)
│
├── analysis/                       # Post-processing outputs
│   ├── preproc_traj/               # Preprocessed protein trajectories
│   │   └── <system>/<variant>/v<speed>/rep<rep>/<start_id>/
│   │       ├── traj_protein_dt<dt>.xtc
│   │       └── start_protein.gro
│   └── movies_out/                 # Rendered movies
│       └── <system>__<variant>__v<speed>__rep<rep>__<start_id>__traj_protein_dt<dt>.mp4
│
├── manifests/                      # Job manifests
│   ├── systems.csv                 # Build manifest (from make manifest)
│   ├── smd_manifest.csv            # SMD manifest (from make smd-manifest)
│   └── starts/                     # Sampled start frames metadata
│       └── <system>_<variant>__*.csv
│
├── logs/                           # Slurm output logs
│   ├── build_<system>.out          # Build job logs
│   ├── pulls_<array_id>_<task_id>.out  # SMD array task logs
│   └── smdt_<jobid>.out            # Test run logs (from make smd-test-*)
│
├── scripts/                        # Core automation scripts
│   ├── smd_runner.sh               # SMD job runner (invoked by array tasks)
│   ├── smd_submit_new.sh           # SMD submission logic (NEW runs only)
│   ├── smd_build_manifest.py       # Generate smd_manifest.csv
│   ├── preproc_traj.sh             # Single-run trajectory preprocessing
│   ├── batch_preproc_traj.sh       # Batch preprocessing with GNU Parallel
│   ├── make_movies_vmd.sh          # Movie generation wrapper
│   ├── vmd_make_movie_ppm_chains.tcl  # VMD rendering script
│   ├── smd_test_run.sh             # Short performance test wrapper
│   └── ...
│
├── pipelines/                      # Higher-level orchestration
│   ├── build_preflight.sh          # System prep on login node
│   ├── equil_stage.sh              # Staged equilibration
│   ├── posteq_sample_starts.sh     # Sample start frames from NPT
│   └── analyze_eq.sh               # Equilibration analysis
│
└── templates/                      # MDP templates for build stages
    ├── em.mdp
    ├── nvt.mdp
    ├── npt.mdp
    └── npt_final.mdp
```

---

## Key Configuration Options

### `config.yaml` Structure

#### `globals`
- `dt_ps`: MD time step (ps)
- `k_kj_mol_nm2`: SMD spring constant (kJ/mol/nm²)
- `target_extension_nm`: Default extension target (can be overridden per system/variant)
- `n_reps_default`: Default replicate count
- `pull_speeds_nm_per_ns`: Global list of candidate pull speeds (nm/ns)

#### `globals.slurm`
- `partition`: Slurm partition/queue
- `gpus_per_node`: GPU request (e.g., `h100:1`); **comment out for CPU jobs**
- `cpus_per_task`: CPU threads per task
- `nodes`, `ntasks_per_node`: MPI task distribution (relevant for CPU mode)
- `time_limit`: Max walltime for jobs
- `array_cap`: Concurrent array task limit
- `gromacs_module`: Module load string; **toggle between CUDA and non-CUDA builds**
- `vmd_module` (optional): VMD module hint for movie generation

#### `globals.smd.mdrun`
**GPU-only section**: used for GPU runs; ignored for CPU jobs.
- `ntmpi`, `ntomp`: Thread counts
- `nb`, `bonded`, `pme`, `update`: Device offload (`gpu` for GPU runs)
- `npme`, `ntomp_pme`: PME-specific tuning

#### `systems[]`
Template block; adjust per your systems:
- `name`: System identifier
- `pdb`: Path to source PDB
- `box`: Build box dimensions [x, y, z] in nm
- `perf_ns_per_day`: Estimated SMD throughput (tune after performance tests)
- `SOL_index`: Solvent group index in `index.ndx` (for preprocessing)

#### `systems[].variants[]`
- `id`: Variant identifier
- `anchor`: `{chain, res}` defining anchored residues
- `pulled`: `{chain, res}` defining pulled residues
- `speeds`: Variant-specific pull speeds (nm/ns)
- `n_reps`: Variant-specific replicate count
- `target_extension_nm`: Variant-specific extension target
- `flip_pull_sign` (optional): Set `true` to negate pull direction

---

## Common Tasks

### Check Job Status
```bash
make status
# or:
squeue -u $USER
```

### Tail Logs
```bash
# Build log for a system:
make tail-build SYS=fimA_WT

# SMD array task log:
make tail-pull A=<array_id> a=<task_id>

# Test run log (latest):
tail -f $(ls -t logs/smdt_*.out | head -n 1)
```

### Cancel Jobs
```bash
make cancel A=<job_or_array_id>
# or:
scancel <job_id>
```

### Clean Ledger
Remove completed runs from the submission ledger:
```bash
make smd-clean-ledger
```

### Reset System
Wipe manifests and outputs for a single system:
```bash
make reset SYS=fimA_WT
```

### Reset All
Wipe all manifests and outputs:
```bash
make reset-all
```

---

## Troubleshooting

### Equilibration Issues
- **System not equilibrated**: Check plots from `make analyze-eq`. If energy/RMSD are unstable, extend `equilibrium_md.length_ns` in `config.yaml` or adjust restraint forces.
- **Build failures**: Check `logs/build_<system>.out` for grompp/mdrun errors.

### SMD Submission Issues
- **No jobs submitted**: Ensure `manifests/smd_manifest.csv` exists and has rows with `status=NEW`. Run `make smd-manifest` if needed.
- **Walltime too short**: Tune `perf_ns_per_day` in `config.yaml` based on test runs.
- **GPU/CPU mismatch**: Verify `gromacs_module` and `gpus_per_node` settings match your intended run mode.

### Preprocessing Issues
- **GROMACS module not loaded**: Ensure `globals.slurm.gromacs_module` is set correctly. The preprocessing script auto-loads it if available.
- **trjconv group errors**: Adjust `PROT_GROUP` env var if your index numbering differs (default: `1` for Protein).

### Movie Generation Issues
- **VMD not found**: Add VMD to PATH or load the appropriate module before running `make movies`. Optionally set `globals.slurm.vmd_module` in `config.yaml` for a helpful prompt.
- **ffmpeg errors**: Ensure ffmpeg is installed and in PATH.
- **Headless rendering**: On compute nodes without X11, the script uses `xvfb-run` automatically if available.

---

## Performance Tuning

### Estimating `perf_ns_per_day`
1. Run short tests: `make smd-test-gpu LINES=2,5 CLEAN=1`
2. Extract performance from logs: `grep "Performance" logs/smdt_*.out`
3. Update `perf_ns_per_day` for each system in `config.yaml` to match observed throughput.
4. Re-run dry-run to verify walltime estimates: `SMD_DEBUG=1 SLURM_TEST_ONLY=1 make smd-submit-new`

### GPU vs CPU Trade-offs
- **GPU**: 10-50× faster for typical box sizes; recommended for equilibration and SMD.
- **CPU**: Useful when GPU partitions are congested or for very large systems where CPU scaling compensates.
- **Hybrid**: Run equilibration on GPU, then switch to CPU for SMD if needed (update `config.yaml` accordingly).

---

## Advanced Options

### Environment Variable Overrides

**Preprocessing**:
- `PREPROC_JOBS=8`: Control parallel job count
- `DX_NM=0.02`: Target extension per frame (nm)
- `PROT_GROUP=1`: Protein group index for `gmx trjconv`
- `GMX=gmx_mpi`: Override GROMACS command

**Movie generation**:
- `STRIDE=2`: Frame sampling step (default: 1 = every frame)
- `FPS=30`: Movie playback speed (frames per second)
- `WIDTH=1600 HEIGHT=900`: Render resolution
- `RANGES_FILE=/path/to/chain_ranges.txt`: Manual chain range specification

**SMD runner**:
- `DRY_RUN=1`: Skip mdrun (grompp only; useful for testing)
- `MAXH_HOURS=<float>`: Override `-maxh` for mdrun

---

## Make Targets Reference

### Build & Equilibration
- `make manifest`: Generate build manifest
- `make prep SYS=<name>`: Prepare system on login node (up to ions.gro)
- `make build`: Submit staged build chains (GPU recommended)
- `make analyze-eq SYS=<name>`: Plot equilibration metrics

### SMD
- `make smd-manifest`: Generate SMD manifest with sampled starts
- `make smd-dry-run-one`: Test first manifest row locally (grompp only)
- `make smd-test-gpu LINES=<N[,N2,...]> [CLEAN=1]`: Short GPU performance test
- `make smd-test-cpu LINES=<N[,N2,...]> [CLEAN=1]`: Short CPU performance test
- `make smd-submit-new`: Submit NEW SMD runs (GPU mode)
- `make smd-submit-new-cpu`: Submit NEW SMD runs (CPU mode)
- `make smd-clean-ledger`: Remove completed runs from ledger

### Post-Processing
- `make preproc-traj [J=<jobs>]`: Preprocess all trajectories in parallel
- `make preproc-traj-system SYS=<name> [VAR=<var>] [J=<jobs>]`: Preprocess specific system/variant
- `make movies`: Generate all movies
- `make movies-system SYS=<name> [VAR=<var>]`: Generate movies for specific system/variant

### Utilities
- `make status`: Show your Slurm queue
- `make tail-build SYS=<name>`: Tail build log
- `make tail-pull A=<array_id> a=<task_id>`: Tail SMD array task log
- `make cancel A=<job_id>`: Cancel Slurm job
- `make reset SYS=<name>`: Reset one system
- `make reset-all`: Reset all systems

---

## Example End-to-End Session

```bash
# 1. Clone and configure
git clone <repo-url>
cd fimbriae-stretch
# Edit config.yaml: set systems, GPU partition, CUDA module

# 2. Build and equilibrate
make manifest
make build
make status  # wait for completion

# 3. Verify equilibration
make analyze-eq SYS=fimA_WT
make analyze-eq SYS=mfa1_WT

# 4. Sample starting frames
make posteq-sample

# 5. Generate SMD manifest
make smd-manifest

# 6. Performance test (GPU)
make smd-test-gpu LINES=2,5 CLEAN=1
tail -f $(ls -t logs/smdt_*.out | head -n 1)
grep "Performance" logs/smdt_*.out
# Update perf_ns_per_day in config.yaml

# 7. Dry-run SMD submission
SMD_DEBUG=1 SLURM_TEST_ONLY=1 make smd-submit-new

# 8. Submit SMD arrays
make smd-submit-new
make status

# 9. Wait for completion, then preprocess
make preproc-traj J=8

# 10. Generate movies
make movies

# 11. Enjoy your fimbrial dynamics!
```

---

## Contributing

Improvements, bug reports, and feature requests are welcome. Open an issue or submit a pull request.

---

## Contact

For questions or issues, contact the repository maintainer.
