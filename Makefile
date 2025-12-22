SHELL := /bin/bash
ROOT  := $(shell pwd)

help:
	@echo "Targets (run on the cluster login node):"
	@echo "  pull                 - git pull (rebase) latest code"
	@echo "  validate             - sanity check config.yaml"
	@echo "  manifest             - (re)generate manifests"
	@echo "  prep SYS=...         - run system preparation on login node (up to ions.gro)"
	@echo "  prep-em SYS=...      - prep + EM grompp compile (no run)"
	@echo "  build                - submit staged build chain (queued)"
	@echo "  pulls                - submit SMD array over PENDING rows"
	@echo "  status               - show your Slurm queue"
	@echo "  tail-build SYS=...   - tail build log for a system"
	@echo "  tail-pull A=ID a=TID - tail a pull array task log"
	@echo "  cancel A=ID          - cancel a Slurm job/array"
	@echo "  reset SYS=...        - wipe manifests + build outputs for one system"
	@echo "  reset-all            - wipe manifests + build outputs for all systems"
	@echo "  smd-submit           - submit SMD arrays (GPU mode, default)"
	@echo "  smd-submit-cpu       - submit SMD arrays (CPU mode)"
	@echo "  smd-submit-new       - submit only NEW SMD rows (GPU mode)"
	@echo "  smd-submit-new-cpu   - submit only NEW SMD rows (CPU mode)"
	@echo "  smd-clean-ledger     - remove completed runs from smd_submitted.csv"
	@echo ""
	@echo "Note: smd-submit-new ignores the ledger by default (only skips completed runs)."
	@echo "      Cancel running jobs manually before resubmitting (e.g., scancel <jobid>)."

pull:
	git pull --rebase

validate:
	python3 scripts/validate_config.py

manifest:
	python3 scripts/generate_manifest.py

# --- Prep-only on login node ---
prep:
	@if [ -z "$(SYS)" ]; then echo "Usage: make prep SYS=<SYSTEM_NAME>"; exit 2; fi
	bash pipelines/build_preflight.sh $(SYS)

prep-em:
	@if [ -z "$(SYS)" ]; then echo "Usage: make prep-em SYS=<SYSTEM_NAME>"; exit 2; fi
	bash pipelines/build_preflight.sh $(SYS) --em-grompp

# --- Submit queued jobs ---
build: manifest
	bash scripts/submit_builds.sh

# --- Submit SMD array jobs over PENDING rows ---
pulls: manifest
	bash scripts/submit_arrays.sh

status:
	squeue -u $$USER | egrep "JOBID|fimA|build_|pulls_" || true

tail-build:
	@if [ -z "$(SYS)" ]; then echo "Usage: make tail-build SYS=fimA_WT"; exit 2; fi
	@tail -f logs/build_$(SYS).out

tail-pull:
	@if [ -z "$(A)" ] || [ -z "$(a)" ]; then echo "Usage: make tail-pull A=<array_id> a=<task_id>"; exit 2; fi
	@tail -f logs/pulls_$(A)_$(a).out

cancel:
	@if [ -z "$(A)" ]; then echo "Usage: make cancel A=<job_or_array_id>"; exit 2; fi
	scancel $(A)

# --- Reset helpers ---
reset:
	@if [ -z "$(SYS)" ]; then echo "Usage: make reset SYS=<SYSTEM_NAME>"; exit 2; fi
	@echo "Resetting system: $(SYS)"
	@rm -f manifests/*.csv
	@rm -rf systems/$(SYS)/00_build systems/$(SYS)/20_pulls indices/$(SYS)

reset-all:
	@echo "Resetting ALL systems + manifests"
	@rm -f manifests/*.csv
	@for d in systems/*; do \
		if [ -d "$$d" ]; then \
			rm -rf "$$d/00_build" "$$d/20_pulls"; \
		fi; \
	done
	@rm -rf indices/*

.PHONY: analyze-eq
analyze-eq:
	bash pipelines/analyze_eq.sh $(SYS)

.PHONY: posteq-sample smd-manifest smd-submit smd-submit-cpu smd-qc

# 1) Sample start frames from final NPT (per system/variant), random picks per variant
posteq-sample:
	@bash pipelines/posteq_sample_starts.sh

# 2) Build SMD manifest (systems × variants × speeds × replicates; each replicate gets one start)
smd-manifest:
	@python3 scripts/smd_build_manifest.py

# 3) Submit one Slurm array per system (capped by globals.slurm.array_cap), no requeue
smd-submit: smd-manifest
	@bash pipelines/smd_array_submit.sh

# 3b) Submit one Slurm array per system (CPU mode)
smd-submit-cpu: smd-manifest
	@CPU=1 bash pipelines/smd_array_submit.sh

# 4) Optional QC pass over completed runs (TSV and PNG if enabled)
smd-qc:
	@python3 scripts/smd_qc.py

.PHONY: smd-dry-run-one
# Runs the *first* manifest row locally in DRY_RUN=1 mode: exercises sampling, manifest,
# direction/init computation, index discovery, and grompp — but skips mdrun.
smd-dry-run-one: smd-manifest
	@LINE="$$(awk -F, 'NR==2{print $$0}' manifests/smd_manifest.csv)"; \
	if [ -z "$$LINE" ]; then echo "[smd-dry-run] No rows in manifests/smd_manifest.csv"; exit 2; fi; \
	echo "[smd-dry-run] Testing with: $$LINE"; \
	module purge; module load $$(python3 -c 'import yaml;print(yaml.safe_load(open("config.yaml"))["globals"]["slurm"]["gromacs_module"])'); \
	DRY_RUN=1 bash scripts/smd_runner.sh "$$LINE"

.PHONY: indices-test
indices-test:
	@sys=$$(awk -F, 'NR==2{print $$1}' manifests/smd_manifest.csv); \
	var=$$(awk -F, 'NR==2{print $$2}' manifests/smd_manifest.csv); \
	if [ -z "$$sys" ] || [ -z "$$var" ]; then echo "[indices-test] Need manifests/smd_manifest.csv with at least 1 row."; exit 2; fi; \
	echo "[indices-test] Running builder for $$sys $$var"; \
	PYTHONUNBUFFERED=1 timeout 300s stdbuf -oL -eL python3 scripts/build_indices_and_posres.py $$sys $$var

.PHONY: smd-submit-new
smd-submit-new:  ## Submit only NEW SMD rows (grouped by system+speed; array-capped; GPU mode)
	@bash scripts/smd_submit_new.sh manifests/smd_manifest.csv

.PHONY: smd-submit-new-cpu
smd-submit-new-cpu:  ## Submit only NEW SMD rows (CPU mode)
	@CPU=1 bash scripts/smd_submit_new.sh manifests/smd_manifest.csv

.PHONY: smd-clean-ledger
smd-clean-ledger:  ## Remove completed runs from smd_submitted.csv ledger
	@bash scripts/smd_clean_ledger.sh

.PHONY: preproc-traj
preproc-traj:  ## Preprocess all SMD trajectories (protein-only, centered, PBC-fixed, parallel). Override jobs with J=8
	@bash scripts/batch_preproc_traj.sh "" "" "$(J)"

.PHONY: preproc-traj-system
preproc-traj-system:  ## Preprocess trajectories for a system: make preproc-traj-system SYS=fimA_WT [VAR=AtoD] [J=8]
	@if [ -z "$(SYS)" ]; then echo "Usage: make preproc-traj-system SYS=<SYSTEM> [VAR=<VARIANT>] [J=<jobs>]"; exit 2; fi
	@bash scripts/batch_preproc_traj.sh "$(SYS)" "$(VAR)" "$(J)"

.PHONY: movies
movies:  ## Generate movies from analysis/preproc_traj into analysis/movies_out
	@bash scripts/make_movies_vmd.sh analysis/preproc_traj analysis/movies_out

.PHONY: movies-system
movies-system:  ## Generate movies for a specific system/variant; set SYS=<system> VAR=<variant>
	@if [ -z "$(SYS)" ]; then echo "Usage: make movies-system SYS=<SYSTEM> [VAR=<VARIANT>]"; exit 2; fi
	@ROOT=analysis/preproc_traj/$(SYS); \
	if [ -n "$(VAR)" ]; then ROOT=$$ROOT/$(VAR); fi; \
	bash scripts/make_movies_vmd.sh $$ROOT analysis/movies_out


