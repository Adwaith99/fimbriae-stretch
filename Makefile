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

.PHONY: posteq-sample smd-manifest smd-submit smd-qc

# 1) Sample start frames from final NPT (per system/variant), random picks per variant
posteq-sample:
	@bash pipelines/posteq_sample_starts.sh

# 2) Build SMD manifest (systems × variants × speeds × replicates; each replicate gets one start)
smd-manifest:
	@python3 scripts/smd_build_manifest.py

# 3) Submit one Slurm array per system (capped by globals.slurm.array_cap), no requeue
smd-submit: smd-manifest
	@bash pipelines/smd_array_submit.sh

# 4) Optional QC pass over completed runs (TSV and PNG if enabled)
smd-qc:
	@python3 scripts/smd_qc.py

