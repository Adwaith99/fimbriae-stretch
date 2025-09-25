SHELL := /bin/bash
ROOT  := $(shell pwd)

help:
	@echo "Targets (run on the cluster login node):"
	@echo "  pull             - git pull (rebase) latest code"
	@echo "  validate         - sanity check config.yaml"
	@echo "  manifest         - (re)generate manifests"
	@echo "  build            - submit system build jobs (only if needed)"
	@echo "  pulls            - submit SMD array over PENDING rows"
	@echo "  status           - show your Slurm queue"
	@echo "  tail-build SYS=...      - tail build log for a system"
	@echo "  tail-pull  A=JOBID a=TASKID - tail a pull array task log"
	@echo "  cancel A=JOBID    - cancel a Slurm job/array"

pull:
	git pull --rebase

validate:
	python3 scripts/validate_config.py

manifest:
	python3 scripts/generate_manifest.py

build: manifest
	bash scripts/submit_builds.sh

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
