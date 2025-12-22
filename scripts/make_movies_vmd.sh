#!/usr/bin/env bash
set -euo pipefail

# Generate MP4 movies from preprocessed protein trajectories using VMD+ffmpeg.
# Defaults are aligned with this repo:
# - Input root: analysis/preproc_traj
# - Output dir: analysis/movies_out
#
# Usage:
#   scripts/make_movies_vmd.sh [ROOT] [OUTDIR]
#
# Env overrides:
#   STRIDE (default 1), FPS (default 30), WIDTH (1600), HEIGHT (900)
#   VMD_BIN (default vmd), XVFB_RUN (default xvfb-run -a if available else empty)
#   RANGES_FILE (optional path; if set and exists, passed to TCL script)
#
# Requires: vmd, ffmpeg; xvfb-run recommended on headless nodes.

ROOT="${1:-analysis/preproc_traj}"
OUTDIR="${2:-analysis/movies_out}"

STRIDE="${STRIDE:-1}"
FPS="${FPS:-30}"
WIDTH="${WIDTH:-1600}"
HEIGHT="${HEIGHT:-900}"
VMD_BIN="${VMD_BIN:-vmd}"
RANGES_FILE="${RANGES_FILE:-}"

XVFB_RUN="${XVFB_RUN:-}"
if [[ -z "${XVFB_RUN}" ]]; then
  if command -v xvfb-run >/dev/null 2>&1; then
    XVFB_RUN="xvfb-run -a"
  else
    XVFB_RUN=""
  fi
fi

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/_frames"

# Find leaf run dirs that contain start_protein.gro
mapfile -t run_dirs < <(find "$ROOT" -type f -name start_protein.gro -print0 \
  | xargs -0 -n1 dirname | sort -u)

if [[ ${#run_dirs[@]} -eq 0 ]]; then
  echo "No runs found under $ROOT (looking for start_protein.gro)"
  exit 0
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TCL_SCRIPT="$script_dir/vmd_make_movie_ppm_chains.tcl"
if [[ ! -f "$TCL_SCRIPT" ]]; then
  echo "ERROR: Missing TCL script: $TCL_SCRIPT" >&2
  exit 1
fi

for d in "${run_dirs[@]}"; do
  gro="$d/start_protein.gro"

  shopt -s nullglob
  xtcs=( "$d"/traj_protein_dt*.xtc )
  shopt -u nullglob

  if [[ ${#xtcs[@]} -eq 0 ]]; then
    echo "Skip (no xtc): $d"
    continue
  fi

  # Build name from hierarchy: system_variant_speed_rep_sXXX + traj base
  rel="${d#"$ROOT"/}"                 # e.g. fimA_WT/AtoD/v0.020/rep1/s035
  rel_slug="$(echo "$rel" | tr '/ ' '__')"   # safe filename-ish

  for xtc in "${xtcs[@]}"; do
    base="$(basename "$xtc" .xtc)"           # e.g. traj_protein_dt1000
    slug="${rel_slug}__${base}"

    framesdir="$OUTDIR/_frames/$slug"
    mkdir -p "$framesdir"

    outpre="$framesdir/$slug"               # VMD will make: outpre_00000.ppm + outpre.mp4
    final_mp4="$OUTDIR/${slug}.mp4"

    # Skip if movie already exists and looks non-empty
    if [[ -f "$final_mp4" && -s "$final_mp4" ]]; then
      echo "Skip (exists): $(basename "$final_mp4")"
      continue
    fi

    echo "==> $rel :: $base -> $(basename "$final_mp4")"

    # Optional ranges file
    ranges_arg=""
    if [[ -n "$RANGES_FILE" && -f "$RANGES_FILE" ]]; then
      ranges_arg="$RANGES_FILE"
    fi

    # Run VMD (headless if xvfb-run is available)
    set +e
    ${XVFB_RUN} "$VMD_BIN" -dispdev opengl -e "$TCL_SCRIPT" \
      -args "$gro" "$xtc" "$outpre" "$STRIDE" "$FPS" "$WIDTH" "$HEIGHT" "$ranges_arg"
    status=$?
    set -e
    if [[ $status -ne 0 ]]; then
      echo "WARNING: VMD/ffmpeg failed for $xtc (status $status)"
      # leave framesdir for inspection
      continue
    fi

    # Move the final mp4 into OUTDIR, then clean frames
    if [[ -f "${outpre}.mp4" ]]; then
      mv -f "${outpre}.mp4" "$final_mp4"
      rm -rf "$framesdir"
    else
      echo "WARNING: missing mp4 for $d ($xtc)"
    fi
  done
done

echo "DONE. Movies are in: $OUTDIR"