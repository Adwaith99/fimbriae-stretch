#!/usr/bin/env bash
set -euo pipefail

# Submit pulling jobs as ONE ARRAY PER SYSTEM instead of a monolithic array.
# Benefits:
#  * Allows per-system walltime selection (heterogeneous performance / size)
#  * Can cap concurrency per system using a global array_cap
#  * Easier cancellation / monitoring per system (one job id each)

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
MAN="$ROOT/manifests/manifest.csv"
CFG="$ROOT/config.yaml"

[[ -f "$MAN" ]] || { echo "Manifest not found: $MAN"; exit 2; }

# Read global array cap (fallback 5) and an optional pulling performance map.
ARRAY_CAP=$(python3 - "$CFG" <<'PY'
import yaml, sys
try:
  c=yaml.safe_load(open(sys.argv[1]))
  cap=c.get('globals',{}).get('slurm',{}).get('array_cap',5)
  print(int(cap))
except Exception:
  print(5)
PY
)

# Optional config extension (not yet documented): per-system expected ns_per_day for pulling (to derive walltime).
# If absent we skip dynamic walltime override and use pulls_array_submit.sh defaults.

echo "DEBUG: Filtering manifest rows for status=='pending' (case-insensitive) or empty."
readarray -t PENDING < <(awk -F, 'BEGIN{IGNORECASE=1} NR==1{hdr=$0;next} {
  for(i=1;i<=NF;i++){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", $i) }
  print "DEBUG: Row status=", $7 > "/dev/stderr"
  status=$7; if(status=="" || tolower(status)=="pending") print $0 }' "$MAN")

if (( ${#PENDING[@]} == 0 )); then
  echo "No pending pull rows in $MAN"; exit 0
fi

echo "Total pending rows: ${#PENDING[@]}"
echo "DEBUG: Pending rows:"
for row in "${PENDING[@]}"; do
  echo "$row"
done

# Group rows per system -> write one filtered manifest per system -> submit separate arrays.
TMPDIR="$ROOT/manifests/_pending_split"; mkdir -p "$TMPDIR"
HEADER=$(head -n1 "$MAN")
declare -A SYS_ROWS
for line in "${PENDING[@]}"; do
  # Split space-separated fields: uid system variant speed_nm_per_ns k_kj_mol_nm2 rep status
  IFS=' ' read -r uid system variant speed_nm_per_ns k_kj_mol_nm2 rep status <<< "$line"
  sys="$system"
  [[ -z "$sys" ]] && continue
  SYS_ROWS["$sys"]+="$line"$'\n'
done

for sys in "${!SYS_ROWS[@]}"; do
  rows=${SYS_ROWS[$sys]}
  [[ -z "$rows" ]] && continue
  RUNSTAMP=$(date +'%Y%m%d-%H%M%S')
  RUNMAN="$TMPDIR/manifest.${sys}.pending.$RUNSTAMP.csv"
  { echo "$HEADER"; printf '%s' "$rows"; } > "$RUNMAN"
  n_rows=$(awk 'NR>1 && NF>0{c++} END{print c+0}' "$RUNMAN")
  if (( n_rows == 0 )); then
    continue
  fi

  # Estimate walltime (placeholder: could be refined). We compute time to target extension
  # for the *slowest* speed in this system subset to guarantee coverage.
  WALL=$(python3 - "$CFG" "$RUNMAN" "$sys" <<'PY'
import yaml, csv, sys, math
cfg=yaml.safe_load(open(sys.argv[1]))
g=cfg.get('globals',{})
target=g.get('target_extension_nm',5.0)
global_perf=g.get('perf_ns_per_day')  # may be None if user removes it
safety=g.get('walltime_safety',1.25)
max_h=g.get('max_wall_hours',168)
systems=cfg.get('systems',[])
sys_name=sys.argv[3]
# find system entry
sys_entry=next((s for s in systems if s.get('name')==sys_name), {})
sys_perf=None
# Allow either `perf_ns_per_day` at top level of system or nested pulling: { perf_ns_per_day: X }
if 'perf_ns_per_day' in sys_entry:
    sys_perf=sys_entry['perf_ns_per_day']
elif isinstance(sys_entry.get('pulling'), dict) and 'perf_ns_per_day' in sys_entry['pulling']:
    sys_perf=sys_entry['pulling']['perf_ns_per_day']
perf = sys_perf if sys_perf is not None else (global_perf if global_perf is not None else 10.0)
system_rows=list(csv.DictReader(open(sys.argv[2])))
if not system_rows:
    print('24:00:00'); sys.exit(0)
speeds=sorted({float(r['speed_nm_per_ns']) for r in system_rows}) or [0.1]
slowest=min(speeds)
needed_ns=target/slowest
days=(needed_ns/perf)*safety
hours_total=max(1, days*24.0)
hours_ceil=min(int(math.ceil(hours_total)), int(max_h))
print(f"{hours_ceil:02d}:00:00")
PY
)

  # Cap concurrency per system
  CAP_USED=$(( ARRAY_CAP>0 ? ARRAY_CAP : 1 ))
  echo "Submitting system=${sys} rows=${n_rows} walltime=${WALL} cap=${CAP_USED} manifest=$(basename "$RUNMAN")"
  sbatch --export=ALL,FIM_PULL_MANIFEST="$RUNMAN" \
         --job-name="pulls_${sys}" \
         --time="$WALL" \
         --array=1-${n_rows}%${CAP_USED} \
         pipelines/pulls_array_submit.sh
done

echo "Submission complete (one array per system)."
