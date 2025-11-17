#!/usr/bin/env bash
set -euo pipefail

# Report completion status for all runs under smd/
# Outputs TSV: system\tvariant\tspeed\trep\tstart\tcompleted\tlast_step\texpected_nsteps\tremaining_steps\tns_remaining\test_hours

find_root(){
  local p="$PWD"
  while [[ "$p" != "/" ]]; do
    [[ -f "$p/config.yaml" ]] && { printf "%s" "$p"; return 0; }
    p=$(dirname "$p")
  done
  return 1
}

ROOT=$(find_root)
if [[ -z "${ROOT:-}" ]]; then
  echo "ERROR: cannot find repo root (config.yaml)" >&2
  exit 2
fi

SMD_DIR="$ROOT/smd"
if [[ ! -d "$SMD_DIR" ]]; then
  echo "No smd/ directory found under repo root ($ROOT)." >&2
  exit 0
fi

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "system" "variant" "speed" "rep" "start" "completed" "last_step" "expected_nsteps" "remaining_steps" "ns_remaining" "est_hours" > "$TMP"

# load system perf map and globals (safety)
export ROOT
readarray -t PERF_MAP < <(python3 - <<PY
import yaml, os
cfg=yaml.safe_load(open(os.path.join(os.environ['ROOT'],'config.yaml')))
global_perf=float(cfg.get('globals',{}).get('perf_ns_per_day', 10.0))
for s in cfg.get('systems',[]):
  name=s.get('name')
  perf=float(s.get('perf_ns_per_day', global_perf))
  print(f"{name}:{perf}")
print(f"__GLOBAL__:{global_perf}")
print(f"__SAFETY__:{cfg.get('globals',{}).get('walltime_safety',1.20)}")
PY
)
declare -A system_perf
GLOBAL_PERF=10.0
SAFETY=1.2
for entry in "${PERF_MAP[@]}"; do
  k=${entry%%:*}
  v=${entry#*:}
  if [[ "$k" == "__GLOBAL__" ]]; then
    GLOBAL_PERF="$v"
  elif [[ "$k" == "__SAFETY__" ]]; then
    SAFETY="$v"
  else
    system_perf["$k"]="$v"
  fi
done

# iterate runs via expected_nsteps.txt presence
while IFS= read -r -d $'\0' file; do
  dir=$(dirname "$file")
  rel=${dir#${SMD_DIR}/}
  IFS='/' read -r system variant vs rep startid <<< "$rel"
  speed="${vs#v}"

  expected=$(cat "$file" 2>/dev/null || echo "")
  last_step=""
  if [[ -f "$dir/pull.log" ]]; then
    export DIR="$dir"
    last_step=$(python3 - <<PY
import re,os
f=os.path.join(os.environ['DIR'],'pull.log')
max_step=0
try:
  with open(f) as fh:
    for line in fh:
      m=re.search(r'(?:^Step\s+|Statistics over\s+)(\d+)', line)
      if m:
        v=int(m.group(1))
        if v>max_step:
          max_step=v
except Exception:
  pass
print(max_step if max_step>0 else "")
PY
  fi

  # determine dt_ps (from pull.mdp or config global dt_ps)
  dt_ps=""
  if [[ -f "$dir/pull.mdp" ]]; then
    dt_ps=$(awk 'tolower($1)=="dt"{print $3}' "$dir/pull.mdp" | tail -n1)
  fi
  if [[ -z "$dt_ps" ]]; then
    dt_ps=$(python3 - <<PY
import yaml, os
cfg=yaml.safe_load(open(os.path.join(os.environ['ROOT'],'config.yaml')))
print(cfg.get('globals',{}).get('dt_ps',0.002))
PY
)

  fi

  # numericize
  expected_n=0
  remaining_steps=0
  if [[ "$expected" =~ ^[0-9]+$ ]]; then
    expected_n=$expected
  fi
  if [[ "$last_step" =~ ^[0-9]+$ ]]; then
    remaining_steps=$(( expected_n - last_step ))
    if (( remaining_steps < 0 )); then remaining_steps=0; fi
  else
    remaining_steps=$expected_n
  fi

  ns_remaining="0"
  if [[ "$remaining_steps" =~ ^[0-9]+$ ]]; then
    # dt_ps is in ps -> ns = (steps * dt_ps) / 1000
    ns_remaining=$(awk -v steps="$remaining_steps" -v dt="$dt_ps" 'BEGIN{ if (dt==""||dt==0) {print "0"} else {ns=(steps*dt)/1000.0; printf "%.6f", ns}}')
  fi

  # estimate hours using system perf
  perf=${system_perf[$system]:-$GLOBAL_PERF}
  est_hours=""
  if [[ -n "$ns_remaining" && $(echo "$ns_remaining > 0" | bc -l) -eq 1 ]]; then
    est_hours=$(awk -v ns="$ns_remaining" -v perf="$perf" -v pad="$SAFETY" 'BEGIN{ if (ns==""||ns==0||perf==0) print 0; else {hrs=(ns/perf)*24.0*pad; printf "%.3f", hrs}}')
  else
    est_hours=0
  fi

  completed=0
  if [[ "$last_step" =~ ^[0-9]+$ ]] && [[ "$expected_n" =~ ^[0-9]+$ ]] && (( last_step >= expected_n )); then
    completed=1
    remaining_steps=0
    ns_remaining=0
    est_hours=0
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$system" "$variant" "$speed" "${rep#rep}" "$startid" "$completed" "$last_step" "$expected_n" "$remaining_steps" "$ns_remaining" "$est_hours" >> "$TMP"

done < <(find "$SMD_DIR" -type f -name expected_nsteps.txt -print0)

# Summarize and print human-friendly report
total=$(($(wc -l < "$TMP") - 1))
completed=$(awk -F"\t" 'NR>1{c+=$6}END{print c+0}' "$TMP")
incomplete=$(( total - completed ))
sum_hours=$(awk -F"\t" 'NR>1{h+=($11+0)}END{printf "%.3f", h+0}' "$TMP")

printf "\nReport for smd/ (found %d runs)\n" "$total"
printf "Completed: %d\n" "$completed"
printf "Incomplete: %d\n" "$incomplete"
printf "Estimated wall-hours remaining (sum over incomplete): %s h\n\n" "$sum_hours"

# Print table: try column for pretty output, else raw TSV
if command -v column >/dev/null 2>&1; then
  column -t -s$'\t' "$TMP"
else
  cat "$TMP"
fi

exit 0
