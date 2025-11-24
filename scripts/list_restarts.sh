#!/usr/bin/env bash
set -euo pipefail

# Find repo root (look for config.yaml)
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

# seen dirs map
declare -A _seen_dirs
printf "system\tvariant\tspeed\trep\tstart\trestart_count\trestarts\tlatest_restart\n"

# Use find to iterate unique dirs containing restart markers
while IFS= read -r -d $'\0' marker; do
  dir=$(dirname "$marker")
  # Skip if we've already processed this directory (we iterate over markers, so dedupe)
  if [[ -n "${_seen_dirs[$dir]:-}" ]]; then
    continue
  fi
  _seen_dirs[$dir]=1

  # Determine run components from path: smd/<system>/<variant>/v<speed>/rep<rep>/<start>
  rel=${dir#${SMD_DIR}/}
  IFS='/' read -r system variant vs rep startid <<< "$rel"
  speed="${vs#v}"

  # Collect restart markers in this dir
  markers=$(ls -1 "$dir"/.restart_* 2>/dev/null || true)
  if [[ -z "$markers" ]]; then
    continue
  fi
  # basename and sort
  markers_list=$(printf "%s\n" $markers | xargs -n1 basename | sort -V | tr '\n' ',' | sed 's/,$//')
  restart_count=$(echo "$markers_list" | awk -F, '{ if ($0=="") print 0; else print NF }')
  # latest marker time
  latest_file=$(ls -1t "$dir"/.restart_* 2>/dev/null | head -n1 || true)
  latest_ts=""
  if [[ -n "$latest_file" ]]; then
    latest_ts=$(stat -c '%y' "$latest_file" 2>/dev/null || echo "")
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$system" "$variant" "$speed" "${rep#rep}" "$startid" "$restart_count" "$markers_list" "$latest_ts"

done < <(find "$SMD_DIR" -type f -name ".restart_*" -print0)

exit 0
