#!/usr/bin/env bash
set -euo pipefail

# Minimal test for find_pbcatom.py
# Requires an existing systems/*/00_build directory with start.gro and index.ndx

echo "[test-find-pbcatom] Running basic sanity checks..."

# Check if the script exists
if [[ ! -f "scripts/find_pbcatom.py" ]]; then
  echo "[test-find-pbcatom] ERROR: scripts/find_pbcatom.py not found" >&2
  exit 1
fi

# Try to parse help
python3 scripts/find_pbcatom.py --help > /dev/null 2>&1 && {
  echo "[test-find-pbcatom] ✓ find_pbcatom.py --help works"
} || {
  echo "[test-find-pbcatom] ERROR: find_pbcatom.py --help failed" >&2
  exit 1
}

# Try running with invalid arguments (should error gracefully)
python3 scripts/find_pbcatom.py --gro nonexistent.gro --ndx nonexistent.ndx --group Test --resid 1 --atom CA 2>&1 | grep -q "ERROR" && {
  echo "[test-find-pbcatom] ✓ find_pbcatom.py errors gracefully on missing files"
} || {
  echo "[test-find-pbcatom] WARNING: Expected error message not found" >&2
}

echo "[test-find-pbcatom] Basic checks passed."
