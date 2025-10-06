#!/usr/bin/env python3
import sys, subprocess, re, os

if len(sys.argv) < 3:
    print("usage: detect_energy_index.py <file.edr> <Label>", file=sys.stderr)
    sys.exit(2)

edr = sys.argv[1]
want = sys.argv[2].lower().strip()

# Accept common aliases
ALIASES = {
    "temperature": ["temperature", "temp.", "temp"],
    "pressure":    ["pressure", "pres.", "press"],
    "density":     ["density", "mass density", "rho"],
    "potential":   ["potential", "potential energy"],
}
candidates = ALIASES.get(want, [want])

# quick sanity: non-empty file
try:
    if os.path.getsize(edr) == 0:
        print(f"ERROR: {edr} is zero-length", file=sys.stderr)
        sys.exit(5)
except OSError as e:
    print(f"ERROR: cannot stat {edr}: {e}", file=sys.stderr)
    sys.exit(5)

# dump the energy terms
try:
    out = subprocess.check_output(["gmx", "dump", "-e", edr], text=True, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
    print(e.output, file=sys.stderr)
    sys.exit(3)

pairs = []
# collect any N 'Label' pairs across the whole dump (not just a single line)
for m in re.finditer(r"(\d+)\s+'([^']+)'", out):
    idx = int(m.group(1))
    lab = m.group(2).strip()
    pairs.append((idx, lab))

# prefer exact match among aliases, then substring
low = [(i, s, s.lower()) for i, s in pairs]
best = None
for _, _, lablow in low:
    if any(lablow == a for a in candidates):
        best = _; break
if best is None:
    for i, _, lablow in low:
        if any(a in lablow for a in candidates):
            best = i; break

if best is None:
    print(f"ERROR: Could not find energy term matching {candidates} in {edr}", file=sys.stderr)
    sys.exit(4)

print(best + 1)  # gmx energy menu is 1-based
