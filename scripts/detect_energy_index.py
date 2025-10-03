#!/usr/bin/env python3
import sys, subprocess, re

if len(sys.argv) < 3:
    print("usage: detect_energy_index.py <file.edr> <Label>", file=sys.stderr)
    sys.exit(2)

edr = sys.argv[1]
want = sys.argv[2].lower()

try:
    out = subprocess.check_output(["gmx", "dump", "-e", edr], text=True, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
    print(e.output, file=sys.stderr)
    sys.exit(3)

pairs = []
for line in out.splitlines():
    if "energies:" in line:
        for m in re.finditer(r"(\d+)\s+'([^']+)'", line):
            idx = int(m.group(1))
            lab = m.group(2)
            pairs.append((idx, lab))

# Prefer exact (case-insensitive), else substring
best = None
for idx, lab in pairs:
    if lab.lower() == want:
        best = idx; break
if best is None:
    for idx, lab in pairs:
        if want in lab.lower():
            best = idx; break

if best is None:
    print(f"ERROR: Could not find energy term matching '{want}' in {edr}", file=sys.stderr)
    sys.exit(4)

# gmx energy’s interactive menu is 1-based
print(best + 1)
