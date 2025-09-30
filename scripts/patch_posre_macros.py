#!/usr/bin/env python3
"""
Patch all posre*.itp files in CWD so their [ position_restraints ] entries
use the macro POSRE_FC for fcx/fcy/fcz. Idempotent.

Usage (run in systems/<SYS>/00_build):
    python3 patch_posre_macros.py
"""
import re
from pathlib import Path

def patch_one(p: Path) -> bool:
    text = p.read_text()
    out = []
    in_posre = False
    changed = False
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith('['):
            in_posre = (s.lower().startswith('[ position_restraints ]'))
            out.append(ln)
            continue
        if not in_posre or not s or s.startswith(';') or s.startswith('#'):
            out.append(ln)
            continue
        # Data line inside [ position_restraints ]
        cols = s.split()
        # Expect at least: i funct fcx fcy fcz
        try:
            _ = int(cols[0]); _ = int(cols[1])
        except Exception:
            out.append(ln); continue
        # Already macro-ized?
        if len(cols) >= 5 and (cols[2] == 'POSRE_FC' and cols[3] == 'POSRE_FC' and cols[4] == 'POSRE_FC'):
            out.append(ln); continue
        cols = [cols[0], cols[1], 'POSRE_FC', 'POSRE_FC', 'POSRE_FC']
        out.append((' '.join(cols)))
        changed = True
    if changed:
        p.write_text('\n'.join(out) + '\n')
    return changed

def main():
    cwd = Path('.')
    files = sorted([p for p in cwd.glob('posre*.itp') if p.is_file()])
    if not files:
        print("No posre*.itp files found here.")
        return
    total = 0
    for p in files:
        if patch_one(p):
            total += 1
            print(f"Patched {p}")
        else:
            print(f"Already OK: {p}")
    print(f"Done. Patched {total}/{len(files)} posre files.")

if __name__ == "__main__":
    main()
