#!/usr/bin/env python3
"""
Build per-system variant indices and a local anchor posres ITP,
and robustly rename selected groups to [ Anchor ] and [ Pulled ].

Usage:
    python3 scripts/build_indices_and_posres.py <SYSTEM> <VARIANT>

Outputs:
    indices/<SYSTEM>/<VARIANT>_pull.ndx
    indices/<SYSTEM>/<VARIANT>_posre_anchor_local.itp
"""

import os
import re
import sys
import yaml
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple

ROOT = Path(__file__).resolve().parents[1]
CFG  = ROOT / "config.yaml"

def sh(cmd, cwd=None, check=True, capture=False):
    if capture:
        return subprocess.run(cmd, cwd=cwd, check=check, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        return subprocess.run(cmd, cwd=cwd, check=check)

def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

# ---------- NDX helpers ----------

_and_tokens = re.compile(r'\b(and|&&|&)\b', re.I)
_res_tokens = re.compile(r'\b(residue|resid|resi|res|r)\b', re.I)
_chain_pat  = re.compile(r'\bchain[_\s]*([A-Za-z])\b', re.I)

def normalize_label(s: str) -> str:
    s = s.strip()
    s = s.replace('(', '').replace(')', '')
    s = _and_tokens.sub('_&_', s)
    s = _chain_pat.sub(lambda m: f"ch{m.group(1)}", s)
    s = _res_tokens.sub('r', s)
    s = re.sub(r'[\s]+', '_', s)
    s = re.sub(r'__+', '_', s)
    s = re.sub(r'\br_?([0-9])', r'r_\1', s)
    s = s.lower()
    return s

def load_ndx(path: Path) -> Tuple[List[str], List[List[int]]]:
    groups, atoms = [], []
    cur = None
    for line in path.read_text().splitlines():
        if line.strip().startswith('[') and line.strip().endswith(']'):
            name = line.strip()[1:-1].strip()
            groups.append(name)
            atoms.append([])
            cur = len(groups) - 1
        elif cur is not None:
            for tok in line.split():
                if tok.isdigit():
                    atoms[cur].append(int(tok))
    return groups, atoms

def write_ndx(path: Path, groups: List[str], atoms: List[List[int]]):
    out = []
    for name, lst in zip(groups, atoms):
        out.append(f"[ {name} ]")
        line = []
        for i, a in enumerate(lst, 1):
            line.append(str(a))
            if i % 15 == 0:
                out.append(" ".join(line))
                line = []
        if line:
            out.append(" ".join(line))
        out.append("")  # blank line between groups
    path.write_text("\n".join(out))

def choose_group(groups: List[str], pattern: str) -> int:
    norm_pat = normalize_label(pattern)
    candidates = []
    for i, g in enumerate(groups):
        ng = normalize_label(g)
        if ng == norm_pat:
            candidates.append((i, g, 'exact'))
        elif (norm_pat in ng) or (ng in norm_pat):
            candidates.append((i, g, 'fuzzy'))

    if not candidates:
        # relaxed stripped compare
        strip = lambda x: re.sub(r'[^a-z0-9]+', '', normalize_label(x))
        sp = strip(pattern)
        for i, g in enumerate(groups):
            if strip(g) == sp:
                candidates.append((i, g, 'stripped'))

    if not candidates:
        die(f"No group matched pattern: '{pattern}' (normalized: '{norm_pat}')")

    exact = [c for c in candidates if c[2] == 'exact']
    if len(exact) == 1:
        return exact[0][0]
    if len(candidates) == 1:
        return candidates[0][0]

    msg = ["Pattern is ambiguous. Candidates:"]
    for i, g, kind in candidates:
        msg.append(f"  - [{i}] '{g}' (match={kind}, norm='{normalize_label(g)}')")
    die("\n".join(msg))

# ---------- Config helpers ----------

def load_patterns(cfg, system: str, variant: str):
    # Variant-level first, then system-level fallback
    sys_cfg = cfg.get("systems", {}).get(system, {})
    var_cfg = (sys_cfg.get("variants") or {}).get(variant, {})

    anchor = var_cfg.get("anchor_pattern") or sys_cfg.get("anchor_pattern")
    pulled = var_cfg.get("pulled_pattern") or sys_cfg.get("pulled_pattern")
    if not anchor or not pulled:
        die(f"Missing anchor/pulled patterns for {system}/{variant}. "
            f"Provide anchor_pattern & pulled_pattern in config.yaml (variant or system level).")
    return anchor, pulled

# ---------- Posres anchor (local) ----------

def write_local_posre_itp(outfile: Path, atom_indices: List[int], fc_kj_mol_nm2: float = 1000.0):
    """
    Write a minimal posre file for the anchor atoms in *local* coordinates.
    """
    lines = [
        "; Local posre for Anchor group (generated)",
        "[ position_restraints ]",
        "; ai  funct      fc_x      fc_y      fc_z"
    ]
    for ai in atom_indices:
        lines.append(f"{ai:6d}     1   {fc_kj_mol_nm2:.1f}   {fc_kj_mol_nm2:.1f}   {fc_kj_mol_nm2:.1f}")
    outfile.write_text("\n".join(lines) + "\n")

# ---------- Main ----------

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <SYSTEM> <VARIANT>", file=sys.stderr)
        sys.exit(2)

    system  = sys.argv[1]
    variant = sys.argv[2]

    cfg = yaml.safe_load(CFG.read_text())
    anchor_pat, pulled_pat = load_patterns(cfg, system, variant)

    build_dir = ROOT / "systems" / system / "00_build"
    if not (build_dir / "npt_final.gro").exists():
        die(f"Build artifact missing: {build_dir/'npt_final.gro'}")

    indices_dir = ROOT / "indices" / system
    indices_dir.mkdir(parents=True, exist_ok=True)

    tmp_ndx   = indices_dir / f"{variant}_tmp.ndx"
    pull_ndx  = indices_dir / f"{variant}_pull.ndx"
    posre_itp = indices_dir / f"{variant}_posre_anchor_local.itp"

    # 1) Start from a fresh index (Protein / non-Protein, etc.)
    # We rely on gromacs default groups from the TPR.
    tpr = build_dir / "npt_final.tpr"
    if not tpr.exists():
        # create tpr quickly from gro/top if needed (rare)
        die(f"Missing {tpr} (expected from final NPT).")

    # Create a base ndx from the tpr (includes standard groups)
    # gmx make_ndx -f <tpr> -o tmp.ndx <<< "q"
    sh(["gmx", "make_ndx", "-f", str(tpr), "-o", str(tmp_ndx)], check=True, capture=True).stdout

    # 2) Load and rename to Anchor/Pulled
    groups, atoms = load_ndx(tmp_ndx)

    # Create NonProtein if missing (helper for templates)
    try:
        idx_np = groups.index("Non-Protein")
    except ValueError:
        # Build "non-protein" as all atoms not in Protein
        try:
            idx_prot = groups.index("Protein")
        except ValueError:
            idx_prot = None
        if idx_prot is not None:
            all_atoms = set(a for lst in atoms for a in lst)
            prot_atoms = set(atoms[idx_prot])
            nonprot = sorted(all_atoms - prot_atoms)
            groups.append("NonProtein")
            atoms.append(nonprot)

    # If Anchor/Pulled exist already, we’ll replace names anyway to ensure consistency.
    # Choose matches using robust normalization.
    idx_anchor = choose_group(groups, anchor_pat)
    idx_pulled = choose_group(groups, pulled_pat)

    old_anchor, old_pulled = groups[idx_anchor], groups[idx_pulled]
    groups[idx_anchor] = "Anchor"
    groups[idx_pulled] = "Pulled"

    write_ndx(tmp_ndx, groups, atoms)

    # 3) Persist final .ndx under the expected name
    shutil.copyfile(tmp_ndx, pull_ndx)
    print(f"[OK] Renamed groups: '{old_anchor}'→'Anchor', '{old_pulled}'→'Pulled'")
    print(f"[OK] Wrote {pull_ndx}")

    # 4) Write a local posre for Anchor (force constant from config if present)
    fc = (cfg.get("globals", {}) or {}).get("anchor_posre_k_kj_mol_nm2", 1000.0)
    write_local_posre_itp(posre_itp, atoms[idx_anchor], float(fc))
    print(f"[OK] Wrote {posre_itp} (fc={fc} kJ/mol/nm^2)")

if __name__ == "__main__":
    main()
