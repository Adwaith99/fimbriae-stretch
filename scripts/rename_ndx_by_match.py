#!/usr/bin/env python3
"""
Rename a group in a GROMACS index (.ndx) file by *pattern-matching* its label.

Example:
    python3 rename_ndx_by_match.py input.ndx "chain A & r 135-150" Pulled

This script normalizes both the search pattern and existing group labels so that
things like "chain A & r 135-150" match labels such as "chA_&_r_135-150".
"""

import re
import sys
from pathlib import Path
from typing import List, Tuple

def die(msg: str, code: int = 1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def load_ndx(path: Path) -> Tuple[List[str], List[List[int]]]:
    groups, atoms = [], []
    cur = None
    for line in path.read_text().splitlines():
        if line.strip().startswith('[') and line.strip().endswith(']'):
            name = line.strip()[1:-1].strip()
            groups.append(name)
            atoms.append([])
            cur = len(groups) - 1
        else:
            if cur is not None:
                for token in line.split():
                    if token.isdigit():
                        atoms[cur].append(int(token))
    return groups, atoms

# --- Normalization helpers ----------------------------------------------------

_and_tokens = re.compile(r'\b(and|&&|&)\b', re.I)
_res_tokens = re.compile(r'\b(residue|resid|resi|res)\b', re.I)
_chain_pat  = re.compile(r'\bchain[_\s]*([A-Za-z])\b', re.I)

def normalize_label(s: str) -> str:
    """
    Make 'human' patterns and .ndx labels comparable.

    Examples:
        "chain A & r 135-150"  -> "cha_&_r_135-150"
        "chA_&_r_135-150"      -> "cha_&_r_135-150"
        "CHAIN B and RES 5"    -> "chb_&_r_5"
    """
    s = s.strip()
    s = s.replace('(', '').replace(')', '')
    # unify logical AND to &_ (keep explicit marker)
    s = _and_tokens.sub('_&_', s)
    # chain X -> chX (lowercase 'ch' + letter)
    s = _chain_pat.sub(lambda m: f"ch{m.group(1)}", s)
    # residue/resid/res/resi -> r_
    s = _res_tokens.sub('r', s)
    # unify spaces/hyphens
    s = re.sub(r'[\s]+', '_', s)
    # unify multiple underscores
    s = re.sub(r'__+', '_', s)
    # allow "r135" vs "r_135"
    s = re.sub(r'\br_?([0-9])', r'r_\1', s)
    # lowercase everything
    s = s.lower()
    return s

def choose_group(groups: List[str], pattern: str) -> int:
    norm_pat = normalize_label(pattern)
    candidates = []
    for i, g in enumerate(groups):
        ng = normalize_label(g)
        if ng == norm_pat:
            candidates.append((i, g, 'exact'))
        elif ng.find(norm_pat) != -1 or norm_pat.find(ng) != -1:
            candidates.append((i, g, 'fuzzy'))

    if not candidates:
        # as a last resort, try removing non-alphanumerics and compare
        strip = lambda x: re.sub(r'[^a-z0-9]+', '', normalize_label(x))
        sp = strip(pattern)
        last = []
        for i, g in enumerate(groups):
            if strip(g) == sp:
                last.append((i, g, 'stripped'))
        candidates = last

    if not candidates:
        die(f"No group matched pattern: '{pattern}' (normalized: '{norm_pat}')")

    # prefer exact over fuzzy/stripped
    exact = [c for c in candidates if c[2] == 'exact']
    if len(exact) == 1:
        return exact[0][0]
    if len(candidates) == 1:
        return candidates[0][0]

    # ambiguous
    msg = ["Pattern is ambiguous. Candidates:"]
    for i, g, kind in candidates:
        msg.append(f"  - [{i}] '{g}' (match={kind}, norm='{normalize_label(g)}')")
    die("\n".join(msg))

def write_ndx(path: Path, groups: List[str], atoms: List[List[int]]):
    out_lines = []
    for name, lst in zip(groups, atoms):
        out_lines.append(f"[ {name} ]")
        # wrap 15 integers per line for readability
        line = []
        for i, a in enumerate(lst, 1):
            line.append(str(a))
            if i % 15 == 0:
                out_lines.append(" ".join(line))
                line = []
        if line:
            out_lines.append(" ".join(line))
        out_lines.append("")  # blank line between groups
    path.write_text("\n".join(out_lines))

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.ndx> <pattern to match> <new name>", file=sys.stderr)
        sys.exit(2)

    ndx_in = Path(sys.argv[1])
    pattern = sys.argv[2]
    new_name = sys.argv[3]

    if not ndx_in.exists():
        die(f"Index file not found: {ndx_in}")

    groups, atoms = load_ndx(ndx_in)
    idx = choose_group(groups, pattern)

    old = groups[idx]
    groups[idx] = new_name
    write_ndx(ndx_in, groups, atoms)

    print(f"Renamed group: '{old}'  →  '{new_name}' in {ndx_in}")

if __name__ == "__main__":
    main()
