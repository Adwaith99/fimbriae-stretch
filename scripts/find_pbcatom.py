#!/usr/bin/env python3
"""
Find a specific atom in a .gro file by (resid, atomname) within a given index group.

Usage:
  python3 find_pbcatom.py --gro start.gro --ndx index.ndx --group Anchor --resid 218 --atom CA

Output:
  prints the 1-based atom index (as in GROMACS/GRO) to stdout
  or exits with error if not found.
"""

import argparse
import sys
import subprocess
import re


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gro", required=True, help="Path to .gro file")
    p.add_argument("--ndx", required=True, help="Path to index.ndx file")
    p.add_argument("--group", required=True, help="Index group name (e.g., Anchor, Pulled)")
    p.add_argument("--resid", type=int, required=True, help="Residue number (1-based)")
    p.add_argument("--atom", default="CA", help="Atom name (default: CA)")
    return p.parse_args()


def parse_gro(gro_file):
    """
    Parse .gro file and return list of dicts with keys:
      - index (1-based, as in .gro)
      - resid (int)
      - resname (str)
      - atomname (str)
      - x, y, z (float)
    """
    atoms = []
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        raise ValueError("Invalid .gro file: too few lines")
    
    # Line 0: title, Line 1: natoms
    try:
        natoms = int(lines[1].strip())
    except ValueError:
        raise ValueError("Invalid .gro file: cannot parse natoms from line 2")
    
    # Lines 2 onwards: atom records
    for i in range(natoms):
        if 2 + i >= len(lines):
            raise ValueError(f"Invalid .gro file: missing atom line {i+1}")
        
        line = lines[2 + i]
        # .gro format (fixed-width):
        # residue number (5 pos, int)
        # residue name (5 char, str)
        # atom name (5 char, str)
        # atom number (5 pos, int)
        # position (in nm, x y z in 3 columns, each 8 pos with 4 decimal places)
        # velocity (in km/s (or m/s), x y z in 3 columns, each 8 pos with 4 decimal places)
        
        try:
            resid = int(line[0:5].strip())
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atom_index_gro = int(line[15:20].strip())
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            
            atoms.append({
                'index': atom_index_gro,
                'resid': resid,
                'resname': resname,
                'atomname': atomname,
                'x': x,
                'y': y,
                'z': z
            })
        except (ValueError, IndexError) as e:
            raise ValueError(f"Cannot parse atom line {i+1}: {line.rstrip()}\n{e}")
    
    return atoms


def get_group_atoms(ndx_file, group_name):
    """
    Extract atom indices from index.ndx for a given group.
    Returns a set of 1-based atom indices.
    """
    atoms = set()
    in_group = False
    with open(ndx_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            # Check if this is the group header
            if re.match(rf'^\[\s*{re.escape(group_name)}\s*\]', line):
                in_group = True
                continue
            
            # If we see another group header, stop
            if in_group and re.match(r'^\[\s*\w+\s*\]', line):
                break
            
            # Parse atom indices from this group
            if in_group and line.strip():
                # Split by whitespace and parse integers
                tokens = line.split()
                for tok in tokens:
                    try:
                        atoms.add(int(tok))
                    except ValueError:
                        pass
    
    if not in_group:
        raise ValueError(f"Group '{group_name}' not found in {ndx_file}")
    
    return atoms


def find_pbcatom(gro_file, ndx_file, group_name, target_resid, target_atomname):
    """
    Find the atom index for (group, resid, atomname).
    Returns 1-based atom index.
    """
    atoms = parse_gro(gro_file)
    group_atoms = get_group_atoms(ndx_file, group_name)
    
    # Filter atoms in the group with matching resid and atomname
    matches = [
        a for a in atoms
        if a['index'] in group_atoms
        and a['resid'] == target_resid
        and a['atomname'] == target_atomname
    ]
    
    if len(matches) == 0:
        raise ValueError(
            f"No atom found: group={group_name}, resid={target_resid}, "
            f"atomname={target_atomname}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Multiple atoms found ({len(matches)}): group={group_name}, "
            f"resid={target_resid}, atomname={target_atomname}. "
            f"Expected exactly 1. Indices: {[m['index'] for m in matches]}"
        )
    
    return matches[0]['index']


def main():
    args = parse_args()
    try:
        idx = find_pbcatom(args.gro, args.ndx, args.group, args.resid, args.atom)
        print(idx)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
