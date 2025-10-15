#!/usr/bin/env python3
"""
build_indices_and_posres.py <SYSTEM> <VARIANT>

Pipeline:
  A) read config.yaml → anchor/pulled: chain + residue ranges
  B) make_ndx (clean.pdb):
       - commands:
           chain <ANCHOR_CHAIN> & r <A_START>-<A_END>
           chain <PULLED_CHAIN> & r <P_START>-<P_END>
         => creates groups named chX_&_r_<start-end>
       - post-process index.ndx: rename those two headers to [ Anchor ] and [ Pulled ]
  C) make_ndx (npt_final.tpr):
       - find Protein group id
       - run "! <protein_id>" to create complement
       - post-process: rename last group header to [ NonProtein ]
  D) build posre for anchor chain residue range by parsing anchor chain ITP [ atoms ]
       - write systems/<SYSTEM>/00_build/posre_anchor_chain_<CHAIN>.itp
       - force constants from globals.k_kj_mol_nm2
       - insert include right after anchor chain itp include in topol.top

All steps are non-interactive, noisy, and fail-fast.
"""

import os, sys, subprocess, shlex, re, yaml, shutil

# ------------------------ small utils ------------------------

def die(msg, code=2):
    print(f"[build-indices] ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def info(msg):
    print(f"[build-indices] {msg}", flush=True)

def run_cmd(cmd, cwd=None, input_str=None):
    """Run command, capture stdout/stderr, raise on nonzero."""
    info(f"RUN: {' '.join(cmd) if isinstance(cmd, (list,tuple)) else cmd}")
    res = subprocess.run(
        cmd if isinstance(cmd, (list,tuple)) else shlex.split(cmd),
        cwd=cwd, input=input_str, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if res.returncode != 0:
        if res.stdout: print(res.stdout)
        if res.stderr: print(res.stderr, file=sys.stderr)
        die(f"Command failed rc={res.returncode}: {cmd}")
    # print last few lines for visibility
    out = res.stdout.strip().splitlines()
    tail = "\n".join(out[-8:])
    if tail: print(tail)
    return res.stdout

def find_repo_root():
    p = os.path.abspath(os.getcwd())
    while True:
        if os.path.isfile(os.path.join(p, "config.yaml")):
            return p
        np = os.path.dirname(p)
        if np == p:
            return None
        p = np

def parse_range(rng: str):
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", rng)
    if m:
        a,b = int(m.group(1)), int(m.group(2))
        return (a,b) if a<=b else (b,a)
    m = re.match(r"^\s*(\d+)\s*$", rng)
    if m:
        x = int(m.group(1)); return (x,x)
    die(f"Bad residue range: '{rng}'")

# ------------------------ loaders ------------------------

def load_cfg(root):
    with open(os.path.join(root, "config.yaml")) as f:
        return yaml.safe_load(f)

def get_variant(cfg, system, variant_id):
    for s in cfg["systems"]:
        if s["name"] == system:
            for v in s.get("variants", []):
                if v["id"] == variant_id:
                    return v
    die(f"Variant '{variant_id}' in system '{system}' not found in config.yaml")

# ------------------------ index helpers ------------------------

def replace_header(ndx_text: str, old: str, new_label: str):
    """Replace header [ old ] → [ new_label ] (exact match)."""
    pattern = rf"^\s*\[\s*{re.escape(old)}\s*\]\s*$"
    return re.sub(pattern, f"[ {new_label} ]", ndx_text, flags=re.MULTILINE)

def last_group_block(ndx_path: str):
    last = []
    curr = []
    with open(ndx_path) as f:
        for line in f:
            if re.match(r"^\s*\[.*\]\s*$", line):
                if curr: last = curr
                curr = [line]
            else:
                curr.append(line)
    if curr: last = curr
    return last

def write_text(path, text):
    with open(path, "w") as f:
        f.write(text)

# ------------------------ posre helpers ------------------------

def parse_atoms_from_itp(itp_text: str):
    """Return list of (ai, resnr) from [ atoms ] section."""
    in_atoms = False
    out = []
    for line in itp_text.splitlines():
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        if s.startswith("["):
            in_atoms = "atoms" in s
            continue
        if in_atoms:
            parts = s.split()
            if len(parts) >= 6 and parts[0].isdigit() and parts[2].isdigit():
                ai    = int(parts[0])
                resnr = int(parts[2])
                out.append((ai,resnr))
    return out

def insert_after_anchor_include(top_path: str, anchor_itp_basename: str, posre_basename: str):
    txt = open(top_path,"r").read()
    inc_line = f'#include "{anchor_itp_basename}"'
    posre_line = f'#include "{posre_basename}"'
    if posre_line in txt:
        info(f"posre include already present in {os.path.basename(top_path)}")
        return
    lines = txt.splitlines()
    out = []
    inserted = False
    for line in lines:
        out.append(line)
        if not inserted and line.strip() == inc_line:
            out.append(posre_line)
            inserted = True
    if not inserted:
        info(f"WARNING: anchor include line not found; appending posre include at end.")
        out.append(posre_line)
    write_text(top_path, "\n".join(out) + "\n")
    info(f"Inserted posre include → {posre_basename}")

# ------------------------ main ------------------------

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)

    if not shutil.which("gmx"):
        die("gmx not found in PATH. Load your configured module first.")

    system, variant_id = sys.argv[1], sys.argv[2]
    root = find_repo_root()
    if not root:
        die("Could not locate repo root (config.yaml). Run from inside the repo.")

    cfg = load_cfg(root)
    G = cfg["globals"]
    SMD = G.get("smd", {})
    k_kj = float(G.get("k_kj_mol_nm2", 100.0))

    var = get_variant(cfg, system, variant_id)
    a_chain = var["anchor"]["chain"]
    a_rng   = var["anchor"]["res"]
    p_chain = var["pulled"]["chain"]
    p_rng   = var["pulled"]["res"]

    a1,a2 = parse_range(a_rng)
    p1,p2 = parse_range(p_rng)

    build_dir = os.path.join(root, "systems", system, "00_build")
    clean_pdb = os.path.join(build_dir, "clean.pdb")
    if not os.path.isfile(clean_pdb):
        die(f"clean.pdb not found: {clean_pdb}")
    topol_top = os.path.join(build_dir, "topol.top")
    if not os.path.isfile(topol_top):
        die(f"topol.top not found: {topol_top}")

    # Anchor chain itp from template
    tpl = SMD.get("anchor_chain_template", "systems/{system}/00_build/topol_Protein_chain_{chain}.itp")
    anchor_itp = os.path.join(root, tpl.format(system=system, chain=a_chain))
    if not os.path.isfile(anchor_itp):
        die(f"Anchor chain ITP not found: {anchor_itp}")

    # --- A) make_ndx on clean.pdb to build composite groups in one go
    ndx_path = os.path.join(build_dir, "index.ndx")
    info(f"SYSTEM={system} VARIANT={variant_id}")
    info(f"Anchor: chain {a_chain} & r {a1}-{a2}")
    info(f"Pulled: chain {p_chain} & r {p1}-{p2}")
    info("STEP A: make_ndx (clean.pdb) → composite groups")
    cmd = ["gmx","make_ndx","-f",clean_pdb,"-o",ndx_path,"-quiet"]
    script = f"""chain {a_chain} & r {a1}-{a2}
chain {p_chain} & r {p1}-{p2}
q
"""
    run_cmd(cmd, cwd=build_dir, input_str=script)

    # Rename composite headers to Anchor / Pulled
    anchor_hdr = f"ch{a_chain}_&_r_{a1}-{a2}"
    pulled_hdr = f"ch{p_chain}_&_r_{p1}-{p2}"
    info("STEP A2: renaming composite headers → [ Anchor ], [ Pulled ]")
    ndx_txt = open(ndx_path,"r").read()
    before = ndx_txt
    ndx_txt = replace_header(ndx_txt, anchor_hdr, "Anchor")
    ndx_txt = replace_header(ndx_txt, pulled_hdr, "Pulled")
    # Warn if not found
    if before == ndx_txt:
        info("WARNING: Could not find composite headers exactly; keeping file but you should verify.")
    write_text(ndx_path, ndx_txt)
    info(f"Wrote {os.path.relpath(ndx_path, build_dir)}")

    # --- B) make_ndx on npt_final.tpr to add NonProtein
    # We prefer the canonical name you used earlier:
    npt_tpr_guess = os.path.join(build_dir, "npt_final.tpr")
    if not os.path.isfile(npt_tpr_guess):
        # try a loose search within 00_build
        matches = []
        for rootdir,dirs,files in os.walk(build_dir):
            for f in files:
                if f.endswith(".tpr") and "npt" in f.lower() and "final" in f.lower():
                    matches.append(os.path.join(rootdir,f))
        if matches:
            npt_tpr_guess = sorted(matches)[-1]
        else:
            die("Could not locate npt_final.tpr under 00_build; required for NonProtein creation.")

    info(f"STEP B: make_ndx (npt_final.tpr) → find Protein group id")
    # First run: list groups (we capture stdout) and quit
    out = run_cmd(["gmx","make_ndx","-f",npt_tpr_guess,"-n",ndx_path,"-o",ndx_path,"-quiet"],
                  cwd=build_dir, input_str="q\n")
    # Parse Protein group id from printed group table
    protein_id = None
    for line in out.splitlines():
        # format:   1 Protein     : (atoms...)
        m = re.match(r"^\s*(\d+)\s+Protein\b", line)
        if m:
            protein_id = int(m.group(1)); break
    if protein_id is None:
        info("WARNING: Could not find 'Protein' group in listing; NonProtein will be skipped.")
    else:
        info(f"Protein group id = {protein_id}")
        info("STEP B2: create complement → NonProtein")
        run_cmd(["gmx","make_ndx","-s",npt_tpr_guess,"-n",ndx_path,"-o",ndx_path,"-quiet"],
                cwd=build_dir, input_str=f"! {protein_id}\nq\n")
        # Rename last group to NonProtein
        block = last_group_block(ndx_path)
        if block:
            if re.match(r"^\s*\[.*\]\s*$", block[0]):
                block[0] = "[ NonProtein ]\n"
            ndx_all = open(ndx_path,"r").read().rstrip() + "\n" + "".join(block)
            write_text(ndx_path, ndx_all)
            info("Added [ NonProtein ]")

    # --- C) Build posre for anchor chain residue range
    info("STEP C: generate posre for anchor chain residue range from anchor chain ITP")
    itp_txt = open(anchor_itp,"r").read()
    atoms = parse_atoms_from_itp(itp_txt)
    if not atoms:
        die(f"No atoms parsed from [ atoms ] in {os.path.basename(anchor_itp)}")
    sel_ai = [ai for (ai,resnr) in atoms if a1 <= resnr <= a2]
    if not sel_ai:
        die(f"No atoms found in residue range {a1}-{a2} for chain {a_chain} in anchor ITP")
    posre_path = os.path.join(build_dir, f"posre_anchor_chain_{a_chain}.itp")
    with open(posre_path,"w") as f:
        f.write("[ position_restraints ]\n")
        f.write("; ai  funct  fc_x   fc_y   fc_z\n")
        for ai in sel_ai:
            f.write(f"{ai:6d}   1   {k_kj:.3f}  {k_kj:.3f}  {k_kj:.3f}\n")
    info(f"Wrote {os.path.relpath(posre_path, build_dir)} ({len(sel_ai)} atoms; k={k_kj:.3f})")

    # Insert posre include after anchor itp in topol.top
    anchor_itp_base = os.path.basename(anchor_itp)
    posre_base = os.path.basename(posre_path)
    insert_after_anchor_include(topol_top, anchor_itp_base, posre_base)
    info("DONE.")

if __name__ == "__main__":
    main()
