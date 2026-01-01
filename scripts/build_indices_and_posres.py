#!/usr/bin/env python3
"""
build_indices_and_posres.py <SYSTEM> <VARIANT>   (INDEX-ONLY MODE)

Now ONLY builds index groups:
  A) make_ndx (clean.pdb):
       - "chain <ANCHOR_CHAIN> & r <A_START>-<A_END>" → composite named chX_&_r_a-b
       - "chain <PULLED_CHAIN> & r <P_START>-<P_END>" → composite
       - post-process: rename those headers to [ Anchor ] and [ Pulled ]
  B) make_ndx (npt_final.tpr):
       - find 'Protein' group id
       - create complement "! <protein_id>" and rename last group header to [ NonProtein ]

NO posres generation and NO edits to topol.top anymore.
"""

import os, sys, subprocess, shlex, re, yaml, shutil

def die(msg, code=2):
    print(f"[build-indices] ERROR: {msg}", file=sys.stderr); sys.exit(code)
def info(msg): print(f"[build-indices] {msg}", flush=True)
def run_cmd(cmd, cwd=None, input_str=None):
    info(f"RUN: {' '.join(cmd) if isinstance(cmd,(list,tuple)) else cmd}")
    p = subprocess.run(cmd if isinstance(cmd,(list,tuple)) else shlex.split(cmd),
                       cwd=cwd, input=input_str, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        if p.stdout: print(p.stdout)
        if p.stderr: print(p.stderr, file=sys.stderr)
        die(f"Command failed rc={p.returncode}: {cmd}")
    tail = "\n".join(p.stdout.strip().splitlines()[-8:])
    if tail: print(tail)
    return p.stdout

def find_repo_root():
    p = os.path.abspath(os.getcwd())
    while True:
        if os.path.isfile(os.path.join(p, "config.yaml")): return p
        np = os.path.dirname(p)
        if np == p: return None
        p = np

def parse_range(rng: str):
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", rng)
    if m:
        a,b = int(m.group(1)), int(m.group(2));  return (a,b) if a<=b else (b,a)
    m = re.match(r"^\s*(\d+)\s*$", rng)
    if m: x = int(m.group(1)); return (x,x)
    die(f"Bad residue range: '{rng}'")

def load_cfg(root): return yaml.safe_load(open(os.path.join(root, "config.yaml")))
def get_variant(cfg, system, variant_id):
    for s in cfg["systems"]:
        if s["name"]==system:
            for v in s.get("variants", []):
                if v["id"]==variant_id: return v
    die(f"Variant '{variant_id}' in system '{system}' not found")

def replace_header(ndx_text: str, old: str, new_label: str):
    pat = rf"^\s*\[\s*{re.escape(old)}\s*\]\s*$"
    return re.sub(pat, f"[ {new_label} ]", ndx_text, flags=re.MULTILINE)

def last_group_block(ndx_path: str):
    last=[]; cur=[]
    for line in open(ndx_path):
        if re.match(r"^\s*\[.*\]\s*$", line):
            if cur: last=cur
            cur=[line]
        else:
            cur.append(line)
    if cur: last=cur
    return last

def main():
    if len(sys.argv)!=3:
        print(__doc__); sys.exit(2)
    if not shutil.which("gmx"): die("gmx not found in PATH. Load your module first.")

    system, variant_id = sys.argv[1], sys.argv[2]
    root = find_repo_root() or die("Could not find repo root (config.yaml)")
    cfg = load_cfg(root)
    var = get_variant(cfg, system, variant_id)
    a_chain, a_rng = var["anchor"]["chain"], var["anchor"]["res"]
    p_chain, p_rng = var["pulled"]["chain"], var["pulled"]["res"]
    a1,a2 = parse_range(a_rng)
    # support multiple pulled ranges separated by commas (e.g. "239-341,400-450")
    p_ranges = [r.strip() for r in str(p_rng).split(',') if r.strip()]
    parsed_p_ranges = [parse_range(r) for r in p_ranges]

    build_dir = os.path.join(root, "systems", system, "00_build")
    clean_pdb = os.path.join(build_dir, "clean.pdb")
    if not os.path.isfile(clean_pdb): die(f"clean.pdb not found: {clean_pdb}")
    ndx_path = os.path.join(build_dir, "index.ndx")

    # A) composites on clean.pdb in one go
    info(f"SYSTEM={system} VARIANT={variant_id}")
    info(f"STEP A: make_ndx(clean.pdb) – composite groups")
    # Build make_ndx script: anchor + one or more pulled ranges
    script_lines = [f"chain {a_chain} & r {a1}-{a2}"]
    for pr in parsed_p_ranges:
        script_lines.append(f"chain {p_chain} & r {pr[0]}-{pr[1]}")
    script_lines.append("q")
    script = "\n".join(script_lines) + "\n"
    run_cmd(["gmx","make_ndx","-f",clean_pdb,"-o",ndx_path,"-quiet"], cwd=build_dir, input_str=script)

    # Rename composite headers
    info("STEP A2: renaming composite headers → [ Anchor ] / [ Pulled ]")
    ndx_txt = open(ndx_path).read()
    # Rename Anchor header (single range)
    ndx_txt2 = replace_header(ndx_txt, f"ch{a_chain}_&_r_{a1}-{a2}", "Anchor")

    # For Pulled: find all composite headers created for each parsed range,
    # collect their member indices and merge into a single [ Pulled ] group.
    import itertools
    pulled_block_re = re.compile(rf"^\s*\[\s*ch{re.escape(p_chain)}_&_r_(\d+)-(\d+)\s*\]\s*$", flags=re.MULTILINE)
    lines = ndx_txt2.splitlines(keepends=True)

    # Parse into blocks and collect pulled indices
    blocks = []
    cur = []
    headers = []
    for ln in lines:
        if re.match(r"^\s*\[.*\]\s*$", ln):
            if cur:
                blocks.append((headers[-1] if headers else None, cur))
            cur = [ln]
            hdr = ln.strip()
            headers.append(hdr)
        else:
            cur.append(ln)
    if cur:
        blocks.append((headers[-1] if headers else None, cur))

    pulled_indices = []
    new_blocks = []
    for hdr, blk in blocks:
        if hdr and pulled_block_re.match(hdr):
            # gather indices from blk (skip header line)
            for l in blk[1:]:
                for tok in l.strip().split():
                    if tok.isdigit(): pulled_indices.append(int(tok))
            # skip adding this block (we will replace with merged Pulled)
            continue
        new_blocks.append(blk)

    if pulled_indices:
        # deduplicate while preserving order
        seen = set(); uniq = []
        for x in pulled_indices:
            if x not in seen:
                seen.add(x); uniq.append(x)
        # format indices into lines (15 per line for readability)
        idx_lines = []
        for i in range(0, len(uniq), 15):
            idx_lines.append(" ".join(str(x) for x in uniq[i:i+15]) + "\n")
        pulled_block = [f"[ Pulled ]\n"] + idx_lines
        new_blocks.append(pulled_block)
    else:
        info("WARNING: no pulled composite blocks found; verify index.ndx")

    # write back
    with open(ndx_path, "w") as fh:
        for blk in new_blocks:
            fh.writelines(blk)
    info(f"Wrote {os.path.relpath(ndx_path, build_dir)}")
    info(f"Wrote {os.path.relpath(ndx_path, build_dir)}")

    # B) NonProtein from npt_final.tpr
    info("STEP B: make_ndx(npt_final.tpr) – NonProtein complement of Protein")
    # find npt_final.tpr
    npt_tpr = os.path.join(build_dir, "npt_final.tpr")
    if not os.path.isfile(npt_tpr):
        cand=[]
        for r,_,fs in os.walk(build_dir):
            for f in fs:
                if f.endswith(".tpr") and "npt" in f.lower() and "final" in f.lower():
                    cand.append(os.path.join(r,f))
        if cand: npt_tpr = sorted(cand)[-1]
        else: die("npt_final.tpr not found under 00_build")

    out = run_cmd(["gmx","make_ndx","-f",npt_tpr,"-n",ndx_path,"-o",ndx_path,"-quiet"], cwd=build_dir, input_str="q\n")
    prot_id=None
    for line in out.splitlines():
        m=re.match(r"^\s*(\d+)\s+Protein\b", line)
        if m: prot_id=int(m.group(1)); break
    if prot_id is None:
        info("WARNING: 'Protein' group not found; skipping NonProtein")
        info("DONE."); return
    run_cmd(["gmx","make_ndx","-f",npt_tpr,"-n",ndx_path,"-o",ndx_path,"-quiet"], cwd=build_dir, input_str=f"! {prot_id}\nq\n")
    # rename last block to NonProtein
    block = last_group_block(ndx_path)
    if block:
        if re.match(r"^\s*\[.*\]\s*$", block[0]): block[0]="[ NonProtein ]\n"
        txt = open(ndx_path).read().rstrip()+"\n"+"".join(block)
        open(ndx_path,"w").write(txt)
        info("Added [ NonProtein ]")
    info("DONE.")

if __name__=="__main__":
    main()
