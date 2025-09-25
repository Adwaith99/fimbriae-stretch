#!/usr/bin/env python3
"""
Builds final pull index and local-numbered posre for a system+variant.

- Step A: make_ndx on clean.pdb with selections -> pull_base.ndx
- Step A.5: rename headers matching the selection text to [ Anchor ] / [ Pulled ]
- Step B: import onto ions.gro with -n; create complement of Protein; post-rename to [ NonProtein ]
- Step C: map Anchor global atoms -> local indices in the anchor molecule .itp; write posre_anchor_local.itp

USAGE: build_indices_and_posres.py <system> <variant_id>
"""
import sys, yaml, subprocess, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG  = yaml.safe_load(open(ROOT/"config.yaml"))
SYSTEM = sys.argv[1]
VARID  = sys.argv[2]

syscfg = next(s for s in CFG["systems"] if s["name"]==SYSTEM)
variant = next(v for v in syscfg["variants"] if v["id"]==VARID)

build = ROOT/f"systems/{SYSTEM}/00_build"
clean = build/"clean.pdb"
ionsg = build/"ions.gro"      # exists after build; else use npt_303K_0.gro for atom table
final_gro = build/"npt_303K_0.gro"
topol_top = build/"topol.top"
anchor_itp = Path(syscfg["anchor_molecule_itp"])

idx_dir = ROOT/f"indices/{SYSTEM}"
idx_dir.mkdir(parents=True, exist_ok=True)
base_ndx = idx_dir/f"{VARID}_pull_base.ndx"
base_named = idx_dir/f"{VARID}_pull_base_named.ndx"
final_ndx = idx_dir/f"{VARID}_pull.ndx"
posre_itp = idx_dir/f"{VARID}_posre_anchor_local.itp"

# ---- Step A: base ndx on clean.pdb ----
selA = f"chain {variant['anchor']['chain']} & r {variant['anchor']['res']}"
selP = f"chain {variant['pulled']['chain']} & r {variant['pulled']['res']}"

if not base_ndx.exists():
    cmd = f"gmx make_ndx -f {clean} -o {base_ndx}"
    inp = f"{selA}\n{selP}\nq\n"
    subprocess.run(cmd.split(), input=inp.encode(), check=True)

# ---- Step A.5: rename headers by match substrings ----
from rename_ndx_by_match import main as _  # allow both direct & import
# If imported main not available, do inline:
def rename_headers(fin, fout, ma, mp):
    text = Path(fin).read_text()
    headers = list(re.finditer(r'(?m)^\[\s*(.+?)\s*\]\s*$', text))
    out=[]; last=0; ra=rp=False
    for m in headers:
        out.append(text[last:m.start()])
        head=m.group(1)
        if (not ra) and (ma in head):
            out.append("[ Anchor ]\n"); ra=True
        elif (not rp) and (mp in head):
            out.append("[ Pulled ]\n"); rp=True
        else:
            out.append(text[m.start():m.end()])
        last=m.end()
    out.append(text[last:])
    Path(fout).write_text(''.join(out))
    if not (ra and rp):
        print("WARN: header rename incomplete (check matches).")
rename_headers(base_ndx, base_named, selA, selP)

# ---- Step B: import onto ions.gro and append NonProtein ----
# Import using -n; then add "protein" and "! protein"
cmd = f"gmx make_ndx -f {ionsg} -n {base_named} -o {final_ndx}"
inp = "protein\n! protein\nq\n"
subprocess.run(cmd.split(), input=inp.encode(), check=True)

# Post-rename the complement header (contains '!' and 'Protein') to [ NonProtein ]
txt = Path(final_ndx).read_text().splitlines()
for i,l in enumerate(txt):
    if l.strip().startswith('[') and '!' in l and 'Protein' in l:
        txt[i] = "[ NonProtein ]"
Path(final_ndx).write_text("\n".join(txt)+"\n")

# ---- Step C: local-numbered posre for anchor molecule ----
# C.1 parse final_ndx for [ Anchor ] atom list (global indices)
def parse_ndx_group(ndx_path, name):
    atoms=[]; cur=None
    for line in Path(ndx_path).read_text().splitlines():
        if line.startswith('['):
            cur=line.strip()[1:-1].strip()
            continue
        if cur==name and line.strip():
            atoms += [int(x) for x in line.split()]
    return atoms
anchor_globals = parse_ndx_group(final_ndx, "Anchor")
if not anchor_globals:
    raise SystemExit("No Anchor atoms found in final ndx.")

# C.2 build global atom table from final_gro
def parse_gro_atoms(gro_path):
    lines = Path(gro_path).read_text().splitlines()
    body = lines[2:-1]  # strip title, natoms, box
    table=[]
    gi=0
    for ln in body:
        gi+=1
        # GRO fixed columns: resnr(5) resname(5) atom(5) atomnr(5) ...
        resnr = int(ln[0:5])
        resnm = ln[5:10].strip()
        atom  = ln[10:15].strip()
        table.append((gi,resnr,resnm,atom))
    return table
globtab = parse_gro_atoms(final_gro)
gmap = {gi:(resnr,resnm,atom) for gi,resnr,resnm,atom in globtab}

# C.3 parse local [ atoms ] from anchor molecule .itp
def parse_itp_atoms(itp_path):
    lines = Path(itp_path).read_text().splitlines()
    in_atoms=False; table=[]
    for ln in lines:
        s=ln.strip()
        if not s or s.startswith(';'): continue
        if s.startswith('['):
            in_atoms = (s.lower().startswith('[ atoms ]'))
            continue
        if in_atoms:
            # fields: nr type resnr resname atom cgnr charge mass
            cols=s.split()
            try:
                nr = int(cols[0]); resnr=int(cols[2]); resnm=cols[3]; atom=cols[4]
            except Exception:
                continue
            table.append((nr,resnr,resnm,atom))
    return table
localtab = parse_itp_atoms(anchor_itp)
# Build lookup by (resnr,resnm,atom) -> local index
lookup={}
for nr,resnr,resnm,atom in localtab:
    lookup.setdefault((resnr,resnm,atom), []).append(nr)

# C.4 map global anchor atoms -> local indices
local_idxs=[]
ambiguous=[]
missing=[]
for gi in anchor_globals:
    key = gmap.get(gi,None)
    if not key:
        missing.append(gi); continue
    cand = lookup.get(key, [])
    if len(cand)==1:
        local_idxs.append(cand[0])
    elif len(cand)>1:
        ambiguous.append((gi,key,cand))
    else:
        # try atom-only fallback (rare)
        miss = True
        for k,v in lookup.items():
            if (k[1]==key[1]) and (k[2]==key[2]):
                local_idxs.append(v[0]); miss=False; break
        if miss: missing.append((gi,key))

if missing:
    print(f"WARN: {len(missing)} anchor atoms not mapped to local indices (check numbering/names).")
if ambiguous:
    print(f"WARN: {len(ambiguous)} ambiguous mappings; picking first local index for each.")

local_idxs = sorted(set(local_idxs))

# C.5 write posre_anchor_local.itp
fcx,fcy,fcz = CFG["globals"]["posre"]["force"]
lines = ["[ position_restraints ]", "; i  funct  fcx  fcy  fcz"]
for li in local_idxs:
    lines.append(f"{li:6d}     1  {fcx}  {fcy}  {fcz}")
Path(posre_itp).write_text("\n".join(lines)+"\n")

# Manifest (optional)
(Path(idx_dir/f"{VARID}_pull_index_manifest.txt")
 .write_text(f"AnchorGlobals={len(anchor_globals)}\nLocalIndices={len(local_idxs)}\n"))
print(f"OK: wrote {final_ndx} and {posre_itp}")
