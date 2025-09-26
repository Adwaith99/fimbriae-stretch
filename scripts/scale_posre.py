#!/usr/bin/env python3
"""
Scale fcx/fcy/fcz in a posre.itp and write a new file.
Usage: scale_posre.py <in_posre.itp> <out_posre.itp> <fc_value>
"""
import sys, re, pathlib
if len(sys.argv)!=4:
    print("Usage: scale_posre.py in.itp out.itp 500"); sys.exit(2)
inp, outp, fc = sys.argv[1], sys.argv[2], int(sys.argv[3])
text = pathlib.Path(inp).read_text()
out_lines=[]
in_section=False
for ln in text.splitlines():
    s=ln.strip()
    if s.lower().startswith("[ position_restraints ]"):
        in_section=True
        out_lines.append(ln)
        continue
    if in_section:
        if s.startswith('['):  # next section
            in_section=False
            out_lines.append(ln); continue
        if not s or s.startswith(';'):
            out_lines.append(ln); continue
        cols = s.split()
        # posre rows: i funct fcx fcy fcz
        try:
            i, funct = int(cols[0]), int(cols[1])
            # replace fcx/fcy/fcz by fc
            cols = [cols[0], cols[1], str(fc), str(fc), str(fc)]
            out_lines.append(" ".join(cols))
        except Exception:
            out_lines.append(ln)
    else:
        out_lines.append(ln)
pathlib.Path(outp).write_text("\n".join(out_lines)+"\n")
print(f"Wrote {outp} with fc={fc}")
