#!/usr/bin/env python3
"""
rename_ndx_by_match.py <in.ndx> <out.ndx> <match_anchor> <match_pulled>
Renames the first header containing match_anchor -> [ Anchor ]
and the first header containing match_pulled -> [ Pulled ]
Keeps everything else identical (including default groups if present).
"""
import sys, re
if len(sys.argv)!=5:
    print("Usage: rename_ndx_by_match.py in.ndx out.ndx 'chain A & r 128-1200' 'chain D & r 300-360'")
    sys.exit(2)
fin,fout,ma,mp = sys.argv[1:]
text = open(fin).read()
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
open(fout,"w").write(''.join(out))
print(f"Wrote {fout}; AnchorRenamed={ra} PulledRenamed={rp}")
