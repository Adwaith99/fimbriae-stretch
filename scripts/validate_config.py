#!/usr/bin/env python3
import sys, yaml
try:
    cfg=yaml.safe_load(open("config.yaml"))
except Exception as e:
    print("YAML error:", e); sys.exit(1)
assert "globals" in cfg and "systems" in cfg, "config.yaml must have globals and systems"
print("config.yaml OK")
