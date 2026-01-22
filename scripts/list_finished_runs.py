#!/usr/bin/env python3
import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from is_finished_run import is_finished_run  # noqa: E402


def find_repo_root(start_dir: str):
    p = os.path.abspath(start_dir)
    while True:
        if os.path.isfile(os.path.join(p, "config.yaml")):
            return p
        parent = os.path.dirname(p)
        if parent == p:
            return None
        p = parent


def resolve_roots(root_arg):
    if root_arg:
        root_abs = os.path.abspath(root_arg)
        if os.path.isfile(os.path.join(root_abs, "config.yaml")):
            repo_root = root_abs
            smd_root = os.path.join(repo_root, "smd")
            return repo_root, smd_root
        if os.path.basename(root_abs) == "smd":
            smd_root = root_abs
            repo_root = find_repo_root(os.path.dirname(root_abs))
            return repo_root, smd_root
        repo_root = find_repo_root(root_abs)
        if repo_root:
            smd_root = os.path.join(repo_root, "smd")
            return repo_root, smd_root
    repo_root = find_repo_root(os.getcwd())
    if not repo_root:
        return None, None
    return repo_root, os.path.join(repo_root, "smd")


def norm_speed_filter(speed: str | None):
    if not speed:
        return None
    if speed.startswith("v"):
        return speed
    try:
        return f"v{float(speed):.3f}"
    except Exception:
        return speed


def norm_rep_filter(rep: str | None):
    if not rep:
        return None
    if rep.startswith("rep"):
        return rep
    if rep.isdigit():
        return f"rep{rep}"
    return rep


def matches_filter(value: str, filt: str | None):
    if not filt:
        return True
    return value == filt


def matches_prefix_equiv(value: str, filt: str | None, prefix: str):
    if not filt:
        return True
    if value == filt:
        return True
    if value.startswith(prefix) and filt.startswith(prefix):
        return value == filt
    if value.startswith(prefix) and not filt.startswith(prefix):
        return value[len(prefix) :] == filt
    if filt.startswith(prefix) and not value.startswith(prefix):
        return filt[len(prefix) :] == value
    return False


def main(argv):
    parser = argparse.ArgumentParser(description="List finished SMD run dirs.")
    parser.add_argument("root", nargs="?", help="Repo root or smd/ directory (optional)")
    parser.add_argument("--system")
    parser.add_argument("--variant")
    parser.add_argument("--speed", help="Speed (e.g., 0.020 or v0.020)")
    parser.add_argument("--rep", help="Replicate (e.g., 1 or rep1)")
    parser.add_argument("--start", help="Start id (e.g., s035 or 035)")
    args = parser.parse_args(argv[1:])

    system = args.system or os.environ.get("SYS") or ""
    variant = args.variant or os.environ.get("VAR") or ""
    speed = args.speed or os.environ.get("SPEED") or ""
    rep = args.rep or os.environ.get("REP") or ""
    start = args.start or os.environ.get("START") or ""

    repo_root, smd_root = resolve_roots(args.root)
    if not repo_root or not smd_root or not os.path.isdir(smd_root):
        print("ERROR: could not locate repo root/config.yaml or smd/ directory.", file=sys.stderr)
        return 2

    speed = norm_speed_filter(speed) if speed else None
    rep = norm_rep_filter(rep) if rep else None
    start = start or None

    pattern = os.path.join(smd_root, "*", "*", "v*", "rep*", "s*")
    candidates = [p for p in sorted(glob(pattern)) if os.path.isdir(p)]

    out = []
    for run_dir in candidates:
        rel = os.path.relpath(run_dir, smd_root)
        parts = rel.split(os.sep)
        if len(parts) != 5:
            continue
        sysname, varname, vtag, rep_tag, start_id = parts
        if system and not matches_filter(sysname, system):
            continue
        if variant and not matches_filter(varname, variant):
            continue
        if speed and not matches_filter(vtag, speed):
            continue
        if rep and not matches_prefix_equiv(rep_tag, rep, "rep"):
            continue
        if start and not matches_prefix_equiv(start_id, start, "s"):
            continue
        if is_finished_run(run_dir):
            out.append(os.path.relpath(run_dir, repo_root))

    for p in out:
        print(p)
    return 0


def glob(pattern):
    import glob as _g
    return _g.glob(pattern)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
