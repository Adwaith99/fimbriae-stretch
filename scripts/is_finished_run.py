#!/usr/bin/env python3
import os
import re
import sys


def _read_expected_nsteps(run_dir: str):
    path = os.path.join(run_dir, "expected_nsteps.txt")
    if not os.path.isfile(path):
        return None
    try:
        return int(open(path, "r").read().strip())
    except Exception:
        return None


def _max_step_from_log(log_path: str) -> int:
    max_step = 0
    try:
        with open(log_path, "r", errors="ignore") as fh:
            for line in fh:
                m = re.search(r"(?:^Step\s+|Statistics over\s+)(\d+)", line)
                if m:
                    v = int(m.group(1))
                    if v > max_step:
                        max_step = v
    except Exception:
        return 0
    return max_step


def is_finished_run(run_dir: str) -> bool:
    log_path = os.path.join(run_dir, "pull.log")
    expected_nsteps = _read_expected_nsteps(run_dir)
    if expected_nsteps is None or not os.path.isfile(log_path):
        return False
    max_step = _max_step_from_log(log_path)
    return max_step >= expected_nsteps


def main(argv):
    if len(argv) != 2:
        print("Usage: is_finished_run.py <run_dir>", file=sys.stderr)
        return 2
    run_dir = argv[1]
    return 0 if is_finished_run(run_dir) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
