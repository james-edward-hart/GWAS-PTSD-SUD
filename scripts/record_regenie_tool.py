#!/usr/bin/env python3

import argparse
import hashlib
import os
import shutil
import subprocess
import sys


def die(message):
    raise SystemExit(f"ERROR: {message}")


def resolve_executable(tool):
    if not tool:
        die("regenie tool path is empty")
    resolved = tool if os.sep in tool else shutil.which(tool)
    if not resolved or not os.path.exists(resolved):
        die(f"regenie executable not found: {tool}")
    if not os.access(resolved, os.X_OK):
        die(f"regenie is not executable: {resolved}")
    return os.path.realpath(resolved)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def regenie_version(path):
    proc = subprocess.run([path, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        die(f"regenie --version failed with exit code {proc.returncode}")
    return " | ".join(line.strip() for line in proc.stdout.splitlines() if line.strip()) or "NA"


def main():
    parser = argparse.ArgumentParser(description="Record regenie executable provenance.")
    parser.add_argument("--tool", default="regenie", help="regenie executable name or path")
    parser.add_argument("--out", required=True, help="output key/value TSV")
    args = parser.parse_args()

    resolved = resolve_executable(args.tool)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    rows = [
        ("tool_regenie_path", resolved),
        ("tool_regenie_sha256", sha256_file(resolved)),
        ("tool_regenie_version", regenie_version(resolved)),
    ]
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
