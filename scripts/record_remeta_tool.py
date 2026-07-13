#!/usr/bin/env python3

import argparse
import hashlib
import os
import shutil
import subprocess


def die(message):
    raise SystemExit(f"ERROR: {message}")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Record ReMeta executable provenance.")
    parser.add_argument("--tool", default="remeta")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    resolved = args.tool if os.sep in args.tool else shutil.which(args.tool)
    if not resolved or not os.path.isfile(resolved) or not os.access(resolved, os.X_OK):
        die(f"ReMeta executable not found or not executable: {args.tool}")
    resolved = os.path.realpath(resolved)
    proc = subprocess.run([resolved, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        die(f"remeta --version failed with exit code {proc.returncode}")
    version = " | ".join(line.strip() for line in proc.stdout.splitlines() if line.strip()) or "NA"
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        handle.write(f"tool_remeta_path\t{resolved}\n")
        handle.write(f"tool_remeta_sha256\t{sha256_file(resolved)}\n")
        handle.write(f"tool_remeta_version\t{version}\n")


if __name__ == "__main__":
    main()
