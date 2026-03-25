#!/usr/bin/env python3
"""Minimal test script to verify gwaslab, bcftools, and tabix are available."""

import sys
import subprocess
import argparse
from importlib.metadata import version
import gwaslab as gl

parser = argparse.ArgumentParser(description="Test gwaslab installation")
parser.add_argument("--output", type=str, required=True, help="Output file prefix")
args = parser.parse_args()

gwaslab_version = version("gwaslab")
print(f"Python version: {sys.version}")
print(f"gwaslab version: {gwaslab_version}")

for tool in ["bcftools", "tabix"]:
    result = subprocess.run([tool, "--version"], capture_output=True, text=True)
    first_line = result.stdout.splitlines()[0] if result.stdout else result.stderr.splitlines()[0]
    print(first_line)

print("All tools loaded successfully!")

output_file = f"{args.output}.munged.sumstats.gz"
with open(output_file, "w") as f:
    f.write(f"gwaslab_version\t{gwaslab_version}\n")
print(f"Wrote stub output: {output_file}")
