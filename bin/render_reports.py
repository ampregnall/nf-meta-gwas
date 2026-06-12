#!/usr/bin/env python3
"""Generate QC and results Quarto HTML reports for a single phenotype.

Writes a params.yml, copies both .qmd templates into the working directory,
then invokes quarto render for each report.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def parse_args():
    p = argparse.ArgumentParser(description="Render Quarto HTML reports")
    p.add_argument("--phenotype", required=True)
    p.add_argument("--results_dir", required=True, help="Path to top-level results/ directory")
    p.add_argument("--locus_annotated", required=True)
    p.add_argument("--credsets", nargs="+", required=True)
    p.add_argument("--filter_stats", required=True)
    p.add_argument("--ldsc_h2", required=True)
    p.add_argument("--has_vep", action="store_true", default=False)
    p.add_argument("--has_hyprcoloc", action="store_true", default=False)
    p.add_argument("--has_flames", action="store_true", default=False)
    p.add_argument("--templates_dir", required=True,
                   help="Directory containing qc_report.qmd and results_report.qmd")
    p.add_argument("--output_prefix", required=True)
    return p.parse_args()


def render(template_name: str, params_file: Path, output_name: str) -> None:
    cmd = [
        "quarto", "render", template_name,
        "--execute-params", str(params_file),
        "--output", output_name,
        "--to", "html",
    ]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"ERROR: quarto render failed for {template_name}", file=sys.stderr)
        sys.exit(result.returncode)


def main():
    args = parse_args()

    params = {
        "phenotype": args.phenotype,
        "results_dir": str(Path(args.results_dir).resolve()),
        "locus_annotated": str(Path(args.locus_annotated).resolve()),
        "credsets": [str(Path(f).resolve()) for f in args.credsets],
        "filter_stats": str(Path(args.filter_stats).resolve()),
        "ldsc_h2": str(Path(args.ldsc_h2).resolve()),
        "has_vep": args.has_vep,
        "has_hyprcoloc": args.has_hyprcoloc,
        "has_flames": args.has_flames,
    }

    params_file = Path("params.yml")
    with open(params_file, "w") as f:
        yaml.dump(params, f, default_flow_style=False)

    templates_dir = Path(args.templates_dir)
    prefix = args.output_prefix

    for template, outfile in [
        ("qc_report.qmd",      f"{prefix}.qc.html"),
        ("results_report.qmd", f"{prefix}.results.html"),
    ]:
        shutil.copy(templates_dir / template, template)
        render(template, params_file, outfile)


if __name__ == "__main__":
    main()
