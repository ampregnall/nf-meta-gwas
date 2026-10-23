#!/usr/bin/env python3
"""Convert a credible set file to VEP's 5-column headerless input format."""

import argparse
import gwaslab as gl


def parse_args():
    p = argparse.ArgumentParser(description="Format credible set for VEP input")
    p.add_argument("--credset", required=True, help="Credible set TSV from ABF fine-mapping")
    p.add_argument("--output",  required=True, help="Output VEP input file path")
    return p.parse_args()


def main():
    args = parse_args()
    ss = gl.Sumstats(sumstats = args.credset, tab_fmt="\t", fmt="gwaslab", build="38")
    ss.basic_check()
    ss.to_format(path = args.output, fmt = "vep")


if __name__ == "__main__":
    main()
