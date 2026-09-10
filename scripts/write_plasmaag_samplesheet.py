#!/usr/bin/env python3
"""Write a PlasMAAG --reads_and_assembly_dir / --reads_and_contigs samplesheet:
one row per sample that contributed reads to an assembly group, all pointing
at that group's shared assembly (its own plasmaag_input/ 3-file bundle for
SPAdes, or final.contigs.fa directly for MEGAHIT -- PlasMAAG's own README
distinguishes these two input modes by the third column's name).
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--r1", nargs="+", required=True)
parser.add_argument("--r2", nargs="+", required=True)
parser.add_argument("--assembly-path", required=True,
                     help="plasmaag_input/ dir (SPAdes) or final.contigs.fa (MEGAHIT)")
parser.add_argument("--third-column", required=True, choices=["assembly_dir", "contigs"])
parser.add_argument("--out", required=True)
args = parser.parse_args()

if len(args.r1) != len(args.r2):
    raise ValueError(f"r1/r2 count mismatch: {len(args.r1)} vs {len(args.r2)}")

with open(args.out, "w") as f:
    f.write(f"read1\tread2\t{args.third_column}\n")
    for r1, r2 in zip(args.r1, args.r2):
        f.write(f"{r1}\t{r2}\t{args.assembly_path}\n")
