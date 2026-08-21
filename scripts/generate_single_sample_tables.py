#!/usr/bin/env python3
"""
Generate assemblies.tsv / mappings.tsv for single-sample assembly groups,
given a list of participants.

For each Tobin in {Yes, Potential} sample of the requested participants:
  - one assembly_group (named spa{sample_id without underscore}), reusing an
    existing metaSPAdes assembly under --prebuilt-dir if present, else left
    for metagenome_assemble.smk to assemble from raw reads.
  - mapping rows: bowtie2 for the assembly's own sample, strobealign for
    same-participant chronological neighbors and for --other-pool-n samples
    drawn from distinct OTHER participants (reproducible via --seed).

    Same-participant window: normally --window before + --window after
    (capped by how many that participant actually has on each side). If
    one side is completely empty (the sample is that participant's first
    or last chronologically), the other side gets --window + 1 instead of
    --window, capped by availability -- a small compensating bump, not an
    attempt to force every sample to the same total neighbor count (a
    participant with few total samples simply has fewer to draw from,
    that's not a bug to work around).

    Other-participant pool: NOT a single fixed pool shared identically
    across every group (that breaks once every participant is a target --
    nobody would be left "outside" to draw a shared background pool from,
    and even patched, participants who land in that shared pool would get
    fewer than everyone else for their own assembly's mapping table).
    Instead each assembly group independently draws --other-pool-n samples
    from --other-pool-n distinct other participants (excluding itself, so
    always exactly --other-pool-n, no edge case) -- still fully
    reproducible via --seed, just not literally identical across groups.
"""
import argparse
import os
import random
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tobin-tsv", default="/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/to_bin_sample_host_date.tsv")
    p.add_argument("--prebuilt-dir", default="/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/B02_spades_single_fixed")
    p.add_argument("--raw-reads-dir", default="/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/002_HostReadRemoval")
    p.add_argument("--participants", required=True, help="Comma-separated list of Subject ids to build assemblies for, e.g. 292,431,470")
    p.add_argument("--window", type=int, default=3, help="Chronological same-participant neighbors before/after to cross-map (default 3)")
    p.add_argument("--other-pool-n", type=int, default=5, help="Number of distinct other participants (1 sample each), drawn independently per assembly group (default 5)")
    p.add_argument("--seed", type=int, default=42, help="Random seed for the other-participant pool (default 42)")
    p.add_argument("--assembler-abbrev", default="spa")
    p.add_argument("--out-assemblies", required=True)
    p.add_argument("--out-mappings", required=True)
    return p.parse_args()


def raw_fastqs(raw_reads_dir, sample_id):
    d = os.path.join(raw_reads_dir, sample_id)
    return (
        os.path.join(d, f"{sample_id}_hostremoved_R1.fastq.gz"),
        os.path.join(d, f"{sample_id}_hostremoved_R2.fastq.gz"),
    )


def prebuilt_assembly_path(prebuilt_dir, participant, sample_id, abbrev):
    tag = f"{abbrev}{sample_id.replace('_', '')}"
    path = os.path.join(prebuilt_dir, str(participant), sample_id, f"{tag}scaffolds_fixed.fasta")
    return path if os.path.isfile(path) else None


def main():
    args = parse_args()
    targets = [s.strip() for s in args.participants.split(",") if s.strip()]

    df = pd.read_csv(args.tobin_tsv, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df["Time"] = pd.to_datetime(df["Time"], format="%d/%m/%Y")
    df = df[df["Tobin"].isin(["Yes", "Potential"])].copy()
    df["Subject"] = df["Subject"].astype(str)

    missing_raw = [s for s in df["Sample"] if not os.path.isfile(raw_fastqs(args.raw_reads_dir, s)[0])]
    if missing_raw:
        sys.exit(f"ERROR: missing raw fastqs for samples: {missing_raw}")

    # Per-participant chronologically ordered sample lists (all Tobin samples, for windowing)
    per_participant = {
        subj: grp.sort_values("Time")["Sample"].tolist()
        for subj, grp in df.groupby("Subject")
    }

    rng = random.Random(args.seed)

    asm_rows = []
    map_rows = []
    n_prebuilt, n_new = 0, 0

    for subj in targets:
        if subj not in per_participant:
            sys.exit(f"ERROR: participant {subj} has no Tobin Yes/Potential samples in {args.tobin_tsv}")
        ordered = per_participant[subj]
        for i, sample_id in enumerate(ordered):
            group = f"{args.assembler_abbrev}{sample_id.replace('_', '')}"
            fq1, fq2 = raw_fastqs(args.raw_reads_dir, sample_id)

            prebuilt = prebuilt_assembly_path(args.prebuilt_dir, subj, sample_id, args.assembler_abbrev)
            if prebuilt:
                asm_rows.append({"assembly_group": group, "assembly_path": prebuilt, "fq1": "None", "fq2": "None"})
                n_prebuilt += 1
            else:
                asm_rows.append({"assembly_group": group, "assembly_path": "None", "fq1": fq1, "fq2": fq2})
                n_new += 1

            map_rows.append({"sample_id": sample_id, "assembly_group": group, "fq1": fq1, "fq2": fq2, "mapping_tool": "bowtie2"})

            # Same-participant chronological window: normally `window` before
            # + `window` after, capped by availability on each side. If one
            # side is completely empty (this is that participant's first or
            # last sample), bump the other side to `window + 1` instead --
            # a small compensating nudge, not an attempt to force every
            # sample to the same total (participants with few samples just
            # have fewer neighbors available, full stop).
            n_before_avail, n_after_avail = i, len(ordered) - 1 - i
            if n_before_avail == 0:
                n_before, n_after = 0, min(args.window + 1, n_after_avail)
            elif n_after_avail == 0:
                n_before, n_after = min(args.window + 1, n_before_avail), 0
            else:
                n_before, n_after = min(args.window, n_before_avail), min(args.window, n_after_avail)

            neighbor_indices = list(range(i - n_before, i)) + list(range(i + 1, i + 1 + n_after))
            for j in neighbor_indices:
                nb_id = ordered[j]
                nb_fq1, nb_fq2 = raw_fastqs(args.raw_reads_dir, nb_id)
                map_rows.append({"sample_id": nb_id, "assembly_group": group, "fq1": nb_fq1, "fq2": nb_fq2, "mapping_tool": "strobealign"})

            # Other-participant pool: drawn independently for THIS group,
            # excluding this group's own participant -- always exactly
            # other_pool_n distinct other participants, no edge case,
            # reproducible via --seed (rng state advances per group, so a
            # full re-run is deterministic even though each group's draw
            # differs from every other group's).
            other_candidates = [p for p in per_participant if p != subj]
            other_subjs = rng.sample(other_candidates, args.other_pool_n)
            for other_subj in other_subjs:
                other_sample = rng.choice(per_participant[other_subj])
                o_fq1, o_fq2 = raw_fastqs(args.raw_reads_dir, other_sample)
                map_rows.append({"sample_id": other_sample, "assembly_group": group, "fq1": o_fq1, "fq2": o_fq2, "mapping_tool": "strobealign"})

    pd.DataFrame(asm_rows).to_csv(args.out_assemblies, sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(args.out_mappings, sep="\t", index=False)

    print(f"Wrote {len(asm_rows)} assembly groups ({n_prebuilt} pre-built, {n_new} to be SPAdes-assembled) -> {args.out_assemblies}")
    print(f"Wrote {len(map_rows)} mapping rows -> {args.out_mappings}")

    map_df = pd.DataFrame(map_rows)
    rows_per_group = map_df.groupby("assembly_group").size()
    for subj in targets:
        n = len(per_participant[subj])
        groups = [f"{args.assembler_abbrev}{s.replace('_', '')}" for s in per_participant[subj]]
        counts = rows_per_group[groups]
        print(f"  participant {subj}: {n} groups, {counts.min()}-{counts.max()} mapping rows/group (exact, edge-adaptive)")


if __name__ == "__main__":
    main()
