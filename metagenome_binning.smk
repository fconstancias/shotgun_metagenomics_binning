include: "common.smk"

METABAT2_CFG = config.get("metabat2", {})
METABAT2_MIN_CONTIG = int(METABAT2_CFG.get("min_contig", 2000))
METABAT2_MAX_EDGES = int(METABAT2_CFG.get("max_edges", 500))
METABAT2_MIN_CV = float(METABAT2_CFG.get("min_cv", 1.0))
METABAT2_MIN_CV_SUM = float(METABAT2_CFG.get("min_cv_sum", 1.0))
METABAT2_MAX_P = int(METABAT2_CFG.get("max_p", 95))
METABAT2_MIN_S = int(METABAT2_CFG.get("min_s", 60))
METABAT2_MIN_CLS_SIZE = int(METABAT2_CFG.get("min_cls_size", 200000))
METABAT2_SAVE_CLS = bool(METABAT2_CFG.get("save_cls", True))
METABAT2_VERBOSE = bool(METABAT2_CFG.get("verbose", True))

CONCOCT_CLUSTERS = config.get("concoct_clusters", [10, 15, 20])

# SemiBin2's -b (BAM) and -a (strobealign-aemb) inputs are mutually exclusive
# in a single invocation, and -a mode requires SemiBin2-specific split-contig
# abundances (not plain whole-contig depth). So SemiBin2 only runs for groups
# where every mapped sample is bowtie2 (typically co-assemblies), using -b
# directly on the real BAMs; mixed bowtie2+strobealign groups skip SemiBin2.
semibin2_capable_groups = [
    g for g in binnable_groups
    if all(MAPPING_TOOL_FOR[(g, s)] == "bowtie2" for s in ASM_TO_MAPPED_SAMPLES[g])
]

# Anvi'o contigs-db + profile-db, independent of CONCOCT clustering (see
# build_anvio_profile_db below) -- unlike concoct_capable_groups (>=2
# bowtie2 samples, since CONCOCT needs a *merged* profile), this only needs
# >=1: groups with exactly one bowtie2 sample (the common case for a
# single-sample-assembly project, where that one sample is the assembly's
# own self-sample) skip anvi-merge entirely and use that one profile
# directly, since anvi-merge itself refuses to merge a single profile.
anvio_profile_capable_groups = [
    g for g in binnable_groups
    if any(MAPPING_TOOL_FOR[(g, s)] == "bowtie2" for s in ASM_TO_MAPPED_SAMPLES[g])
]

def get_profile_db_target(group):
    """Final anvi'o profile-db target for one group: the merged profile
    (a file -- rule anvi_merge's own declared output) if it has >=2
    bowtie2-mapped samples, else the single sample's own profile directory
    directly (anvi-merge requires at least 2 inputs). Must target rule
    anvi_profile's actual declared output (the directory() itself, not a
    file path inside it) -- targeting "PROFILE_{sample}/PROFILE.db"
    confuses Snakemake's wildcard matching (the {sample_id} wildcard
    greedily swallows the /PROFILE.db suffix, since Snakemake doesn't know
    the directory's real contents, only that it's a directory() output)."""
    bt2_samples = [s for s in ASM_TO_MAPPED_SAMPLES[group] if MAPPING_TOOL_FOR[(group, s)] == "bowtie2"]
    if len(bt2_samples) >= 2:
        return f"{OUT}/concoct/{group}/MERGED/PROFILE.db"
    return f"{OUT}/concoct/{group}/PROFILE_{bt2_samples[0]}"

############################################
# Helpers
############################################

def get_contigs_for_binning(wildcards):
    """Return reformatted contigs if anvi_reformat else original assembly."""
    assembler = ASSEMBLER_FOR[wildcards.assembly_group]
    base = f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.fa"
    if ANVI_REFORMAT:
        return f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.reformatted.fa"
    return base

def get_needed_depth_files(wildcards):
    samples = ASM_TO_MAPPED_SAMPLES[wildcards.assembly_group]
    tools = {MAPPING_TOOL_FOR[(wildcards.assembly_group, s)] for s in samples}
    paths = []
    if "bowtie2" in tools:
        paths.append(f"{OUT}/mapping/{wildcards.assembly_group}/bowtie2_depth.txt")
    # strobealign cross-sample depth: one aemb TSV per sample
    if "strobealign" in tools:
        paths += [f"{OUT}/aemb/{wildcards.assembly_group}/{s}.tsv"
                  for s in samples
                  if MAPPING_TOOL_FOR[(wildcards.assembly_group, s)] == "strobealign"]
    return paths

def get_all_bin_source_dirs(wildcards):
    """Dependencies for aggregate_bins_local: each binner already renames its
    own bins in place (assembly_group + binner baked into the filename), so
    this only needs to force correct scheduling order, not enumerate files —
    the actual discovery happens at run time inside aggregate_bins_local
    itself. MetaBAT2 uses the binning checkpoint (its bin count is dynamic);
    VAMB/SemiBin2 use their own concrete tracked outputs."""
    paths = []
    for group in binnable_groups:
        paths.append(checkpoints.binning.get(assembly_group=group).output.out_dir)
    if "vamb" in BINNERS:
        paths += [f"{OUT}/vamb/{g}/vae_clusters_unsplit.tsv" for g in binnable_groups]
    if "semibin2" in BINNERS:
        paths += [f"{OUT}/semibin/{g}/output_bins/" for g in semibin2_capable_groups]
    return paths

def get_bowtie2_samples(wildcards):
    """Return list of samples mapped with bowtie2 to this assembly_group."""
    return [s for s in ASM_TO_MAPPED_SAMPLES[wildcards.assembly_group]
            if MAPPING_TOOL_FOR[(wildcards.assembly_group, s)] == "bowtie2"]

############################################
# Rule all
############################################

def _all_targets(wildcards):
    targets = []
    if binnable_groups:
        targets.append(f"{OUT}/renamed_bins/.bins_aggregated.done")
    if "semibin2" in BINNERS:
        targets += [f"{OUT}/semibin/{g}/output_bins/" for g in semibin2_capable_groups]
    if "vamb" in BINNERS:
        targets += [f"{OUT}/vamb/{g}/vae_clusters_unsplit.tsv" for g in binnable_groups]
    if RUN_CONCOCT:
        for g in concoct_capable_groups:
            for nc in CONCOCT_CLUSTERS:
                targets.append(f"{OUT}/concoct/{g}/.concoct_{nc}.done")
    # Anvi'o contigs-db + profile-db as their own target, independent of
    # CONCOCT clustering above -- e.g. for downstream anvi'o tools (gene
    # coverage/detection export, DESMAN's SNV extraction, anvi-interactive)
    # that only need a profile-db, not a full CONCOCT collection. Contigs-db
    # itself is built in metagenome_assemble.smk (shared with the optional
    # SCG-taxonomy/stats step there) -- already exists for every group
    # whenever anvi_taxonomy_and_stats: true was set for the assemble run,
    # independent of this flag.
    if config.get("build_anvio_profile_db", False):
        targets += [f"{OUT}/concoct/{g}/.hmms.done" for g in anvio_profile_capable_groups]
        targets += [get_profile_db_target(g) for g in anvio_profile_capable_groups]
    return targets

rule all:
    input: _all_targets

############################################
# Read mapping: bowtie2
############################################

rule bowtie2_build:
    input:
        assembly = get_contigs_for_binning
    output:
        multiext(f"{OUT}/assembly/{{assembler}}/{{assembly_group}}/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        conda_env("mapping", "envs/mapping.yaml")
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "bowtie2-build --threads {threads} {input.assembly} {OUT}/assembly/{wildcards.assembler}/{wildcards.assembly_group}/index"

rule bowtie2_map_and_sort:
    input:
        r1 = lambda w: MAPPING_READS[w.sample_id][0],
        r2 = lambda w: MAPPING_READS[w.sample_id][1],
        index = lambda w: multiext(f"{OUT}/assembly/{ASSEMBLER_FOR[w.assembly_group]}/{w.assembly_group}/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        bam = f"{OUT}/mapping/{{assembly_group}}/{{sample_id}}.bowtie2.sorted.bam",
        bai = f"{OUT}/mapping/{{assembly_group}}/{{sample_id}}.bowtie2.sorted.bam.bai"
    conda:
        conda_env("mapping", "envs/mapping.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    params:
        idx_prefix = lambda w: f"{OUT}/assembly/{ASSEMBLER_FOR[w.assembly_group]}/{w.assembly_group}/index"
    shell:
        """
        bowtie2 -x {params.idx_prefix} -1 {input.r1} -2 {input.r2} --no-unal -p {threads} | samtools view -u -b - | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

############################################
# Depth tables
############################################

rule bowtie2_depth:
    input:
        bams = lambda w: [f"{OUT}/mapping/{w.assembly_group}/{s}.bowtie2.sorted.bam" for s in ASM_TO_MAPPED_SAMPLES[w.assembly_group] if MAPPING_TOOL_FOR[(w.assembly_group, s)] == "bowtie2"],
        bais = lambda w: [f"{OUT}/mapping/{w.assembly_group}/{s}.bowtie2.sorted.bam.bai" for s in ASM_TO_MAPPED_SAMPLES[w.assembly_group] if MAPPING_TOOL_FOR[(w.assembly_group, s)] == "bowtie2"]
    output:
        depth = f"{OUT}/mapping/{{assembly_group}}/bowtie2_depth.txt"
    conda:
        conda_env("binning", "envs/binning.yaml")
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams}"

rule combine_depths:
    input:
        depths   = get_needed_depth_files,
        assembly = get_contigs_for_binning
    output:
        depth = f"{OUT}/mapping/{{assembly_group}}/depth.txt"
    run:
        def parse_assembly_info(fasta_path):
            order = []
            lengths = {}
            current = None
            current_len = 0
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if current is not None:
                            lengths[current] = current_len
                        current = line[1:].strip().split()[0]
                        if current not in lengths:
                            order.append(current)
                        current_len = 0
                    else:
                        current_len += len(line.strip())
            if current is not None:
                lengths[current] = current_len
            return order, lengths

        def normalize_depth_df(depth_path):
            with open(depth_path) as fh:
                first_line = fh.readline()
            header_cols = first_line.rstrip("\n").split("\t")
            # strobealign --aemb output has no header row at all (just
            # "contig\tdepth" per line), unlike jgi_summarize_bam_contig_depths
            # output which has a proper header including 'contigName'.
            if "contigName" in header_cols or "Name" in header_cols:
                df = pd.read_csv(depth_path, sep="\t")
            else:
                sample_name = os.path.splitext(os.path.basename(depth_path))[0]
                df = pd.read_csv(depth_path, sep="\t", header=None, names=["contigName", sample_name])
                # MetaBAT2 requires every sample depth column to be paired with
                # a "<sample>-var" column; strobealign --aemb has no variance
                # estimate, so supply zeros to keep the abundance file format valid.
                df[f"{sample_name}-var"] = 0.0
            # strobealign --aemb uses 'Name' instead of 'contigName'
            if "Name" in df.columns and "contigName" not in df.columns:
                df = df.rename(columns={"Name": "contigName"})
            df = df.drop(columns=["contigLen", "totalAvgDepth"], errors="ignore")
            value_cols = [c for c in df.columns if c != "contigName"]
            if value_cols:
                df[value_cols] = df[value_cols].apply(pd.to_numeric, errors="coerce")
                # Guard against duplicate rows for the same contig from tool-specific length mismatches.
                df = df.groupby("contigName", as_index=False)[value_cols].mean()
            else:
                df = df[["contigName"]].drop_duplicates()
            return df.set_index("contigName")

        asm_order, asm_lengths = parse_assembly_info(input.assembly)
        merged = None
        for depth_path in input.depths:
            dfi = normalize_depth_df(depth_path)
            merged = dfi if merged is None else merged.join(dfi, how="outer")

        if merged is None:
            merged = pd.DataFrame(index=pd.Index([], name="contigName"))

        merged = merged.reindex(asm_order).fillna(0)
        merged.insert(0, "contigLen", [asm_lengths[c] for c in merged.index])

        sample_cols = [c for c in merged.columns if c != "contigLen" and not c.endswith("-var")]
        # Insert at position 1 (right after contigLen) so the final column
        # order matches MetaBAT2's expected "contigName, contigLen,
        # totalAvgDepth, sample pairs..." layout, not "contigName,
        # totalAvgDepth, contigLen, ..." (MetaBAT2 reads columns positionally).
        if sample_cols:
            merged.insert(1, "totalAvgDepth", merged[sample_cols].mean(axis=1))
        else:
            merged.insert(1, "totalAvgDepth", 0)

        merged.reset_index(inplace=True)
        merged.to_csv(output.depth, sep="\t", index=False)

############################################
# MetaBAT2
############################################

checkpoint binning:
    input:
        assembly = get_contigs_for_binning,
        depth    = f"{OUT}/mapping/{{assembly_group}}/depth.txt"
    output:
        out_dir = directory(f"{OUT}/metabat2/{{assembly_group}}")
    wildcard_constraints:
        # Without this, Snakemake can mis-resolve a nested bin file path
        # (e.g. requesting "metabat2/spaS276/bin.6.fa" as an input elsewhere)
        # by greedily binding the whole "spaS276/bin.6.fa" string to the
        # unconstrained assembly_group wildcard, since this directory()
        # output otherwise matches any path nested under it.
        assembly_group = "|".join(re.escape(g) for g in binnable_groups) or "none_placeholder"
    conda:
        conda_env("binning", "envs/binning.yaml")
    threads: 4
    params:
        min_contig  = METABAT2_MIN_CONTIG,
        max_edges   = METABAT2_MAX_EDGES,
        min_cv      = METABAT2_MIN_CV,
        min_cv_sum  = METABAT2_MIN_CV_SUM,
        max_p       = METABAT2_MAX_P,
        min_s       = METABAT2_MIN_S,
        min_cls_size = METABAT2_MIN_CLS_SIZE,
        extra_flags = (" --saveCls" if METABAT2_SAVE_CLS else "") + (" -v" if METABAT2_VERBOSE else "")
    resources:
        mem_mb = 32000,
        time   = "08:00:00"
    shell:
        """
        mkdir -p {output.out_dir}
        metabat2 --numThreads {threads} -i {input.assembly} -a {input.depth} -o {output.out_dir}/bin \
          --minContig {params.min_contig} --maxEdges {params.max_edges} --minCV {params.min_cv} \
          --minCVSum {params.min_cv_sum} --maxP {params.max_p} --minS {params.min_s} \
          --minClsSize {params.min_cls_size}{params.extra_flags}

        # Rename+gzip in place so every bin's file name carries its
        # assembly group and binner (matches the same convention applied to
        # VAMB/SemiBin2 below), and is directly identifiable in downstream
        # CheckM/GTDB-Tk/dRep reports. The --saveCls cluster table (named
        # literally "bin", no extension) is untouched by this glob.
        cd {output.out_dir}
        for f in bin.*.fa; do
            [ -e "$f" ] || continue
            bin_id="${{f#bin.}}"
            bin_id="${{bin_id%.fa}}"
            gzip -c "$f" > "{wildcards.assembly_group}_metabat2_bin.${{bin_id}}.fa.gz"
            rm "$f"
        done
        """

############################################
# Aemb-style coverage TSV (SemiBin2's -a use is gone, but VAMB still needs
# one of these per sample)
############################################

def get_aemb_input(wildcards):
    """Route each (assembly_group, sample_id) pair to its own mapping_tool —
    not just sample_id, since the same sample can use different tools for
    different assembly groups (e.g. bowtie2 for its own assembly, strobealign
    when cross-mapped elsewhere)."""
    tool = MAPPING_TOOL_FOR[(wildcards.assembly_group, wildcards.sample_id)]
    if tool == "bowtie2":
        return {
            "bam": f"{OUT}/mapping/{wildcards.assembly_group}/{wildcards.sample_id}.bowtie2.sorted.bam",
            "bai": f"{OUT}/mapping/{wildcards.assembly_group}/{wildcards.sample_id}.bowtie2.sorted.bam.bai",
        }
    return {
        "contigs": get_contigs_for_binning(wildcards),
        "r1": MAPPING_READS[wildcards.sample_id][0],
        "r2": MAPPING_READS[wildcards.sample_id][1],
    }

rule aemb_tsv:
    """Per (assembly_group, sample_id): aemb-style coverage TSV for VAMB.
    strobealign-mapped samples run --aemb fresh; bowtie2-mapped samples derive
    the same (contig, depth) format from their existing BAM instead of
    mapping a second time (values agree with real strobealign --aemb to
    within ~1-2%)."""
    input:
        unpack(get_aemb_input)
    output:
        f"{OUT}/aemb/{{assembly_group}}/{{sample_id}}.tsv"
    conda:
        conda_env("mapping", "envs/mapping.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        time   = "02:00:00"
    run:
        tool = MAPPING_TOOL_FOR[(wildcards.assembly_group, wildcards.sample_id)]
        if tool == "bowtie2":
            shell("jgi_summarize_bam_contig_depths --outputDepth {output}.raw {input.bam}")
            shell("awk 'BEGIN{{OFS=\"\\t\"}} NR>1 {{print $1, $3}}' {output}.raw > {output}")
            shell("rm {output}.raw")
        else:
            shell("strobealign --aemb -t {threads} {input.contigs} {input.r1} {input.r2} > {output}")

############################################
# SemiBin2
############################################

rule semibin2_bin:
    input:
        contigs = get_contigs_for_binning,
        bams = lambda w: expand(
            f"{OUT}/mapping/{w.assembly_group}/{{sample}}.bowtie2.sorted.bam",
            sample=ASM_TO_MAPPED_SAMPLES[w.assembly_group]
        ),
        bais = lambda w: expand(
            f"{OUT}/mapping/{w.assembly_group}/{{sample}}.bowtie2.sorted.bam.bai",
            sample=ASM_TO_MAPPED_SAMPLES[w.assembly_group]
        )
    output:
        directory(f"{OUT}/semibin/{{assembly_group}}/output_bins/")
    params:
        outdir = f"{OUT}/semibin/{{assembly_group}}",
        # SemiBin2 defaults to 15 epochs, ~1h/epoch even on tiny toy data —
        # override for fast toy/test runs; leave unset for real runs.
        epochs_flag = f"--epochs {config['semibin2_epochs']}" if config.get("semibin2_epochs") else ""
    wildcard_constraints:
        assembly_group = "|".join(re.escape(g) for g in semibin2_capable_groups) or "none_placeholder"
    conda:
        conda_env("semibin2", "envs/semibin2.yaml")
    threads: 8
    resources:
        mem_mb = 32000,
        time   = "12:00:00"
    shell:
        """
        SemiBin2 single_easy_bin -i {input.contigs} --engine cpu -b {input.bams} -o {params.outdir} -p {threads} {params.epochs_flag}

        # Rename in place (already gzipped) so each bin carries its
        # assembly group, matching the convention applied to MetaBAT2/VAMB.
        cd {output}
        for f in SemiBin_*.fa.gz; do
            [ -e "$f" ] || continue
            mv "$f" "{wildcards.assembly_group}_$f"
        done
        """

############################################
# VAMB (latest, --abundance_tsv from merged aemb files)
############################################

rule vamb_bin:
    input:
        contigs    = get_contigs_for_binning,
        aemb_files = lambda w: expand(
            f"{OUT}/aemb/{w.assembly_group}/{{sample}}.tsv",
            sample=ASM_TO_MAPPED_SAMPLES[w.assembly_group]
        )
    output:
        f"{OUT}/vamb/{{assembly_group}}/vae_clusters_unsplit.tsv"
    params:
        outdir    = f"{OUT}/vamb/{{assembly_group}}",
        abund_tsv = f"{OUT}/vamb/{{assembly_group}}_abundance.tsv",
        samples   = lambda w: ASM_TO_MAPPED_SAMPLES[w.assembly_group],
        minfasta  = config.get("vamb_min_fasta_size", 200000)
    conda:
        conda_env("vamb", "envs/vamb.yaml")
    threads: 8
    resources:
        mem_mb = 64000,
        time   = "12:00:00"
    run:
        # VAMB 5.x has no --aemb flag; it wants a single TSV with header
        # "contigname\t<sample1>\t<sample2>..." (--abundance_tsv), so merge
        # the per-sample aemb TSVs. --minfasta is required for VAMB to
        # actually write per-bin FASTA files (into outdir/bins/), which
        # Binette needs later; without it VAMB only writes cluster tables.
        merged = None
        for sample, path in zip(params.samples, input.aemb_files):
            df = pd.read_csv(path, sep="\t", header=None, names=["contigname", sample]).set_index("contigname")
            merged = df if merged is None else merged.join(df, how="outer")
        merged = merged.fillna(0).reset_index()
        os.makedirs(os.path.dirname(params.abund_tsv), exist_ok=True)
        merged.to_csv(params.abund_tsv, sep="\t", index=False)
        shell(
            "rm -rf {params.outdir} && "
            "vamb bin default --outdir {params.outdir} --fasta {input.contigs} "
            "--abundance_tsv {params.abund_tsv} -p {threads} --minfasta {params.minfasta} && "
            # Rename+gzip in place, normalizing extension to .fa.gz (same as
            # MetaBAT2/SemiBin2) so every binner's output can be scanned with
            # one glob/one -x extension downstream — no reason to keep VAMB's
            # own ".fna" convention once it's just being read as FASTA.
            "if [ -d {params.outdir}/bins ]; then "
            "  cd {params.outdir}/bins && "
            "  for f in *.fna; do "
            "    [ -e \"$f\" ] || continue; "
            "    id=\"${{f%.fna}}\"; "
            "    gzip -c \"$f\" > \"{wildcards.assembly_group}_vamb_${{id}}.fa.gz\"; "
            "    rm \"$f\"; "
            "  done; "
            "fi"
        )

############################################
# Aggregate final bins (across all binners)
############################################
# Each binner (metabat2/vamb/semibin2 above) already renames its own bins
# in place as "{assembly_group}_{binner}_..." — no separate per-bin rename
# rule needed. This just symlinks all of them into one directory so
# summarise_mags.smk's non-Binette fallback (run_binette: false) can hand
# CheckM/GTDB-Tk/dRep a single --genome_dir covering every binner, not just
# MetaBAT2.

rule aggregate_bins_local:
    input:
        bin_source_dirs = get_all_bin_source_dirs
    output:
        f"{OUT}/renamed_bins/.bins_aggregated.done"
    run:
        os.makedirs(f"{OUT}/renamed_bins", exist_ok=True)

        def _link_all(pattern):
            for f in glob.glob(pattern):
                dest = os.path.join(f"{OUT}/renamed_bins", os.path.basename(f))
                if not os.path.exists(dest):
                    os.symlink(os.path.abspath(f), dest)

        for group in binnable_groups:
            checkpoint_dir = checkpoints.binning.get(assembly_group=group).output.out_dir
            _link_all(os.path.join(checkpoint_dir, "*_metabat2_bin.*.fa.gz"))

        if "vamb" in BINNERS:
            for group in binnable_groups:
                _link_all(f"{OUT}/vamb/{group}/bins/*_vamb_*.fa.gz")

        if "semibin2" in BINNERS:
            for group in semibin2_capable_groups:
                _link_all(f"{OUT}/semibin/{group}/output_bins/{group}_SemiBin_*.fa.gz")

        with open(output[0], "w"):
            pass

############################################
# CONCOCT binning via Anvi'o (separate track)
############################################
# Contigs DB creation (anvi-gen-contigs-database) now lives in
# metagenome_assemble.smk, shared with the optional SCG-taxonomy/stats step
# there — {OUT}/concoct/{group}/{group}.db is a pre-existing input here.

rule anvi_run_hmms:
    input:
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db"
    output:
        f"{OUT}/concoct/{{assembly_group}}/.hmms.done"
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        # contigs-db may already have some/all default HMM sources from the
        # optional anvi_taxonomy_and_stats step in metagenome_assemble.smk,
        # which runs the same default anvi-run-hmms set. Only request the
        # sources that are actually missing: --just-do-it would delete and
        # redo everything requested, including sources already present
        # (confirmed by reading anvio's TablesForHMMHits.check_sources --
        # just_do_it removes then re-searches the full intersection), which
        # would waste ~6min/group across every already-annotated group for
        # no new result. Requesting only the gap also sidesteps the "some
        # HMM sources are already in the database" ConfigError entirely, so
        # --just-do-it is never needed.
        MISSING=$(python3 -c "
import subprocess
from anvio.terminal import SuppressAllOutput
with SuppressAllOutput():
    import anvio.data.hmm as hmm_data
out = subprocess.run(['sqlite3', '{input.db}', 'select source from hmm_hits_info;'], capture_output=True, text=True).stdout
present = set(l.strip() for l in out.splitlines() if l.strip())
missing = sorted(set(hmm_data.sources.keys()) - present)
print(','.join(missing))
")
        if [ -n "$MISSING" ]; then
            anvi-run-hmms -c {input.db} -T {threads} -I "$MISSING"
        fi
        touch {output}
        """

rule anvi_profile:
    input:
        bam = f"{OUT}/mapping/{{assembly_group}}/{{sample_id}}.bowtie2.sorted.bam",
        db  = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db"
    output:
        profile_dir = directory(f"{OUT}/concoct/{{assembly_group}}/PROFILE_{{sample_id}}")
    params:
        sample_name = lambda w: w.sample_id
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        anvi-profile -i {input.bam} \
            -c {input.db} \
            --output-dir {output.profile_dir} \
            --sample-name {params.sample_name} \
            --num-threads {threads} \
            --skip-hierarchical-clustering
        """

rule anvi_merge:
    input:
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db",
        profiles = lambda w: expand(
            f"{OUT}/concoct/{w.assembly_group}/PROFILE_{{s}}",
            s=[s for s in ASM_TO_MAPPED_SAMPLES[w.assembly_group]
               if MAPPING_TOOL_FOR[(w.assembly_group, s)] == "bowtie2"]
        )
    output:
        merged_dir = directory(f"{OUT}/concoct/{{assembly_group}}/MERGED"),
        merged_db  = f"{OUT}/concoct/{{assembly_group}}/MERGED/PROFILE.db"
    params:
        profile_dbs = lambda w, input: " ".join(f"{p}/PROFILE.db" for p in input.profiles)
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        anvi-merge {params.profile_dbs} \
            -o {output.merged_dir} \
            -c {input.db} \
            --skip-hierarchical-clustering -W
        """

rule anvi_cluster_contigs:
    input:
        merged_db = f"{OUT}/concoct/{{assembly_group}}/MERGED/PROFILE.db",
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db"
    output:
        f"{OUT}/concoct/{{assembly_group}}/.cluster_{{nclust}}.done"
    params:
        nclust = lambda w: w.nclust,
        collection_name = lambda w: f"{w.assembly_group}_concoct_{w.nclust}"
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        anvi-cluster-contigs -p {input.merged_db} \
            -c {input.db} \
            -T {threads} \
            --driver concoct \
            --clusters {params.nclust} \
            -C {params.collection_name} \
            --just-do-it
        touch {output}
        """

rule anvi_summarize:
    input:
        merged_db = f"{OUT}/concoct/{{assembly_group}}/MERGED/PROFILE.db",
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db",
        cluster_done = f"{OUT}/concoct/{{assembly_group}}/.cluster_{{nclust}}.done"
    output:
        summary_dir = directory(f"{OUT}/concoct/{{assembly_group}}/concoct_{{nclust}}_summary"),
        bins_dir = directory(f"{OUT}/concoct/{{assembly_group}}/concoct_{{nclust}}_bins")
    params:
        collection_name = lambda w: f"{w.assembly_group}_concoct_{w.nclust}"
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        """
        anvi-summarize -p {input.merged_db} \
            -c {input.db} \
            -C {params.collection_name} \
            -o {output.summary_dir}

        # Create bins directory (anvi-summarize outputs bins inside summary_dir)
        mkdir -p {output.bins_dir}
        """

rule concoct_done_sentinel:
    input:
        summary_dir = f"{OUT}/concoct/{{assembly_group}}/concoct_{{nclust}}_summary",
        bins_dir = f"{OUT}/concoct/{{assembly_group}}/concoct_{{nclust}}_bins"
    output:
        f"{OUT}/concoct/{{assembly_group}}/.concoct_{{nclust}}.done"
    shell:
        "touch {output}"
