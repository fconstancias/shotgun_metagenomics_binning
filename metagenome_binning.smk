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

############################################
# Helpers
############################################

def get_contigs_for_binning(wildcards):
    """Return reformatted contigs if anvi_reformat else original assembly."""
    base = f"{OUT}/assembly/{ASSEMBLER}/{wildcards.assembly_group}/final.contigs.fa"
    if ANVI_REFORMAT:
        return f"{OUT}/assembly/{ASSEMBLER}/{wildcards.assembly_group}/final.contigs.reformatted.fa"
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

def get_all_renamed_bins(wildcards):
    renamed_bins = []
    for group in all_groups:
        if group not in ASM_TO_MAPPED_SAMPLES:
            continue
        checkpoint_dir = checkpoints.binning.get(assembly_group=group).output.out_dir
        found_bins = glob.glob(os.path.join(checkpoint_dir, "bin.*.fa")) + \
                     glob.glob(os.path.join(checkpoint_dir, "bin.*.fasta")) + \
                     glob.glob(os.path.join(checkpoint_dir, "bin.*.fa.gz"))
        for b in found_bins:
            bin_base = os.path.basename(b)
            if bin_base.endswith(".fa.gz") or bin_base.endswith(".fasta.gz"):
                renamed_bins.append(f"{OUT}/renamed_bins/{group}_{bin_base}")
            else:
                renamed_bins.append(f"{OUT}/renamed_bins/{group}_{bin_base}.gz")
    return renamed_bins

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
        targets += [f"{OUT}/semibin/{g}/output_recluster_bins/" for g in semibin2_capable_groups]
    if "vamb" in BINNERS:
        targets += [f"{OUT}/vamb/{g}/vae_clusters_unsplit.tsv" for g in binnable_groups]
    if RUN_CONCOCT:
        for g in concoct_capable_groups:
            for nc in CONCOCT_CLUSTERS:
                targets.append(f"{OUT}/concoct/{g}/.concoct_{nc}.done")
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
        multiext(f"{OUT}/assembly/{ASSEMBLER}/{{assembly_group}}/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    conda:
        conda_env("mapping", "envs/mapping.yaml")
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "bowtie2-build --threads {threads} {input.assembly} {OUT}/assembly/{ASSEMBLER}/{wildcards.assembly_group}/index"

rule bowtie2_map_and_sort:
    input:
        r1 = lambda w: MAPPING_READS[w.sample_id][0],
        r2 = lambda w: MAPPING_READS[w.sample_id][1],
        index = lambda w: multiext(f"{OUT}/assembly/{ASSEMBLER}/{w.assembly_group}/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        bam = f"{OUT}/mapping/{{assembly_group}}/{{sample_id}}.bowtie2.sorted.bam",
        bai = f"{OUT}/mapping/{{assembly_group}}/{{sample_id}}.bowtie2.sorted.bam.bai"
    conda:
        conda_env("mapping", "envs/mapping.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 240
    shell:
        """
        bowtie2 -x {OUT}/assembly/{ASSEMBLER}/{wildcards.assembly_group}/index -1 {input.r1} -2 {input.r2} --no-unal -p {threads} | samtools view -u -b - | samtools sort -@ {threads} -o {output.bam}
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
        out_dir = directory(f"{OUT}/bins/{{assembly_group}}")
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
        directory(f"{OUT}/semibin/{{assembly_group}}/output_recluster_bins/")
    params:
        outdir = f"{OUT}/semibin/{{assembly_group}}"
    wildcard_constraints:
        assembly_group = "|".join(re.escape(g) for g in semibin2_capable_groups) or "none_placeholder"
    conda:
        conda_env("semibin2", "envs/semibin2.yaml")
    threads: 8
    resources:
        mem_mb = 32000,
        time   = "12:00:00"
    shell:
        "SemiBin2 single_easy_bin -i {input.contigs} --engine cpu -b {input.bams} -o {params.outdir} -p {threads}"

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
            "--abundance_tsv {params.abund_tsv} -p {threads} --minfasta {params.minfasta}"
        )

############################################
# Rename and aggregate final bins
############################################

rule rename_bin:
    input:
        fa = lambda w: next(
            p for p in [f"{OUT}/bins/{w.assembly_group}/{w.bin_id}.{ext}" for ext in ("fa", "fasta")]
            if os.path.exists(p)
        )
    output:
        f"{OUT}/renamed_bins/{{assembly_group}}_{{bin_id}}.fa.gz"
    wildcard_constraints:
        assembly_group = "|".join(re.escape(g) for g in binnable_groups)
    shell:
        """
        mkdir -p {OUT}/renamed_bins
        gzip -c {input.fa} > {output}
        """

rule aggregate_bins_local:
    input:
        renamed = get_all_renamed_bins
    output:
        f"{OUT}/renamed_bins/.bins_aggregated.done"
    shell:
        """
        mkdir -p {OUT}/renamed_bins
        touch {output}
        """

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
        anvi-run-hmms -c {input.db} -T {threads}
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
