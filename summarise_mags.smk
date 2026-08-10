include: "common.smk"

CHECKM_DIRNAME = config.get("checkm_outdir_name", "checkm1")
DREP_DIRNAME = config.get("drep_outdir_name", "dRep")
CHECKM_DIR = f"{OUT}/{CHECKM_DIRNAME}"
DREP_DIR = f"{OUT}/{DREP_DIRNAME}"
BINETTE_DIR = f"{OUT}/binette"
BINETTE_RENAMED_DIR = f"{OUT}/binette_renamed_bins"

# rule all must be the first rule in the file (Snakemake's implicit default
# target), so it comes before the helper functions/rules it transitively
# depends on.
rule all:
    input:
        f"{CHECKM_DIR}/checkm.txt",
        f"{OUT}/gtdbtk_classify",
        f"{DREP_DIR}/dRep.genomeInfo",
        f"{DREP_DIR}/dereplicated_genomes"

############################################
# Binette (multi-binner refinement, optional)
############################################

def get_assembly_fasta(wildcards):
    """Return reformatted contigs if anvi_reformat else original assembly
    (same convention as metagenome_binning.smk's get_contigs_for_binning)."""
    assembler = ASSEMBLER_FOR[wildcards.assembly_group]
    base = f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.fa"
    if ANVI_REFORMAT:
        return f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.reformatted.fa"
    return base

def _dir_has_bins(path):
    return os.path.isdir(path) and bool(
        glob.glob(os.path.join(path, "*.fa"))
        + glob.glob(os.path.join(path, "*.fasta"))
        + glob.glob(os.path.join(path, "*.fna"))    # legacy/pre-existing VAMB output (unrenamed)
        + glob.glob(os.path.join(path, "*.fa.gz"))  # all binners' renamed-in-place bins, incl. VAMB
    )

# Optional: extra binner directories for binette --bin_dirs, beyond the
# three this pipeline produces itself. Each entry is a path template
# formatted with output_dir/assembly_group, e.g.:
#   extra_binner_dirs:
#     - "{output_dir}/binnerY/{assembly_group}/bins"
#     - "/absolute/path/binnerZ_results/{assembly_group}"
EXTRA_BINNER_DIRS = config.get("extra_binner_dirs", [])

# Optional: compare bins across separate pipeline runs whose assembly_group
# names don't line up (e.g. a participant's co-assembly run named
# "participantX_co" vs their single-sample run named "participantX") — the
# templated extra_binner_dirs above can't pair those up since it applies the
# same list, substituted with the *current* group's own name, to every group.
# This is a plain per-group dict of already-resolved absolute paths instead:
#   extra_binner_dirs_by_group:
#     participantX_co:
#       - "/abs/path/results_singlesample/metabat2/participantX"
#       - "/abs/path/results_singlesample/vamb/participantX/bins"
#       - "/abs/path/results_singlesample/semibin/participantX/output_bins"
EXTRA_BINNER_DIRS_BY_GROUP = config.get("extra_binner_dirs_by_group", {})

# Optional: auto-derive extra bin dirs for a co-assembly group from a
# *different* pipeline run (e.g. a single-sample run of the same
# participants), instead of listing them by hand. Group names don't need to
# match between runs — samples are matched by fq1 file path (the one thing
# guaranteed to identify "the same physical sample" regardless of naming
# conventions on either side):
#   cross_run_bin_comparisons:
#     spa_co1:                                    # a group in *this* run
#       output_dir:     "/abs/path/results_singlesample"
#       assemblies_tsv: "/abs/path/results_singlesample/../assemblies_singlesample.tsv"
# For each of spa_co1's constituent samples (from *this* run's
# mappings_tsv), we look up that sample's fq1 in the other run's
# assemblies_tsv; whichever assembly_group(s) it matches there are where
# that sample was actually assembled in the other run.
CROSS_RUN_COMPARISONS = config.get("cross_run_bin_comparisons", {})
_cross_run_asm_df_cache = {}

def _cross_run_dirs_for_group(group):
    cfg = CROSS_RUN_COMPARISONS.get(group)
    if not cfg:
        return []
    other_out = cfg["output_dir"]
    other_asm_tsv = cfg["assemblies_tsv"]
    if other_asm_tsv not in _cross_run_asm_df_cache:
        df = pd.read_csv(other_asm_tsv, sep="\t").fillna("None")
        df.columns = [c.strip() for c in df.columns]
        _cross_run_asm_df_cache[other_asm_tsv] = df
    other_asm_df = _cross_run_asm_df_cache[other_asm_tsv]

    dirs = []
    other_groups_matched = set()
    for sample in ASM_TO_MAPPED_SAMPLES.get(group, []):
        r1, _ = MAPPING_READS[sample]
        other_groups_matched |= set(other_asm_df.loc[other_asm_df["fq1"] == r1, "assembly_group"])
    for other_group in other_groups_matched:
        dirs += [
            f"{other_out}/metabat2/{other_group}",
            f"{other_out}/semibin/{other_group}/output_bins",
            f"{other_out}/vamb/{other_group}/bins",
        ]
    return dirs

def get_binette_bin_dirs(wildcards):
    """Return list of binner output dirs for binette --bin_dirs, discovered
    directly on disk. Unlike metagenome_binning.smk (where Binette used to
    live), this workflow runs in a separate Snakemake invocation *after*
    binning has already finished, so there is no checkpoint to query here —
    binning's output directories are just read as plain paths."""
    candidates = [
        f"{OUT}/metabat2/{wildcards.assembly_group}",
        f"{OUT}/semibin/{wildcards.assembly_group}/output_bins",
        f"{OUT}/vamb/{wildcards.assembly_group}/bins",
    ] + [
        tpl.format(output_dir=OUT, assembly_group=wildcards.assembly_group)
        for tpl in EXTRA_BINNER_DIRS
    ] + EXTRA_BINNER_DIRS_BY_GROUP.get(wildcards.assembly_group, []) \
      + _cross_run_dirs_for_group(wildcards.assembly_group)
    return [d for d in candidates if _dir_has_bins(d)]

checkpoint binette_refine:
    input:
        assembly = get_assembly_fasta,
        bin_dirs = get_binette_bin_dirs
    output:
        quality_reports = f"{BINETTE_DIR}/{{assembly_group}}/final_bins_quality_reports.tsv",
        final_bins_dir  = directory(f"{BINETTE_DIR}/{{assembly_group}}/final_bins")
    params:
        outdir     = f"{BINETTE_DIR}/{{assembly_group}}",
        bin_dirs   = lambda w, input: " ".join(input.bin_dirs),
        checkm2_db = config.get("binette_checkm2_db", "")
    conda:
        conda_env("binette", "envs/binette.yaml")
    threads: 8
    resources:
        mem_mb = 32000,
        time   = "12:00:00"
    shell:
        """
        export CHECKM2DB='{params.checkm2_db}'
        binette --bin_dirs {params.bin_dirs} -c {input.assembly} -t {threads} -o {params.outdir}
        """

def get_binette_bin_fasta(wildcards):
    """Locate one Binette-refined bin's FASTA, via the binette_refine
    checkpoint so Snakemake re-evaluates the DAG only after Binette has
    actually run and its final bin count/names are known on disk."""
    final_bins_dir = checkpoints.binette_refine.get(assembly_group=wildcards.assembly_group).output.final_bins_dir
    for ext in ("fa", "fasta"):
        p = os.path.join(final_bins_dir, f"{wildcards.bin_id}.{ext}")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"No {wildcards.bin_id}.{{fa,fasta}} in {final_bins_dir}")

rule rename_binette_bin:
    input:
        fa = get_binette_bin_fasta
    output:
        f"{BINETTE_RENAMED_DIR}/{{assembly_group}}_{{bin_id}}.fa.gz"
    wildcard_constraints:
        # Binette's bin files (e.g. "binette_bin11.fa") contain an underscore
        # themselves, so greedy wildcard matching would otherwise mis-split
        # "spaS144_binette_bin11" as assembly_group="spaS144_binette",
        # bin_id="bin11" instead of assembly_group="spaS144", bin_id="binette_bin11".
        assembly_group = "|".join(re.escape(g) for g in binnable_groups)
    shell:
        """
        mkdir -p {BINETTE_RENAMED_DIR}
        gzip -c {input.fa} > {output}
        """

def get_all_binette_renamed_bins(wildcards):
    """Drives which rename_binette_bin jobs exist. Uses the binette_refine
    checkpoint (not a blind glob) so Snakemake first runs Binette for every
    group, then re-evaluates this function against its real final_bins/
    contents instead of whatever happened to be on disk at DAG-build time."""
    renamed = []
    for group in binnable_groups:
        final_bins_dir = checkpoints.binette_refine.get(assembly_group=group).output.final_bins_dir
        found = glob.glob(os.path.join(final_bins_dir, "*.fa")) + glob.glob(os.path.join(final_bins_dir, "*.fa.gz"))
        for b in found:
            base = os.path.basename(b)
            if base.endswith(".fa.gz") or base.endswith(".fasta.gz"):
                renamed.append(f"{BINETTE_RENAMED_DIR}/{group}_{base}")
            else:
                renamed.append(f"{BINETTE_RENAMED_DIR}/{group}_{base}.gz")
    return renamed

rule aggregate_binette_bins:
    input:
        renamed         = get_all_binette_renamed_bins,
        quality_reports = lambda w: [checkpoints.binette_refine.get(assembly_group=g).output.quality_reports for g in binnable_groups]
    output:
        f"{BINETTE_RENAMED_DIR}/.bins_aggregated.done"
    shell:
        """
        mkdir -p {BINETTE_RENAMED_DIR}
        touch {output}
        """

############################################
# CheckM / dRep / GTDB-Tk
############################################

def bin_source_dir():
    """Bins feeding CheckM/dRep/GTDB-Tk: Binette's refined bins when enabled,
    else the raw per-binner bins produced by metagenome_binning.smk."""
    return BINETTE_RENAMED_DIR if RUN_BINETTE else f"{OUT}/renamed_bins"

def get_all_renamed_bins_local(wildcards):
    """
    Finds preexisting bins inside the renamed directory.
    This enables running this evaluation module directly on precalculated bins.
    """
    bin_dir = bin_source_dir()
    if os.path.exists(bin_dir):
        bins = glob.glob(os.path.join(bin_dir, "*.fa.gz")) + glob.glob(os.path.join(bin_dir, "*.fasta.gz"))
        if bins:
            return bins
    # Fallback to trigger upstream logic if empty
    if RUN_BINETTE:
        return [f"{BINETTE_RENAMED_DIR}/.bins_aggregated.done"]
    return [f"{OUT}/renamed_bins/.bins_aggregated.done"]

rule checkm:
    input:
        bins = get_all_renamed_bins_local
    output:
        report = f"{CHECKM_DIR}/checkm.txt",
        out_dir = directory(CHECKM_DIR)
    conda:
        "envs/checkm.yaml"
    params:
        db_path = config["checkm_db_path"],
        bin_dir = bin_source_dir(),
        ext = config.get("checkm_extension", "gz")
    threads: 16
    shell:
        """
        export CHECKM_DATA_PATH={params.db_path}
        mkdir -p {output.out_dir}
        checkm lineage_wf -t {threads} --pplacer_threads {threads} -x {params.ext} --tab_table -f {output.report} {params.bin_dir} {output.out_dir}
        """

rule generate_drep_info:
    input:
        checkm_report = f"{CHECKM_DIR}/checkm.txt"
    output:
        genome_info = f"{DREP_DIR}/dRep.genomeInfo"
    run:
        import pandas as pd
        df_cm = pd.read_csv(input.checkm_report, sep="\t")
        df_cm.columns = [c.strip() for c in df_cm.columns]
        bin_col = [c for c in df_cm.columns if "bin" in c.lower() or "id" in c.lower()][0]
        comp_col = [c for c in df_cm.columns if "complete" in c.lower()][0]
        cont_col = [c for c in df_cm.columns if "contam" in c.lower()][0]
        
        out_rows = []
        for _, row in df_cm.iterrows():
            bin_id = str(row[bin_col])
            # CheckM strips the whole compound extension (e.g. a bin
            # scanned with -x gz named "foo.fa.gz" is reported as just
            # "foo", not "foo.fa") -- every bin in bin_source_dir() is
            # normalized to ".fa.gz" (see metagenome_binning.smk's in-place
            # per-binner renaming and rename_binette_bin), so re-append
            # exactly that to match the real filename on disk.
            if not bin_id.endswith(".fa.gz"):
                bin_id = f"{bin_id}.fa.gz"
            out_rows.append({
                "genome": bin_id,
                "completeness": row[comp_col],
                "contamination": row[cont_col]
            })
        df_out = pd.DataFrame(out_rows)
        df_out.to_csv(output.genome_info, index=False)

rule drep:
    input:
        bins = get_all_renamed_bins_local,
        genome_info = f"{DREP_DIR}/dRep.genomeInfo"
    output:
        out_dir = directory(f"{DREP_DIR}/dereplicated_genomes")
    conda:
        "envs/drep.yaml"
    params:
        bin_dir = bin_source_dir(),
        drep_out = DREP_DIR,
        comp    = config.get("drep_completeness", 66),
        cont    = config.get("drep_contamination", 20),
        pa      = config.get("drep_primary_ani", 0.95),
        sa      = config.get("drep_secondary_ani", 0.98),
        nc      = config.get("drep_min_overlap", 0.30),
        alg     = config.get("drep_s_algorithm", "ANIn")
    threads: 8
    shell:
        """
        dRep dereplicate {params.drep_out} -g {params.bin_dir}/*.fa.gz --genomeInfo {input.genome_info} -p {threads} -comp {params.comp} -con {params.cont} -pa {params.pa} -sa {params.sa} -nc {params.nc} --S_algorithm {params.alg}
        """

rule gtdbtk:
    input:
        bins = get_all_renamed_bins_local
    output:
        out_dir = directory(f"{OUT}/gtdbtk_classify")
    conda:
        "envs/gtdbtk.yaml"
    params:
        db_path = config["gtdbtk_db_path"],
        bin_dir = bin_source_dir(),
        ext = config.get("gtdbtk_extension", "gz"),
        min_af = config.get("gtdbtk_min_af", 0.5),
        pplacer_cpus = int(config.get("gtdbtk_pplacer_cpus", 8)),
        write_scg_flag = "--write_single_copy_genes" if config.get("gtdbtk_write_single_copy_genes", True) else "",
        keep_intermediates_flag = "--keep_intermediates" if config.get("gtdbtk_keep_intermediates", True) else ""
    threads: 16
    shell:
        """
        export GTDBTK_DATA_PATH={params.db_path}
        mkdir -p {output.out_dir}
        gtdbtk classify_wf \
          --genome_dir {params.bin_dir} \
          -x {params.ext} \
          --out_dir {output.out_dir} \
                    {params.write_scg_flag} \
                    {params.keep_intermediates_flag} \
          --min_af {params.min_af} \
                    --cpus {threads} --pplacer_cpus {params.pplacer_cpus}
        """