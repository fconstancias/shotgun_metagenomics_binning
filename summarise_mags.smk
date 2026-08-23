include: "common.smk"

CHECKM_DIRNAME = config.get("checkm_outdir_name", "checkm1")
DREP_DIRNAME = config.get("drep_outdir_name", "dRep")
CHECKM_DIR = f"{OUT}/{CHECKM_DIRNAME}"
DREP_DIR = f"{OUT}/{DREP_DIRNAME}"
BINETTE_DIR = f"{OUT}/binette"
BINETTE_RENAMED_DIR = f"{OUT}/binette_renamed_bins"

# Binette already computes CheckM2 completeness/contamination for every
# winning bin (final_bins_quality_reports.tsv, per group) -- when true,
# dRep's genomeInfo is built directly from those instead of running CheckM1
# a second time on the same bins. Only takes effect when run_binette is
# also true (nothing to reuse otherwise); ignored/false by default so
# existing configs keep today's CheckM1-always-runs behavior unchanged.
USE_BINETTE_CHECKM2 = config.get("use_binette_checkm2_for_drep", False)
SKIP_CHECKM1 = RUN_BINETTE and USE_BINETTE_CHECKM2

# rule all must be the first rule in the file (Snakemake's implicit default
# target), so it comes before the helper functions/rules it transitively
# depends on.
rule all:
    input:
        *([] if SKIP_CHECKM1 else [f"{CHECKM_DIR}/checkm.txt"]),
        f"{OUT}/gtdbtk_classify",
        f"{DREP_DIR}/dRep.genomeInfo",
        f"{DREP_DIR}/dereplicated_genomes",
        *[f"{OUT}/per_group_drep/{g}/dRep/dereplicated_genomes"
          for g in config.get("per_group_drep_comparisons", {})],
        *[f"{OUT}/per_group_drep/{g}/gtdbtk_classify"
          for g in config.get("per_group_drep_comparisons", {})]

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

def get_binette_bin_dirs(wildcards):
    """Return list of binner output dirs for binette --bin_dirs, discovered
    directly on disk. Unlike metagenome_binning.smk (where Binette used to
    live), this workflow runs in a separate Snakemake invocation *after*
    binning has already finished, so there is no checkpoint to query here —
    binning's output directories are just read as plain paths.

    IMPORTANT: every directory here (default three, extra_binner_dirs,
    extra_binner_dirs_by_group) must contain bins whose contigs come from
    THIS group's own assembly (the same file passed as -c below) — Binette
    requires every input bin's contigs to exist in that one contigs FASTA,
    and errors out otherwise. This only ever holds for another binner run
    on the *same* assembly, never for bins from a genuinely different
    assembly (e.g. a separate co-assembly vs single-sample run of the same
    participants) — for that, see extra_bin_source_dirs below instead,
    which pools already-complete genome bins downstream of Binette, where
    whole-genome comparison (dRep's ANI) is the correct tool, not Binette."""
    candidates = [
        f"{OUT}/metabat2/{wildcards.assembly_group}",
        f"{OUT}/semibin/{wildcards.assembly_group}/output_bins",
        f"{OUT}/vamb/{wildcards.assembly_group}/bins",
    ] + [
        tpl.format(output_dir=OUT, assembly_group=wildcards.assembly_group)
        for tpl in EXTRA_BINNER_DIRS
    ] + EXTRA_BINNER_DIRS_BY_GROUP.get(wildcards.assembly_group, [])
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

def _own_bin_source_dir():
    """This run's own bins, before any cross-run pooling: Binette's refined
    bins when enabled, else the raw per-binner bins from metagenome_binning.smk."""
    return BINETTE_RENAMED_DIR if RUN_BINETTE else f"{OUT}/renamed_bins"

# Optional: compare bins across separate pipeline runs at the CheckM/dRep/
# GTDB-Tk level (not Binette — see the note on get_binette_bin_dirs above
# for why). Each entry is another run's bin_source_dir()-equivalent — its
# own binette_renamed_bins/ (if it used Binette) or renamed_bins/ (if not)
# — pooled in wholesale, no group/sample matching needed here: dRep's
# whole-genome ANI dereplication naturally collapses redundant genomes
# regardless of which run/assembly/participant they came from, and safely
# ignores anything that isn't actually similar.
#   extra_bin_source_dirs:
#     - "/abs/path/results_singlesample/renamed_bins"
#     - "/abs/path/results_singlesample/binette_renamed_bins"
EXTRA_BIN_SOURCE_DIRS = config.get("extra_bin_source_dirs", [])

def bin_source_dir():
    """Bins feeding CheckM/dRep/GTDB-Tk. Pooled (via pool_extra_bins) when
    extra_bin_source_dirs is set, otherwise this run's own bins directly."""
    return f"{OUT}/pooled_bins" if EXTRA_BIN_SOURCE_DIRS else _own_bin_source_dir()

if EXTRA_BIN_SOURCE_DIRS:
    rule pool_extra_bins:
        input:
            f"{_own_bin_source_dir()}/.bins_aggregated.done"
        output:
            f"{OUT}/pooled_bins/.pooled.done"
        run:
            os.makedirs(f"{OUT}/pooled_bins", exist_ok=True)
            for d in [_own_bin_source_dir()] + EXTRA_BIN_SOURCE_DIRS:
                for f in glob.glob(os.path.join(d, "*.fa.gz")) + glob.glob(os.path.join(d, "*.fasta.gz")):
                    dest = os.path.join(f"{OUT}/pooled_bins", os.path.basename(f))
                    if not os.path.exists(dest):
                        os.symlink(os.path.abspath(f), dest)
            with open(output[0], "w"):
                pass

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
    if EXTRA_BIN_SOURCE_DIRS:
        return [f"{OUT}/pooled_bins/.pooled.done"]
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

if SKIP_CHECKM1:
    rule generate_drep_info:
        """Binette already scored every winning bin with CheckM2
        (final_bins_quality_reports.tsv, per group) -- reuse that directly
        instead of running CheckM1 a second time on the same bin set."""
        input:
            f"{BINETTE_RENAMED_DIR}/.bins_aggregated.done"
        output:
            genome_info = f"{DREP_DIR}/dRep.genomeInfo"
        run:
            import pandas as pd
            out_rows = []
            for group in binnable_groups:
                report = f"{BINETTE_DIR}/{group}/final_bins_quality_reports.tsv"
                df_bt = pd.read_csv(report, sep="\t")
                df_bt.columns = [c.strip() for c in df_bt.columns]
                for _, row in df_bt.iterrows():
                    out_rows.append({
                        "genome": f"{group}_{row['name']}.fa.gz",
                        "completeness": row["completeness"],
                        "contamination": row["contamination"]
                    })
            pd.DataFrame(out_rows).to_csv(output.genome_info, index=False)
else:
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
    resources:
        tmpdir = "/scratch/tmp"
    shell:
        """
        # dRep's ANIn/ANImf secondary algorithms call nucmer directly, which
        # cannot read gzip-compressed FASTA -- it silently produces empty
        # alignments (ani=0.0 for every pair, including a genome against
        # itself) rather than erroring, so every genome ends up as its own
        # singleton cluster and dereplication looks like it ran but never
        # actually removes any redundancy. Mash-based algorithms (fastANI,
        # and primary clustering itself) read .fa.gz natively and are
        # unaffected -- only decompress when the chosen secondary algorithm
        # actually needs it.
        if [ "{params.alg}" = "ANIn" ] || [ "{params.alg}" = "ANImf" ]; then
            GENOME_DIR=$(mktemp -d {resources.tmpdir}/drep_decompressed.XXXXXX)
            trap "rm -rf $GENOME_DIR" EXIT
            for f in {params.bin_dir}/*.fa.gz; do
                gunzip -c "$f" > "$GENOME_DIR/$(basename "$f" .gz)"
            done
            # genomeInfo's "genome" column still says "*.fa.gz" -- dRep
            # matches it against the -g basenames verbatim, so a mismatch
            # here makes it think quality info is missing for every genome
            # and fall back to running CheckM itself (which isn't even on
            # PATH in this env, so it hard-errors). Strip .gz to match the
            # decompressed .fa names actually passed to -g.
            GENOME_INFO=$GENOME_DIR/dRep.genomeInfo.fa_names.csv
            sed 's/\.fa\.gz,/.fa,/' {input.genome_info} > $GENOME_INFO
        else
            GENOME_DIR={params.bin_dir}
            GENOME_EXT=".fa.gz"
            GENOME_INFO={input.genome_info}
        fi
        dRep dereplicate {params.drep_out} -g $GENOME_DIR/*${{GENOME_EXT:-.fa}} --genomeInfo $GENOME_INFO -p {threads} -comp {params.comp} -con {params.cont} -pa {params.pa} -sa {params.sa} -nc {params.nc} --S_algorithm {params.alg}
        """

############################################
# Per-participant dRep comparison (optional, alongside the global one above)
############################################
# The global checkm/drep/gtdbtk above pool every bin in this run (plus
# extra_bin_source_dirs, if set) into ONE dereplication. This section
# instead runs a SEPARATE, scoped CheckM+dRep+GTDB-Tk pass per configured
# group — e.g. "compare just spa_co1's own bins against just the matching
# single-sample run's bins for the same participant, with its own
# taxonomy" — without mixing in every other participant. This multiplies
# GTDB-Tk's (already the slowest step) runtime by however many groups are
# configured here — the global gtdbtk rule below still runs once more on
# top of that, on the full dereplicated dataset.
#
#   per_group_drep_comparisons:
#     spa_co1:                                       # a group in *this* run
#       output_dir:     "/abs/path/results_singlesample"
#       assemblies_tsv: "/abs/path/results_singlesample/assemblies_singlesample.tsv"
#
# Matching works the same way as extra_bin_source_dirs' sibling concept:
# spa_co1's constituent samples (from *this* run's mappings_tsv) are
# looked up by fq1 path in the other run's assemblies_tsv; whichever
# assembly_group(s) matched there are pooled in for this comparison —
# group/sample naming can differ freely between the two runs.
PER_GROUP_DREP_COMPARISONS = config.get("per_group_drep_comparisons", {})
_per_group_other_asm_df_cache = {}

def _other_run_groups_for(group):
    cfg = PER_GROUP_DREP_COMPARISONS.get(group)
    if not cfg:
        return []
    other_asm_tsv = cfg["assemblies_tsv"]
    if other_asm_tsv not in _per_group_other_asm_df_cache:
        df = pd.read_csv(other_asm_tsv, sep="\t").fillna("None")
        df.columns = [c.strip() for c in df.columns]
        _per_group_other_asm_df_cache[other_asm_tsv] = df
    other_asm_df = _per_group_other_asm_df_cache[other_asm_tsv]
    matched = set()
    for sample in ASM_TO_MAPPED_SAMPLES.get(group, []):
        r1, _ = MAPPING_READS[sample]
        matched |= set(other_asm_df.loc[other_asm_df["fq1"] == r1, "assembly_group"])
    return sorted(matched)

if PER_GROUP_DREP_COMPARISONS:
    _per_group_constraint = "|".join(re.escape(g) for g in PER_GROUP_DREP_COMPARISONS)

    def _link_prefixed(src_dir, prefix, dest_dir):
        if not os.path.isdir(src_dir):
            return
        for f in glob.glob(os.path.join(src_dir, f"{prefix}_*.fa.gz")) + glob.glob(os.path.join(src_dir, f"{prefix}_*.fasta.gz")):
            dest = os.path.join(dest_dir, os.path.basename(f))
            if not os.path.exists(dest):
                os.symlink(os.path.abspath(f), dest)

    def _per_group_pool_inputs(wildcards):
        inputs = [f"{_own_bin_source_dir()}/.bins_aggregated.done"]
        other_out = PER_GROUP_DREP_COMPARISONS[wildcards.compare_group]["output_dir"]
        for cand in (f"{other_out}/binette_renamed_bins/.bins_aggregated.done", f"{other_out}/renamed_bins/.bins_aggregated.done"):
            if os.path.exists(cand):
                inputs.append(cand)
                break
        return inputs

    rule pool_per_group_bins:
        input:
            _per_group_pool_inputs
        output:
            f"{OUT}/per_group_drep/{{compare_group}}/pooled_bins/.pooled.done"
        wildcard_constraints:
            compare_group = _per_group_constraint
        run:
            pool_dir = f"{OUT}/per_group_drep/{wildcards.compare_group}/pooled_bins"
            os.makedirs(pool_dir, exist_ok=True)
            _link_prefixed(_own_bin_source_dir(), wildcards.compare_group, pool_dir)
            other_out = PER_GROUP_DREP_COMPARISONS[wildcards.compare_group]["output_dir"]
            for other_group in _other_run_groups_for(wildcards.compare_group):
                # Prefer Binette's refined bins for this group; only fall back to the
                # raw per-binner bins if it has none (mirrors _own_bin_source_dir()).
                for cand_dir in (f"{other_out}/binette_renamed_bins", f"{other_out}/renamed_bins"):
                    matches = glob.glob(os.path.join(cand_dir, f"{other_group}_*.fa.gz")) + \
                              glob.glob(os.path.join(cand_dir, f"{other_group}_*.fasta.gz"))
                    if matches:
                        _link_prefixed(cand_dir, other_group, pool_dir)
                        break
            with open(output[0], "w"):
                pass

    rule checkm_per_group:
        input:
            f"{OUT}/per_group_drep/{{compare_group}}/pooled_bins/.pooled.done"
        output:
            report = f"{OUT}/per_group_drep/{{compare_group}}/checkm1/checkm.txt",
            out_dir = directory(f"{OUT}/per_group_drep/{{compare_group}}/checkm1")
        wildcard_constraints:
            compare_group = _per_group_constraint
        conda:
            "envs/checkm.yaml"
        params:
            db_path = config["checkm_db_path"],
            bin_dir = lambda w: f"{OUT}/per_group_drep/{w.compare_group}/pooled_bins",
            ext = config.get("checkm_extension", "gz")
        threads: 16
        shell:
            """
            export CHECKM_DATA_PATH={params.db_path}
            mkdir -p {output.out_dir}
            checkm lineage_wf -t {threads} --pplacer_threads {threads} -x {params.ext} --tab_table -f {output.report} {params.bin_dir} {output.out_dir}
            """

    rule generate_drep_info_per_group:
        input:
            checkm_report = f"{OUT}/per_group_drep/{{compare_group}}/checkm1/checkm.txt"
        output:
            genome_info = f"{OUT}/per_group_drep/{{compare_group}}/dRep/dRep.genomeInfo"
        wildcard_constraints:
            compare_group = _per_group_constraint
        run:
            df_cm = pd.read_csv(input.checkm_report, sep="\t")
            df_cm.columns = [c.strip() for c in df_cm.columns]
            bin_col = [c for c in df_cm.columns if "bin" in c.lower() or "id" in c.lower()][0]
            comp_col = [c for c in df_cm.columns if "complete" in c.lower()][0]
            cont_col = [c for c in df_cm.columns if "contam" in c.lower()][0]
            out_rows = []
            for _, row in df_cm.iterrows():
                bin_id = str(row[bin_col])
                if not bin_id.endswith(".fa.gz"):
                    bin_id = f"{bin_id}.fa.gz"
                out_rows.append({"genome": bin_id, "completeness": row[comp_col], "contamination": row[cont_col]})
            pd.DataFrame(out_rows).to_csv(output.genome_info, index=False)

    rule drep_per_group:
        input:
            pooled = f"{OUT}/per_group_drep/{{compare_group}}/pooled_bins/.pooled.done",
            genome_info = f"{OUT}/per_group_drep/{{compare_group}}/dRep/dRep.genomeInfo"
        output:
            out_dir = directory(f"{OUT}/per_group_drep/{{compare_group}}/dRep/dereplicated_genomes")
        wildcard_constraints:
            compare_group = _per_group_constraint
        conda:
            "envs/drep.yaml"
        params:
            bin_dir = lambda w: f"{OUT}/per_group_drep/{w.compare_group}/pooled_bins",
            drep_out = lambda w: f"{OUT}/per_group_drep/{w.compare_group}/dRep",
            comp    = config.get("drep_completeness", 66),
            cont    = config.get("drep_contamination", 20),
            pa      = config.get("drep_primary_ani", 0.95),
            sa      = config.get("drep_secondary_ani", 0.98),
            nc      = config.get("drep_min_overlap", 0.30),
            alg     = config.get("drep_s_algorithm", "ANIn")
        threads: 8
        resources:
            tmpdir = "/scratch/tmp"
        shell:
            """
            # See rule drep's identical guard above: ANIn/ANImf call nucmer
            # directly, which silently fails on gzip-compressed FASTA.
            if [ "{params.alg}" = "ANIn" ] || [ "{params.alg}" = "ANImf" ]; then
                GENOME_DIR=$(mktemp -d {resources.tmpdir}/drep_decompressed.XXXXXX)
                trap "rm -rf $GENOME_DIR" EXIT
                for f in {params.bin_dir}/*.fa.gz; do
                    gunzip -c "$f" > "$GENOME_DIR/$(basename "$f" .gz)"
                done
                # Also see rule drep: genomeInfo's "genome" column must match
                # the decompressed .fa basenames, not the original .fa.gz.
                GENOME_INFO=$GENOME_DIR/dRep.genomeInfo.fa_names.csv
                sed 's/\.fa\.gz,/.fa,/' {input.genome_info} > $GENOME_INFO
            else
                GENOME_DIR={params.bin_dir}
                GENOME_EXT=".fa.gz"
                GENOME_INFO={input.genome_info}
            fi
            dRep dereplicate {params.drep_out} -g $GENOME_DIR/*${{GENOME_EXT:-.fa}} --genomeInfo $GENOME_INFO -p {threads} -comp {params.comp} -con {params.cont} -pa {params.pa} -sa {params.sa} -nc {params.nc} --S_algorithm {params.alg}
            """

    rule gtdbtk_per_group:
        input:
            f"{OUT}/per_group_drep/{{compare_group}}/pooled_bins/.pooled.done"
        output:
            out_dir = directory(f"{OUT}/per_group_drep/{{compare_group}}/gtdbtk_classify")
        wildcard_constraints:
            compare_group = _per_group_constraint
        conda:
            "envs/gtdbtk.yaml"
        params:
            db_path = config["gtdbtk_db_path"],
            bin_dir = lambda w: f"{OUT}/per_group_drep/{w.compare_group}/pooled_bins",
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