include: "common.smk"

CHECKM_DIRNAME = config.get("checkm_outdir_name", "checkm1")
DREP_DIRNAME = config.get("drep_outdir_name", "dRep")
CHECKM_DIR = f"{OUT}/{CHECKM_DIRNAME}"
DREP_DIR = f"{OUT}/{DREP_DIRNAME}"

def get_all_renamed_bins_local(wildcards):
    """
    Finds preexisting bins inside the renamed directory.
    This enables running this evaluation module directly on precalculated bins.
    """
    bin_dir = f"{OUT}/renamed_bins"
    if os.path.exists(bin_dir):
        bins = glob.glob(os.path.join(bin_dir, "*.fa.gz")) + glob.glob(os.path.join(bin_dir, "*.fasta.gz"))
        if bins:
            return bins
    # Fallback to trigger upstream logic if empty
    return [f"{OUT}/renamed_bins/.bins_aggregated.done"]

rule all:
    input:
        f"{CHECKM_DIR}/checkm.txt",
        f"{OUT}/gtdbtk_classify",
        f"{DREP_DIR}/dRep.genomeInfo",
        f"{DREP_DIR}/dereplicated_genomes"

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
        bin_dir = f"{OUT}/renamed_bins",
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
        bin_dir = f"{OUT}/renamed_bins",
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
        bin_dir = f"{OUT}/renamed_bins",
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