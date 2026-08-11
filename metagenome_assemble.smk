# Include shared data structures and imports
include: "common.smk"

# Symlink creation runs on the head node — no cluster node needed
localrules: use_prebuilt_assembly

def get_assembly_fasta(wildcards):
    """Return reformatted contigs if anvi_reformat else the raw assembly —
    same convention metagenome_binning.smk uses to select its own input."""
    assembler = ASSEMBLER_FOR[wildcards.assembly_group]
    base = f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.fa"
    if ANVI_REFORMAT:
        return f"{OUT}/assembly/{assembler}/{wildcards.assembly_group}/final.contigs.reformatted.fa"
    return base

# Contigs DB creation (anvi-gen-contigs-database) lives here, not in
# metagenome_binning.smk, so it's shared between the CONCOCT track (binning)
# and the optional SCG-taxonomy/contig-stats step below without being built
# twice for the same assembly group.
_anvi_db_targets = []
if RUN_CONCOCT:
    _anvi_db_targets += [f"{OUT}/concoct/{g}/{g}.db" for g in concoct_capable_groups]
if ANVI_REFORMAT and ANVI_TAXONOMY_AND_STATS:
    _anvi_db_targets += [f"{OUT}/concoct/{g}/.taxonomy.done" for g in all_groups]
    _anvi_db_targets += [f"{OUT}/concoct/{g}/{g}_contigs_stats.txt" for g in all_groups]

# Target only assembly outputs — path encodes which assembler was used.
# Target the reformatted contigs directly when anvi_reformat is on, so
# metagenome_binning.smk finds them already prepared either way (it does no
# reformatting itself — that's purely a mapping/binning workflow).
rule all:
    input:
        [f"{OUT}/assembly/{ASSEMBLER_FOR[g]}/{g}/final.contigs.fa" for g in all_groups] + (
            [f"{OUT}/assembly/{ASSEMBLER_FOR[g]}/{g}/final.contigs.reformatted.fa" for g in all_groups]
            if ANVI_REFORMAT else []
        ) + _anvi_db_targets

rule use_prebuilt_assembly:
    input:
        external = lambda w: PREBUILT_ASM[w.assembly_group]
    output:
        contigs = f"{OUT}/assembly/{{assembler}}/{{assembly_group}}/final.contigs.fa"
    wildcard_constraints:
        assembly_group = prebuilt_constraint
    run:
        os.makedirs(os.path.dirname(output.contigs), exist_ok=True)
        if os.path.exists(output.contigs):
            os.remove(output.contigs)
        os.symlink(os.path.abspath(input.external), output.contigs)

rule reformat_contigs:
    """Normalizes headers/filters short contigs via anvi-script-reformat-fasta,
    for both freshly-assembled and pre-built contigs — this is a one-time
    step on the assembly output, not something metagenome_binning.smk (map +
    bin only) needs to concern itself with."""
    input:
        f"{OUT}/assembly/{{assembler}}/{{assembly_group}}/final.contigs.fa"
    output:
        f"{OUT}/assembly/{{assembler}}/{{assembly_group}}/final.contigs.reformatted.fa"
    params:
        min_len = config.get("anvi_min_contig_len", 1000),
        prefix  = lambda w: w.assembly_group
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    shell:
        """
        anvi-script-reformat-fasta {input} \
            --simplify-names \
            --min-len {params.min_len} \
            --prefix {params.prefix} \
            -o {output}
        """

############################################
# Anvi'o contigs DB (CONCOCT track + optional taxonomy/stats below)
############################################

rule anvi_gen_contigs_db:
    input:
        assembly = get_assembly_fasta
    output:
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db"
    params:
        name = lambda w: w.assembly_group
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        """
        mkdir -p {OUT}/concoct/{wildcards.assembly_group}
        anvi-gen-contigs-database -T {threads} \
            -f {input.assembly} \
            -o {output.db} \
            -n {params.name}
        """

############################################
# Optional: SCG taxonomy + contig stats (requires anvi_reformat: true and
# anvi_taxonomy_and_stats: true — assumes `anvi-setup-scg-taxonomy` has
# already been run once for the anvio conda env)
############################################

rule anvi_run_hmms:
    """anvi-run-scg-taxonomy requires SCG HMM hits (e.g. Bacteria_71) to
    already be annotated in the contigs DB."""
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
        anvi-run-hmms -c {input.db} -T {threads} --just-do-it
        touch {output}
        """

rule anvi_run_scg_taxonomy:
    input:
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db",
        hmms_done = f"{OUT}/concoct/{{assembly_group}}/.hmms.done"
    output:
        f"{OUT}/concoct/{{assembly_group}}/.taxonomy.done"
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 120
    shell:
        """
        anvi-run-scg-taxonomy -c {input.db} -T {threads}
        touch {output}
        """

rule anvi_contigs_stats:
    input:
        db = f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}.db"
    output:
        f"{OUT}/concoct/{{assembly_group}}/{{assembly_group}}_contigs_stats.txt"
    conda:
        conda_env("anvio", "envs/anvio.yaml")
    shell:
        "anvi-display-contigs-stats {input.db} --report-as-text -o {output}"

rule spades_assemble:
    input:
        r1 = lambda w: ASM_READS_R1[w.assembly_group],
        r2 = lambda w: ASM_READS_R2[w.assembly_group]
    output:
        contigs = f"{OUT}/assembly/spades/{{assembly_group}}/final.contigs.fa"
    log:
        err = f"{OUT}/logs/assembly/spades/{{assembly_group}}.log"
    conda:
        "envs/assembly.yaml"
    wildcard_constraints:
        assembly_group = built_constraint
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 1440,
        assembly_slots = 1
    params:
        mem_gb = lambda wildcards, resources: max(4, int(resources.mem_mb / 1024) - 3),
        spades_input = lambda wildcards, input: " ".join(
            f"--pe{idx}-1 {r1} --pe{idx}-2 {r2}"
            for idx, (r1, r2) in enumerate(zip(input.r1, input.r2), start=1)
        ),
        tmpdir = f"{OUT}/assembly/spades/{{assembly_group}}/tmp_spades"
    shell:
        """
        rm -rf {params.tmpdir}
        spades.py --meta {params.spades_input} -t {threads} -m {params.mem_gb} \
            -o {params.tmpdir} > {log.err} 2>&1
        mv {params.tmpdir}/scaffolds.fasta {output.contigs}
        rm -rf {params.tmpdir}
        """

rule megahit_assemble:
    input:
        r1 = lambda w: ASM_READS_R1[w.assembly_group],
        r2 = lambda w: ASM_READS_R2[w.assembly_group]
    output:
        contigs = f"{OUT}/assembly/megahit/{{assembly_group}}/final.contigs.fa"
    log:
        err = f"{OUT}/logs/assembly/megahit/{{assembly_group}}.log"
    conda:
        "envs/assembly.yaml"
    wildcard_constraints:
        assembly_group = built_constraint
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 720,
        assembly_slots = 1
    params:
        r1_comma = lambda wildcards, input: ",".join(input.r1),
        r2_comma = lambda wildcards, input: ",".join(input.r2),
        tmpdir = f"{OUT}/assembly/megahit/{{assembly_group}}/tmp_megahit"
    shell:
        """
        rm -rf {params.tmpdir}
        megahit -1 {params.r1_comma} -2 {params.r2_comma} --min-contig-len 500 -m 0.99 \
            -t {threads} --out-dir {params.tmpdir} --out-prefix {wildcards.assembly_group} > {log.err} 2>&1
        mv {params.tmpdir}/{wildcards.assembly_group}.contigs.fa {output.contigs}
        rm -rf {params.tmpdir}
        """