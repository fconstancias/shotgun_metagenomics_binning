# Include shared data structures and imports
include: "common.smk"

# Symlink creation runs on the head node — no cluster node needed
localrules: use_prebuilt_assembly

# Target only assembly outputs — path encodes which assembler was used
rule all:
    input:
        [f"{OUT}/assembly/{ASSEMBLER}/{g}/final.contigs.fa" for g in all_groups]

rule use_prebuilt_assembly:
    input:
        external = lambda w: PREBUILT_ASM[w.assembly_group]
    output:
        contigs = f"{OUT}/assembly/{ASSEMBLER}/{{assembly_group}}/final.contigs.fa"
    wildcard_constraints:
        assembly_group = prebuilt_constraint
    run:
        os.makedirs(os.path.dirname(output.contigs), exist_ok=True)
        if os.path.exists(output.contigs):
            os.remove(output.contigs)
        os.symlink(os.path.abspath(input.external), output.contigs)

rule spades_assemble:
    input:
        r1 = lambda w: ASM_READS_R1[w.assembly_group],
        r2 = lambda w: ASM_READS_R2[w.assembly_group]
    output:
        contigs = f"{OUT}/assembly/spades/{{assembly_group}}/final.contigs.fa"
    log:
        err = f"logs/assembly/spades/{{assembly_group}}.log"
    conda:
        "envs/assembly.yaml"
    wildcard_constraints:
        assembly_group = built_constraint
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 1440
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
        err = f"logs/assembly/megahit/{{assembly_group}}.log"
    conda:
        "envs/assembly.yaml"
    wildcard_constraints:
        assembly_group = built_constraint
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 720
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