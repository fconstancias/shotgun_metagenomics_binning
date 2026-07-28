# Test C — full pipeline from raw reads (mixed assemblers + CONCOCT)

This test exercises the complete 3-workflow pipeline end-to-end, starting
from raw reads rather than pre-built assemblies:

```
reads ──▶ [1] metagenome_assemble.smk ──▶ [2] metagenome_binning.smk ──▶ [3] summarise_mags.smk
          (real assembly, mixed          (mapping + metabat2/semibin2/     (Binette + CheckM +
           assemblers, Anvi'o reformat    vamb + CONCOCT)                   GTDB-Tk + dRep)
           + SCG taxonomy/stats)
```

Files:
- `assemblies_toy_realasm.tsv` — assembly groups + per-group assembler
- `mappings_toy.tsv` — shared with Test B; already matches the sample/tool
  layout this test needs, so it's reused as-is
- `config_assemble_toy_realasm.yaml`
- `config_binning_toy_realasm.yaml`
- `config_summarise_toy_realasm.yaml`

Output lands in `results_toy_realasm/`.

## What makes this test different from Test B (`*_toy_full.yaml`)

Test B symlinks pre-built SPAdes scaffolds (`assemblies_toy.tsv`). This test
(`assemblies_toy_realasm.tsv`) assembles from raw reads instead, and uses a
**different assembler per assembly group** via a per-row `assembler` column
in the TSV (falls back to the config's global `assembler:` key for any
TSV/row that doesn't set it) — a new capability added specifically for this
test (see `ASSEMBLER_FOR` in `common.smk`).

## Assembly group layout

| assembly_group | samples (mapping_tool)                                              | assembler | binners (via config)                  |
|-----------------|----------------------------------------------------------------------|-----------|-----------------------------------------|
| `spaS144`       | S_144, S_178, S_194, S_9, S_265 — **all bowtie2**                    | megahit   | metabat2, semibin2, vamb, **CONCOCT**   |
| `spaS276`       | S_276 (bowtie2, assembled) + S_133, S_86, S_265 (strobealign, cross) | megahit   | metabat2, vamb                          |
| `spaS135`       | S_135 (bowtie2, assembled) + S_119, S_213, S_280 (strobealign, cross)| spades    | metabat2, vamb                          |

- `spaS144` is a real **co-assembly**: all 5 samples' reads are assembled
  together (every sample is bowtie2-mapped in `mappings_toy.tsv`), which is
  what makes it both CONCOCT-capable (needs ≥2 bowtie2 profiles) and
  SemiBin2-capable (needs every sample bowtie2, `-b` direct-BAM mode).
- `spaS276` / `spaS135` are **single-sample assemblies**: only the
  bowtie2-mapped "self" sample's reads go into the assembly; the other
  samples in `mappings_toy.tsv` are strobealign cross-mapped for depth
  signal only, never assembled. SemiBin2 and CONCOCT are skipped for these
  groups (mixed bowtie2+strobealign / single bowtie2 profile respectively) —
  they still get MetaBAT2 + VAMB.

## Config highlights

- `config_assemble_toy_realasm.yaml`: `anvi_reformat: true`,
  `anvi_taxonomy_and_stats: true` (anvi-run-scg-taxonomy +
  anvi-display-contigs-stats — requires `anvi-setup-scg-taxonomy` to have
  been run once for the `anvio` conda env; already done here), `run_concoct:
  true` (builds the Anvi'o contigs DB shared with CONCOCT in workflow 2).
- `config_binning_toy_realasm.yaml`: `binners: [metabat2, semibin2, vamb]` +
  `run_concoct: true`; `semibin2_epochs: 2` (reduced from SemiBin2's default
  15 training epochs — toy/code-path test only, not meant to produce
  meaningful bins).
- `config_summarise_toy_realasm.yaml`: `run_binette: true` — refines the
  per-binner bins (MetaBAT2/SemiBin2/VAMB/CONCOCT-adjacent multi-binner
  track only; CONCOCT itself is a separate manual-curation track, not fed
  into Binette) before CheckM/GTDB-Tk/dRep.

## Running

From this directory, with the `snakemake` conda env active:

```bash
export CONDA_PREFIX=/home/ljc444/.conda/snakemake_envs
export SNAKEMAKE_FLAGS="--use-conda --conda-frontend conda --conda-prefix $CONDA_PREFIX"

# 1) Assembly
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS

# 2) Mapping + binning (after 1 completes)
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS

# 3) Summarise (after 2 completes)
snakemake -s ../summarise_mags.smk --configfile config_summarise_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS
```

Dry-run any of the above first with `-n` to sanity-check the DAG.
