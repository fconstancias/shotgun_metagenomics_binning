# shotgun_metagenomics_binning

A modular Snakemake pipeline for metagenomic assembly, read mapping, binning, and MAG dereplication/taxonomy. The three workflows can be run independently or chained sequentially.

---
## TODO:

- Branch with CoverM for mapping with option to use galah instead of Drep in summarise_mags
- later: use Simka / SimkaMin to select samples for mapping
- later: ONT data ?

## Overview

```
Host-removed reads
   │
   ▼
[1] metagenome_assemble.smk   ─── Assembly (SPAdes or MEGAHIT) or symlink to
   │                               pre-built contigs → optional Anvi'o reformat
   │                               → optional contigs DB + SCG taxonomy/stats
   ▼
[2] metagenome_binning.smk    ─── Read mapping (bowtie2) → strobealign aemb
   │                               coverage → MetaBAT2 / SemiBin2 / VAMB
   │                               binning → optional CONCOCT (via Anvi'o,
   │                               reusing [1]'s contigs DB)
   ▼
[3] summarise_mags.smk        ─── Optional Binette refinement, CheckM quality,
                                   GTDB-Tk taxonomy, dRep dereplication
```

---

## Setup

### Snakemake conda environment

The system conda on esrum (23.3.1) is too old for Snakemake 9+, which requires ≥24.7.1. Create a dedicated `snakemake` conda env with a modern conda bundled inside it, then verify:

```bash
conda create -n snakemake -c conda-forge conda">=24.7" mamba snakemake python=3.12
conda activate snakemake
export CONDA_EXE=~/.conda/envs/snakemake/bin/conda   # add to ~/.bashrc, or repeat each session
$CONDA_EXE --version   # should show 24.x
snakemake --version    # should show 9.x
```

Snakemake reads `CONDA_EXE` to call conda directly, bypassing the shell function that wraps the system conda.

### Conda env prefix for pipeline tool envs

All pipeline tool envs (mapping, binning, semibin2, vamb, binette, …) are installed once into a shared prefix and reused across runs. Every run command below assumes these two exports — pick any writable directory for `CONDA_PREFIX` (first run builds each env there; later runs reuse them):

```bash
export CONDA_PREFIX=/path/to/your/shared/conda_envs   # e.g. ~/.conda/snakemake_envs
export SNAKEMAKE_FLAGS="--use-conda --conda-frontend conda --conda-prefix $CONDA_PREFIX"
```

### SLURM executor plugin (required for `--profile cluster/`)

Snakemake 9+ moved cluster backends into separate executor plugins — install the SLURM one into the `snakemake` env:

```bash
mamba install -n snakemake -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

`cluster/config.yaml` also needs your SLURM account and partition set under `default-resources` (see [Troubleshooting](#troubleshooting) if you hit account/partition errors):

```yaml
default-resources:
  slurm_account:   cbmr           # find yours with: sacctmgr -n -s list user "$USER" format=account
  slurm_partition: standardqueue  # find available partitions with: sinfo -o "%P %l %D %c %m"
```

### External reference databases

Needed for `summarise_mags.smk` (CheckM/GTDB-Tk/Binette) and the optional
`anvi_taxonomy_and_stats` step — none of these are downloaded by the
pipeline itself, and none are small:

```bash
# CheckM v1 (checkm_db_path)
mkdir -p /path/to/checkm_db && cd /path/to/checkm_db
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzf checkm_data_2015_01_16.tar.gz

# GTDB-Tk (gtdbtk_db_path) — the download-db.sh script ships with the
# gtdbtk conda package; run it from inside that env. Large (100GB+).
conda activate <your gtdbtk env>
download-db.sh

# CheckM2 (binette_checkm2_db, used internally by Binette)
conda activate <your binette/checkm2 env>
checkm2 database --download --path /path/to/checkm2_db

# Anvi'o SCG taxonomy (only if anvi_taxonomy_and_stats: true)
conda activate <your anvio env>
anvi-setup-scg-taxonomy
```

Point `checkm_db_path`, `gtdbtk_db_path`, and `binette_checkm2_db` in your
config at wherever you put these.

---

## Input TSV Files

### `assemblies.tsv` — used by workflows 1 and 2

Defines assembly groups and their read inputs. Each row is one sample contributing to an assembly group. Leave `assembly_path` blank (or `None`) to trigger de-novo assembly; provide a path to skip assembly and symlink existing contigs.

| Column           | Required | Description |
|------------------|----------|-------------|
| `assembly_group` | ✅*      | Group name (e.g. `co_assembly_A`, `sample01`). Multiple rows with the same group = co-assembly. Can be left blank when `group_key` (below) is used to auto-generate it. |
| `assembly_path`  | ✅       | Path to pre-built contigs (`.fa`/`.fasta`). Leave empty to assemble from reads. |
| `fq1`            | ✅*      | Path to forward reads (R1). Required if `assembly_path` is empty. |
| `fq2`            | ✅*      | Path to reverse reads (R2). Required if `assembly_path` is empty. |
| `assembler`      | ❌       | Optional per-row override of `assembler:` (megahit/spades) — lets different assembly groups use different assemblers in the same run. Falls back to the config's `assembler:` key if omitted. |
| `group_key`       | ❌       | Optional. See **auto-naming** below. |

Example:
```tsv
assembly_group	assembly_path	fq1	fq2
co_assembly_A		/data/sampleA_R1.fq.gz	/data/sampleA_R2.fq.gz
co_assembly_A		/data/sampleB_R1.fq.gz	/data/sampleB_R2.fq.gz
prebuilt_B	/data/contigs_B.fa	None	None
```

**Co-assembly and single-sample assembly** both work with the same TSV format:
- Co-assembly: multiple rows share the same `assembly_group`.
- Single-sample: one row per `assembly_group` (the group name is typically the sample name).

#### Auto-naming assembly groups (`group_key`)

Instead of inventing `assembly_group` names by hand, add a `group_key` column (e.g. a participant/cohort id) and leave `assembly_group` blank (or `None`/`auto`) on those rows. Rows sharing the same `group_key` are grouped exactly like rows sharing the same `assembly_group` today — the key just controls *grouping*, not the final name:
- **Co-assembly** (a `group_key` shared by ≥2 rows) → named `{assembler_abbrev}_co{n}`, e.g. `spa_co1`, `mh_co2`. `n` increments per assembler, in TSV row order.
- **Single-sample** (a `group_key` used by exactly 1 row) → named `{assembler_abbrev}_{group_key}`, e.g. `spa_S276`.

`assembler_abbrev` defaults to `{"spades": "spa", "megahit": "mh"}`; override/extend via the `assembler_abbrev:` config key (a dict mapping assembler name → short label). Any row with an explicit, non-blank `assembly_group` is left untouched regardless of `group_key` — auto-naming only fills in blanks. TSVs without a `group_key` column are completely unaffected (fully backward compatible).

---

### `mappings.tsv` — used by workflows 2 and 3

Defines which samples are mapped back to which assembly group for depth estimation.

| Column           | Required | Description |
|------------------|----------|-------------|
| `sample_id`      | ✅       | Unique sample identifier. |
| `assembly_group` | ✅       | Assembly group this sample maps to (must match `assemblies.tsv`). |
| `fq1`            | ✅       | Path to forward reads (R1). |
| `fq2`            | ✅       | Path to reverse reads (R2). |
| `mapping_tool`   | ✅       | `bowtie2` or `strobealign`. Resolved per `(sample_id, assembly_group)` row. |

`mapping_tool` is looked up per `(sample_id, assembly_group)` pair, so **the same sample can use a different tool depending on which assembly group it's mapped to**.

**Recommended policy**, driven by which binners you want available for a group (see [Coverage strategy](#coverage-strategy)):
- **Co-assembly groups → `bowtie2` for every sample.** Every sample gets a real BAM, which unlocks SemiBin2 (`-b`, direct BAM input) and CONCOCT (needs ≥2 bowtie2 profiles) in addition to MetaBAT2/VAMB.
- **Single-sample assembly groups → `bowtie2` for the assembled sample itself** (accurate self-depth), **`strobealign` for every other cross-mapped sample** (same participant's other timepoints, or any other samples used purely for cross-sample depth signal). MetaBAT2 and VAMB both handle this mix fine; SemiBin2 and CONCOCT are skipped for these groups (see below).

Example:
```tsv
sample_id	assembly_group	fq1	fq2	mapping_tool
sampleA	co_assembly_A	/data/sampleA_R1.fq.gz	/data/sampleA_R2.fq.gz	bowtie2
sampleB	co_assembly_A	/data/sampleB_R1.fq.gz	/data/sampleB_R2.fq.gz	bowtie2
sampleA	prebuilt_B	/data/sampleA_R1.fq.gz	/data/sampleA_R2.fq.gz	bowtie2
sampleC	prebuilt_B	/data/sampleC_R1.fq.gz	/data/sampleC_R2.fq.gz	strobealign
```

---

## Workflow Details

### 1. `metagenome_assemble.smk`

**Config:** `config_assemble.yaml`

For each assembly group:
- If `assembly_path` is provided → symlinks the existing contigs into `{output_dir}/assembly/{assembler}/{group}/final.contigs.fa`
- If not → runs **SPAdes** (`--meta`) or **MEGAHIT** depending on `assembler:` in config

Then, for **every** group (regardless of source), if `anvi_reformat: true` → reformats via `anvi-script-reformat-fasta` into `final.contigs.reformatted.fa`. This always happens here, never in `metagenome_binning.smk` — that workflow only maps and bins; `anvi_reformat` there just selects which of these two files to use as input, and must match the value used here.

Also builds the Anvi'o contigs DB (`{output_dir}/concoct/{group}/{group}.db`) whenever it's needed — either by CONCOCT (`run_concoct: true`, see below) or by the optional SCG-taxonomy/contig-stats step — so it's created exactly once per group and shared between both, rather than rebuilt separately.

Key config options:
```yaml
assembler: "megahit"   # or "spades"
output_dir: "results"
assemblies_tsv: "assemblies.tsv"
mappings_tsv:   "mappings.tsv"   # needed for concoct_capable_groups — must match config_binning.yaml

# Reformatting — must match the same keys in config_binning.yaml, which
# just selects reformatted-vs-raw contigs as input (see workflow 2 below)
anvi_reformat:       false
anvi_min_contig_len: 1000

# Must match config_binning.yaml's run_concoct — whether to build the
# contigs DB for CONCOCT-capable groups (>=2 bowtie2-mapped samples)
run_concoct: false

# Optional: anvi-run-scg-taxonomy + anvi-display-contigs-stats per group.
# Requires anvi_reformat: true, and anvi-setup-scg-taxonomy to have been
# run once for the anvio conda env.
anvi_taxonomy_and_stats: false

conda_dirs:
  anvio: null   # e.g. /home/ljc444/.conda/envs/anvio-9
```

Run:
```bash
snakemake -s metagenome_assemble.smk \
  --profile cluster/ --configfile config_assemble.yaml \
  $SNAKEMAKE_FLAGS --jobs 15 --rerun-incomplete --latency-wait 60
```

#### Optional: SCG taxonomy + contig stats (`anvi_taxonomy_and_stats: true`)

Runs on the Anvi'o contigs DB, per assembly group:
- **`anvi-run-scg-taxonomy`** — annotates single-copy gene taxonomy directly in the contigs DB
- **`anvi-display-contigs-stats --report-as-text`** — writes basic assembly stats (N50, contig count, etc.) to `{group}_contigs_stats.txt`

Requires `anvi_reformat: true` and `anvi-setup-scg-taxonomy` to have been run once for the `anvio` conda env (one-time setup, downloads GTDB SCG reference data — not managed by this pipeline).

---

### 2. `metagenome_binning.smk`

**Config:** `config_binning.yaml`

#### Coverage strategy

**MetaBAT2 depth** (driven by `mapping_tool` per sample)
- `bowtie2` samples: full alignment → BAM → `jgi_summarize_bam_contig_depths` → `bowtie2_depth.txt`
- `strobealign` samples: `--aemb` per-contig coverage → per-sample TSVs (bowtie2-mapped samples get an equivalent TSV derived from their own BAM instead of mapping a second time, via the same `aemb_tsv` rule — the two are numerically equivalent, within ~1-2%)
- Both merged into `depth.txt` by `combine_depths`

**VAMB** — uses the same per-sample aemb-style TSVs (real or BAM-derived), merged into one `--abundance_tsv` (VAMB 5.x has no `--aemb` flag; requires a single TSV with header `contigname\t<sample>...`). `--minfasta` is set so VAMB actually writes per-bin FASTA files (`vamb/{group}/bins/*.fna`), not just its cluster table.

**SemiBin2** — its `-b` (BAM) and `-a` (strobealign-aemb) inputs are mutually exclusive per run, and `-a` mode requires SemiBin2's own split-contig abundance format (see `SemiBin2 split_contigs`), not plain per-contig depth. Rather than adopt that extra split-contig step, SemiBin2 only runs for groups where **every** sample is bowtie2-mapped — using `-b` directly on the real BAMs, matching SemiBin2's own recommended BAM-based workflow. Mixed bowtie2+strobealign groups (typical single-sample assemblies) skip SemiBin2 automatically — see the `mapping_tool` policy above.

**CONCOCT** has the same ≥2-bowtie2-samples requirement (see below) — so co-assembly groups (all-bowtie2) are where MetaBAT2, VAMB, SemiBin2, and CONCOCT all work together; single-sample assembly groups get MetaBAT2 + VAMB only.

#### Steps

1. **(Optional) Anvi'o contig reformat** — `anvi-script-reformat-fasta --simplify-names --min-len X --prefix {group}`. Runs before `bowtie2_build`, so the index is built from the simplified-header FASTA and everything stays consistent. **Only applies to freshly-assembled (SPAdes/MEGAHIT) contigs** — pre-built assemblies are reformatted earlier, in `metagenome_assemble.smk`'s `reformat_prebuilt_contigs`.
2. **bowtie2-build** index (if any sample uses `bowtie2` for that group)
3. **bowtie2 map + sort** per bowtie2 sample → BAM
4. **bowtie2 depth** → `bowtie2_depth.txt`
5. **strobealign aemb** per `(group, sample)` → `aemb/{group}/{sample}.tsv`
6. **combine_depths** → `depth.txt` (merges bowtie2 depths + strobealign aemb TSVs)
7. **MetaBAT2** (if `metabat2` in `binners`) — uses `depth.txt`
8. **SemiBin2** (if `semibin2` in `binners`) — only for groups where every sample is bowtie2-mapped; uses `-b` directly on the BAMs
9. **VAMB** (if `vamb` in `binners`) — merged `--abundance_tsv` from the per-sample aemb TSVs, latest VAMB API (`vamb bin default --abundance_tsv ... --minfasta ...`)

Multi-binner refinement with Binette now runs later, in `summarise_mags.smk`
(see below) — it reads this workflow's raw per-binner bin directories
directly from `output_dir`, so it can be re-run/tuned without re-running
binning.

Key config options:
```yaml
assemblies_tsv: "assemblies.tsv"
mappings_tsv:   "mappings.tsv"
output_dir:     "results"
assembler:      "spades"   # must match workflow 1

# Binner selection: any combination of metabat2, semibin2, vamb
binners: [metabat2]

# Optional: override SemiBin2's default 15 training epochs (~1h/epoch even
# on small data) — e.g. for fast test runs. Leave unset for real runs.
# semibin2_epochs: 2

# Minimum bin size (bp) for VAMB to write it out as a FASTA file (into
# vamb/{group}/bins/); bins below this are only listed in the cluster table.
vamb_min_fasta_size: 200000

# Selects reformatted vs raw contigs as input (reformatting itself always
# happens in metagenome_assemble.smk) — must match that run's anvi_reformat
anvi_reformat: false

# Optional: pre-installed conda env paths (skip yaml-based install)
# Leave null to install from envs/*.yaml
conda_dirs:
  mapping:  null   # e.g. /home/ljc444/.conda/envs/mapping_env
  binning:  null   # e.g. /home/ljc444/.conda/envs/metabat2
  semibin2: null   # e.g. /home/ljc444/.conda/envs/semibin
  vamb:     null
  anvio:    null   # e.g. /home/ljc444/.conda/envs/anvio-9

# MetaBAT2 tuning (all optional, defaults shown)
metabat2:
  min_contig:   2000
  max_edges:    500
  min_cv:       1.0
  min_cv_sum:   1.0
  max_p:        95
  min_s:        60
  min_cls_size: 200000
  save_cls:     true
  verbose:      true
```

Run:
```bash
snakemake -s metagenome_binning.smk \
  --profile cluster/ --configfile config_binning.yaml \
  $SNAKEMAKE_FLAGS --jobs 15 --rerun-incomplete --latency-wait 60
```

#### Pre-building conda environments (run once before cluster submission)

```bash
snakemake -s metagenome_binning.smk --configfile config_binning.yaml \
  $SNAKEMAKE_FLAGS --conda-create-envs-only --cores 1
```

To skip installation for tools already available, set the relevant path under `conda_dirs:` in the config.

#### CONCOCT binning via Anvi'o (optional, parallel track)

CONCOCT is an alternative binning algorithm that uses Anvi'o for coverage profiling and clustering. It runs **in parallel** to the standard MetaBAT2/SemiBin2/VAMB track and produces independent output, allowing exploration of multiple binning solutions.

**Prerequisite:** CONCOCT is not bundled with Anvi'o and must be built into the `anvio-9` conda env once, per [Anvi'o's own instructions](https://anvio.org/install/ubuntu/#then-you-are-ready-to-install-concoct):
```bash
conda activate anvio-9
pip install cython "setuptools<81" nose   # build + runtime deps missing from anvio-9
mkdir -p ~/github && cd ~/github
git clone https://github.com/merenlab/CONCOCT.git
cd CONCOCT
python setup.py build
python setup.py install
```
Verify with `anvi-cluster-contigs --help` — the `CONCOCT` section should no longer say `[NOT FOUND]`.

**When to use CONCOCT:**
- Exploring alternative binning strategies before final curation
- Datasets where MetaBAT2/Binette results need validation
- Manual bin refinement workflows (each cluster count produces a separate bin set)

**Workflow:**
1. **Anvi'o contigs DB** — created in `metagenome_assemble.smk` (`run_concoct: true` there too), not here — see workflow 1 above
2. **HMM annotations** — adds functional markers to contigs DB
3. **Per-sample profiling** — runs anvi-profile on each bowtie2-mapped BAM
4. **Merge profiles** — combines all samples into one merged profile per assembly_group
5. **CONCOCT clustering** — clusters contigs at multiple thresholds (e.g., 10, 15, 20 clusters)
6. **Summarize & export** — generates bin FASTA files and HTML reports per cluster count

**Config options:**
```yaml
run_concoct: false                    # set to true to enable
concoct_clusters: [10, 15, 20]       # test multiple cluster counts
```

**Output structure:**
```
concoct/
├── {assembly_group}/
│   ├── {assembly_group}.db              # Anvi'o contigs database (from metagenome_assemble.smk)
│   ├── {assembly_group}_contigs_stats.txt  # if anvi_taxonomy_and_stats: true
│   ├── .taxonomy.done                   # if anvi_taxonomy_and_stats: true
│   ├── MERGED/                          # Merged profile DB
│   ├── concoct_{nclust}_summary/       # Anvi'o summary + HTML report
│   └── concoct_{nclust}_bins/          # Bin FASTA files for each cluster count
```

**Key differences from MetaBAT2:**
- Requires **at least 2** bowtie2-mapped samples per assembly group (strobealign-only samples are skipped). `anvi-merge` refuses to merge a single profile, so groups with only one bowtie2 sample are silently excluded from CONCOCT's target list — they still go through the normal MetaBAT2/SemiBin2/VAMB track.
- Runs independently of Binette (not combined with other binners)
- Useful for exploratory/manual curation workflows
- Each cluster count produces a complete bin set (not merged or selected)

---

### 3. `summarise_mags.smk`

**Config:** `config_summarise.yaml`

Operates on `.fa.gz` bins in `{output_dir}/renamed_bins/` (or `{output_dir}/binette_renamed_bins/` when `run_binette: true` — see below). Can be run on pre-existing bins without re-running the full pipeline.

Steps:
1. **(Optional) Binette** — refines the raw per-binner bins from `metagenome_binning.smk` (found directly on disk under `output_dir`: `metabat2/{group}/`, `semibin/{group}/output_bins/`, `vamb/{group}/bins/`), selecting the best bins across binners using CheckM2. Requires ≥2 binners to have actually produced bins for a given assembly group. Refined bins are renamed/gzipped into `binette_renamed_bins/`, which CheckM/GTDB-Tk/dRep then use instead of the raw `renamed_bins/`.
2. **CheckM** (v1) — lineage-aware completeness and contamination estimation
3. **GTDB-Tk** — taxonomic classification using the GTDB reference database
4. **dRep** — dereplication at configurable ANI thresholds, using CheckM scores to select representatives

#### Bin naming convention

`metagenome_binning.smk` renames every binner's bins in place, right after
that binner runs, normalizing to a single `.fa.gz` extension so downstream
tools never need to special-case per-binner formats:

| Binner   | Raw filename        | Renamed to                          |
|----------|----------------------|--------------------------------------|
| MetaBAT2 | `bin.N.fa`           | `{assembly_group}_metabat2_bin.N.fa.gz` |
| VAMB     | `{id}.fna`           | `{assembly_group}_vamb_{id}.fa.gz`     |
| SemiBin2 | `SemiBin_N.fa.gz`    | `{assembly_group}_SemiBin_N.fa.gz`     |

`renamed_bins/` (used when `run_binette: false`) is populated by symlinking
(not copying) every binner's already-renamed bins into one directory, so
CheckM/GTDB-Tk/dRep see all binners, not just MetaBAT2.

#### Using bins from outside this pipeline

`summarise_mags.smk` doesn't require its bins to have come from
`metagenome_binning.smk`:
- **Without Binette** (`run_binette: false`): CheckM/dRep/GTDB-Tk only ever
  read `{output_dir}/renamed_bins/*.fa.gz` — drop any gzipped FASTA bins in
  there directly (naming is not enforced) and run the workflow. A valid
  `assemblies_tsv` is still required (parsed unconditionally by
  `common.smk`), but its content is irrelevant to this path since nothing
  here reads assembly contigs.
- **With Binette**: bins need to land in a directory `get_binette_bin_dirs`
  recognizes for the matching `assembly_group` — the three above, or any
  directory listed in the optional `extra_binner_dirs` config key (path
  templates using `{output_dir}`/`{assembly_group}` placeholders), e.g.:
  ```yaml
  extra_binner_dirs:
    - "{output_dir}/binnerY/{assembly_group}/bins"
    - "/absolute/path/binnerZ_results/{assembly_group}"
  ```
  Binette itself doesn't care about filenames, only extensions
  (`.fa`/`.fasta`/`.fna`, optionally gzipped) — bins named however you like
  (`binXXXX.fa`, `binYYY.fa`, ...) are all picked up.

Note: Binette's own CheckM2-based quality report and the separate CheckM1 run
above both estimate completeness/contamination — CheckM1 always runs
regardless of `run_binette`, so genome quality is computed twice when Binette
is enabled. This mirrors the pipeline's existing behavior (previously across
two separate workflow invocations) and is left as-is.

Key config options:
```yaml
output_dir:          "results"

# Binette (multi-binner refinement, optional)
run_binette:         false
binette_checkm2_db:  "/data/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
conda_dirs:
  binette: null   # e.g. /home/ljc444/.conda/envs/binette_env

checkm_db_path:      "/path/to/checkm_db"
checkm_outdir_name:  "checkm1"
checkm_extension:    "gz"

gtdbtk_db_path:      "/path/to/gtdbtk_db"
gtdbtk_extension:    "gz"
gtdbtk_min_af:       0.5
gtdbtk_pplacer_cpus: 8
gtdbtk_write_single_copy_genes: true
gtdbtk_keep_intermediates:      true

drep_outdir_name:    "dRep"
drep_completeness:   66
drep_contamination:  20
drep_primary_ani:    0.95
drep_secondary_ani:  0.98
drep_min_overlap:    0.30
drep_s_algorithm:    "ANIn"
```

Run:
```bash
snakemake -s summarise_mags.smk \
  --profile cluster/ --configfile config_summarise.yaml \
  $SNAKEMAKE_FLAGS --jobs 15 --rerun-incomplete --latency-wait 60
```

---

## Output Structure

```
results/
├── assembly/
│   └── {megahit|spades}/
│       └── {assembly_group}/
│           ├── final.contigs.fa
│           └── final.contigs.reformatted.fa     # if anvi_reformat: true
├── aemb/
│   └── {assembly_group}/
│       └── {sample_id}.tsv            # aemb-style coverage (real strobealign or BAM-derived) — MetaBAT2 cross-sample + VAMB
├── mapping/
│   └── {assembly_group}/
│       ├── {sample_id}.bowtie2.sorted.bam
│       ├── bowtie2_depth.txt
│       └── depth.txt                  # merged depth table (MetaBAT2 input)
├── metabat2/
│   └── {assembly_group}/
│       └── {assembly_group}_metabat2_bin.N.fa.gz   # renamed+gzipped in place
├── semibin/
│   └── {assembly_group}/              # only for groups where every sample is bowtie2-mapped
│       └── output_bins/
│           └── {assembly_group}_SemiBin_N.fa.gz    # renamed in place (already gzipped)
├── vamb/
│   └── {assembly_group}/
│       ├── vae_clusters_unsplit.tsv
│       └── bins/
│           └── {assembly_group}_vamb_{id}.fa.gz    # renamed+gzipped in place
├── renamed_bins/
│   └── {assembly_group}_{binner}_....fa.gz   # symlinks into the three dirs above (metagenome_binning.smk)
├── binette/                            # produced by summarise_mags.smk, if run_binette: true
│   └── {assembly_group}/
│       ├── final_bins_quality_reports.tsv
│       └── final_bins/                # Binette-refined bins
├── binette_renamed_bins/               # produced by summarise_mags.smk, if run_binette: true
│   └── {assembly_group}_binette_binN.fa.gz  # bins CheckM/GTDB-Tk/dRep actually use
├── checkm/
│   └── checkm.txt
├── gtdbtk_classify/
└── dRep/
    ├── dRep.genomeInfo
    └── dereplicated_genomes/
logs/
```

---

## Dependencies

| Environment          | Tools                                          |
|----------------------|------------------------------------------------|
| `envs/assembly.yaml` | SPAdes, MEGAHIT                                |
| `envs/mapping.yaml`  | bowtie2, samtools, strobealign, MetaBAT2 (`jgi_summarize_bam_contig_depths`, used by `aemb_tsv`'s bowtie2 branch) |
| `envs/binning.yaml`  | MetaBAT2 (`jgi_summarize_bam_contig_depths`)   |
| `envs/semibin2.yaml` | SemiBin2                                       |
| `envs/vamb.yaml`     | VAMB (latest)                                  |
| `envs/binette.yaml`  | Binette, CheckM2                               |
| `envs/checkm.yaml`   | CheckM v1                                      |
| `envs/gtdbtk.yaml`   | GTDB-Tk                                        |
| `envs/drep.yaml`     | dRep                                           |

Anvi'o is not installed via yaml — point `conda_dirs.anvio` in the config to an existing installation (e.g. `/home/ljc444/.conda/envs/anvio-9`).

---

## Toy dataset

A toy dataset (3 SPAdes assemblies, 12 subsampled read sets) is available under `toy/`.
All commands below are run from `toy/`, and assume the `SNAKEMAKE_FLAGS` export from [Setup](#conda-env-prefix-for-pipeline-tool-envs).

Five test scenarios are provided, each with its own config pair and output directory.

---

### Scenario 0 — Baseline: pre-built assemblies, MetaBAT2 only

**Config files:** `config_assemble_toy.yaml` + `config_binning_toy.yaml`
**Output:** `results_toy/`

```bash
# Symlink pre-built SPAdes assemblies
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy.yaml --cores 4

# Pre-build missing conda envs (mapping env includes strobealign — built once)
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy.yaml \
  $SNAKEMAKE_FLAGS --conda-create-envs-only --cores 1

# Run
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --rerun-incomplete --latency-wait 60
```

---

### Scenario A — MEGAHIT co-assembly

**Config files:** `config_assemble_toy_coasm.yaml` + `config_binning_toy_coasm.yaml`
**Output:** `results_toy_coasm/`

Assembly groups:
- `coasm_A` — co-assembly of S_144 + S_276 + S_133 (3 samples, de-novo MEGAHIT)
- `spaS135` — single-sample assembly of S_135 (de-novo MEGAHIT)

Mapping: self-mapped samples use `bowtie2`; additional cross-mapped samples use `strobealign`.

```bash
# 1. Assemble with MEGAHIT
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_coasm.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --latency-wait 60

# 2. Bin (MetaBAT2)
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_coasm.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --rerun-incomplete --latency-wait 60
```

---

### Scenario B — Pre-built assemblies + Anvi'o reformat + all binners + Binette

**Config files:** `config_assemble_toy_full.yaml` + `config_binning_toy_full.yaml` + `config_summarise_toy_full.yaml`
**Output:** `results_toy_full/`

Features tested:
- `anvi_reformat: true` (in both `config_assemble_toy_full.yaml` and `config_binning_toy_full.yaml`) — simplify contig headers and filter <1000 bp; for these pre-built assemblies this actually runs in step 1 (`metagenome_assemble.smk`), not step 2
- `binners: [metabat2, semibin2, vamb]` — MetaBAT2 and VAMB run for all three assembly groups; SemiBin2 only runs where every sample is bowtie2-mapped
- Both `mapping_tool` policies from [`mappings.tsv`](#mappingstsv--used-by-workflows-2-and-3) are exercised in this one scenario: `spaS144` has all 5 samples mapped with `bowtie2` (the co-assembly-style policy — SemiBin2-capable), while `spaS276`/`spaS135` keep the mixed bowtie2-self + strobealign-cross policy (SemiBin2 skipped for them)
- `run_binette: true` (in `config_summarise_toy_full.yaml`) — Binette selects the best bins across all three using CheckM2

Set `binette_checkm2_db` in `config_summarise_toy_full.yaml` to the path of your CheckM2 diamond database before running.

```bash
# 1. Symlink assemblies into results_toy_full and reformat them (anvi_reformat: true)
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_full.yaml --cores 4

# 2. Pre-build missing conda envs (mapping, vamb)
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_full.yaml \
  $SNAKEMAKE_FLAGS --conda-create-envs-only --cores 1

# 3. Dry-run to verify DAG
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_full.yaml --cores 1 -n

# 4. Run binning (metabat2, semibin2, vamb)
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_full.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --rerun-incomplete --latency-wait 60

# 5. Refine with Binette + CheckM/GTDB-Tk/dRep
snakemake -s ../summarise_mags.smk --configfile config_summarise_toy_full.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --rerun-incomplete --latency-wait 60
```

---

### Scenario C — CONCOCT binning via Anvi'o (separate manual curation track)

**Config file:** `config_binning_toy_concoct.yaml`  
**Output:** `results_toy_concoct/`

Features tested:
- `run_concoct: true` — parallel CONCOCT clustering workflow (separate from MetaBAT2)
- `concoct_clusters: [10, 15]` — test multiple cluster counts (exported as separate bin sets)
- Manual curation-ready output (bins exported as `.fa` in `concoct/{group}/concoct_{nclust}_*/`)

CONCOCT runs in parallel to the standard binning track but produces independent output — useful for exploring alternative binning solutions before final curation.

```bash
# 1. Symlink assemblies
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_concoct.yaml --cores 4

# 2. Pre-build missing conda envs
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_concoct.yaml \
  $SNAKEMAKE_FLAGS --conda-create-envs-only --cores 1

# 3. Dry-run to verify DAG
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_concoct.yaml --cores 1 -n

# 4. Run locally
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_concoct.yaml \
  $SNAKEMAKE_FLAGS --cores 8 --rerun-incomplete --latency-wait 60

# 4b. Or test on SLURM cluster with minimal jobs (2) for validation before production
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_concoct.yaml \
  --profile ../cluster/ --jobs 2
```

The `--jobs 2` cluster test submits only 2 SLURM jobs concurrently, which is safe for validating the workflow before scaling to production (where you'd use `--jobs 20-50` depending on cluster load).

---

### Scenario D — Full pipeline from raw reads, mixed assemblers, auto-named co-assembly

**Config files:** `config_assemble_toy_realasm.yaml` + `config_binning_toy_realasm.yaml` + `config_summarise_toy_realasm.yaml`
**Output:** `results_toy_realasm/` — see `toy/README.md` for the full writeup (group layout, directory structure, exact commands).

The most complete toy scenario — assembles from raw reads (not pre-built contigs) and exercises nearly every feature above in one run:
- **Per-group assembler** (`assembler` column in `assemblies_toy_realasm.tsv`): `spa_co1` and `spaS276` use MEGAHIT, `spaS135` uses metaSPAdes.
- **Auto-named co-assembly** (`group_key` column): `spa_co1` is a real 5-sample co-assembly, auto-named from its `group_key` rather than hand-typed.
- **All binners + CONCOCT** on the co-assembly group (all bowtie2-mapped → CONCOCT- and SemiBin2-capable); `spaS276`/`spaS135` (bowtie2-self + strobealign-cross) get MetaBAT2 + VAMB only.
- **`anvi_taxonomy_and_stats: true`** — SCG taxonomy + contig stats on the assemble side.
- **`run_binette: true`** — refines all three binners' bins into one final set, then CheckM/GTDB-Tk/dRep.

```bash
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS --rerun-incomplete --latency-wait 60

snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS --rerun-incomplete --latency-wait 60

snakemake -s ../summarise_mags.smk --configfile config_summarise_toy_realasm.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS --rerun-incomplete --latency-wait 60
```

**The bin-count funnel**, from a real validated run (3 groups: 1 co-assembly + 2 single-sample):

| Stage | Count | What happened |
|---|---|---|
| Raw bins across all binners | 114 | MetaBAT2 + VAMB + SemiBin2 (SemiBin2 only for the co-assembly group), summed across all 3 assembly groups |
| After Binette refinement | 26 | Per assembly group, Binette compares the binners' bin sets for the same contigs and selects/merges the best combination — one refined set per group, not a union |
| After dRep dereplication | 24 | Near-identical genomes (e.g. the same organism recovered redundantly by different binners) are collapsed to one representative per ANI cluster |

Each stage is a real reduction step, not just renaming — the counts above will differ run to run (especially with `semibin2_epochs` reduced for toy speed, which directly affects SemiBin2's bin quality/count).

---

### Running on SLURM

Replace the local `--cores 8` run command with:

```bash
snakemake -s ../metagenome_binning.smk --configfile <config_file> --profile ../cluster/ --jobs 20
```

The cluster profile (`cluster/config.yaml`) handles resource allocation, conda env paths, and SLURM submission automatically.

**For validation before production:** use `--jobs 2` to test the full workflow with minimal cluster load:

```bash
snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_concoct.yaml \
  --profile ../cluster/ --jobs 2
```

Once all jobs complete successfully, scale to `--jobs 20-50` for the real dataset.

---

## Troubleshooting

**`error: argument --executor/-e: invalid choice: 'slurm'`**
The SLURM executor plugin isn't installed — see [SLURM executor plugin](#slurm-executor-plugin-required-for---profile-cluster) above.

**`error: unrecognized arguments: --slurm-account=...`**
Account/partition are **job resources**, not CLI flags. Set them under `default-resources` in `cluster/config.yaml`:
```yaml
default-resources:
  slurm_account:   cbmr
  slurm_partition: standardqueue
```

**`No account was given ... sacct: invalid option -- '1'`**
Harmless known bug in `snakemake-executor-plugin-slurm` (≤2.7.1) — its account auto-detection mis-invokes `sacct` and fails, but falls back to submitting under your default association anyway. Set `slurm_account` explicitly (above) to avoid triggering it.

**`A module that was compiled using NumPy 1.x cannot be run in NumPy 2.x` (e.g. in `binette_refine`)**
A leftover `pip install --user` package in `~/.local/lib/pythonX.Y/site-packages` is shadowing the conda env's own (correctly paired) numpy. Diagnose with:
```bash
<env_path>/bin/python -c "import numpy; print(numpy.__file__, numpy.__version__)"
```
If it resolves to `~/.local/...`, fix per-env (without touching `~/.local`, which other tools may need):
```bash
conda env config vars set -p <env_path> PYTHONNOUSERSITE=1
```

**`concoct: ModuleNotFoundError: No module named 'pkg_resources'` / `DistributionNotFound: The 'nose' distribution was not found`**
Building CONCOCT from source (see [CONCOCT binning via Anvi'o](#concoct-binning-via-anvio-optional-parallel-track)) needs `cython` to build, and its runtime `pkg_resources.require("concoct")` check needs `setuptools<81` (removed in newer releases) and `nose` — none of which ship with `anvio-9` by default. Install all three in the `anvio-9` env before building CONCOCT.

**CONCOCT clustering fails with a `TypeError`**
Known incompatibility between CONCOCT and newer scikit-learn ([merenlab/anvio#2154](https://github.com/merenlab/anvio/issues/2154)). Fix in the `anvio-9` env:
```bash
pip install scikit-learn==1.1.0
```

**`SemiBin2 ... Error: abundances from strobealign-aemb can only be used with at least 5 samples` / `KeyError: "None of [Index(['..._1', '..._2', ...])] are in the [index]"`**
SemiBin2's `-a` (strobealign-aemb) mode needs ≥5 samples *and* a SemiBin2-specific split-contig abundance format (`SemiBin2 split_contigs` + aemb against the split FASTA) — plain per-contig depth doesn't have the required `_1`/`_2` rows. This pipeline sidesteps both issues by only running SemiBin2 for groups where every sample is bowtie2-mapped (see [Coverage strategy](#coverage-strategy)), using `-b` directly on the BAMs. If you see either error, the group has a mix of `bowtie2`/`strobealign` samples — that's expected; SemiBin2 is skipped for it by design.

**`vamb bin default: error: unrecognized arguments: --aemb ...` / `MissingOutputException` for `vamb/{group}/clusters.tsv`**
VAMB 5.x removed `--aemb` (needs a single `--abundance_tsv` with header `contigname\t<sample>...`) and renamed its main output from `clusters.tsv` to `vae_clusters_unsplit.tsv`. Already fixed in `vamb_bin`/`_all_targets`; if you see this on a fork/older checkout, update to the current rule.
