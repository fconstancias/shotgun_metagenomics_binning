# shotgun_metagenomics_binning

A modular Snakemake pipeline for metagenomic assembly, read mapping, binning, and MAG dereplication/taxonomy. The three workflows can be run independently or chained sequentially.

---
## TODO:

- Branch with CoverM for mapping with option to use galah instead of Drep in summarise_mags
- later: use Simka / SimkaMin to select samples for mapping
- later: ONT data ?
- Rename `run_concoct` (config key, both `config_assemble`/`config_binning`,
  read once from `common.smk`) to something like `run_anviodb_concoct`, and
  the `{output_dir}/concoct/{group}/` directory (metagenome_assemble.smk) to
  something like `anvio_concoct/`. Both are named after CONCOCT specifically,
  but the thing they actually gate/hold is the shared Anvi'o contigs DB
  (`anvi_gen_contigs_db`/`anvi_run_hmms`), built whenever either
  `run_concoct` OR `anvi_taxonomy_and_stats` is true — so an assembly-only
  run with `anvi_taxonomy_and_stats: true`, `run_concoct: false` (e.g.
  single-sample groups, which can never be CONCOCT-capable) still produces a
  `concoct/{group}/` folder per group that CONCOCT itself will never touch,
  confusing to read at a glance. Naming both after what they actually are
  (the shared Anvi'o DB prerequisite) rather than one of the two possible
  downstream consumers (CONCOCT vs SCG-taxonomy/stats) would be clearer.
  Real rename, touches every existing config file — do deliberately, not
  as a drive-by.
- `spades_assemble` (metagenome_assemble.smk) only keeps `scaffolds.fasta`
  (renamed to `final.contigs.fa`) and deletes the rest of SPAdes's output
  dir (`rm -rf {tmpdir}` after the `mv`) — including the assembly graph
  (`assembly_graph_with_scaffolds.gfa`/`assembly_graph.fastg`) and the
  unscaffolded `contigs.fasta`. Consider keeping the graph file too (e.g.
  `final.assembly_graph.gfa` alongside `final.contigs.fa`) for graph-aware
  downstream uses (Bandage visualization, repeat/strain resolution,
  graph-based binning) instead of discarding it unconditionally.
- Split `metagenome_binning.smk` into separate mapping and binning
  workflows. Right now one file does both (bowtie2/strobealign mapping
  through to MetaBAT2/VAMB/SemiBin2/CONCOCT), so re-tuning binner
  parameters means the mapping DAG is re-evaluated too (even though it's
  a no-op via Snakemake's up-to-date check, it's still one combined
  workflow to reason about). Same rationale as moving Binette out to
  summarise_mags.smk earlier — decoupling stages that have genuinely
  different iteration cadences (map once, re-bin many times while tuning)
  makes each stage easier to re-run/tune independently.

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
anvi-setup-scg-taxonomy -T 8
```

Point `checkm_db_path`, `gtdbtk_db_path`, and `binette_checkm2_db` in your
config at wherever you put these.

#### Anvi'o SCG taxonomy in detail

`anvi-setup-scg-taxonomy` downloads GTDB's single-copy-gene reference data
and builds a diamond (`.dmnd`) search database for each of the 22 marker
genes anvi'o uses (`Bacteria_71`/`Archaea_76`/ribosomal proteins/etc.) —
this is what `anvi-run-hmms` + `anvi-run-scg-taxonomy` need at runtime.

- **One-time per anvi'o conda env**, not per-project or per-output-dir —
  run it once for whichever env `conda_dirs.anvio` points at, and every
  future run using that same env is covered, no need to repeat it.
- Needs outbound internet access from wherever it runs (login node or a
  compute node with internet, not purely internal cluster storage).
- By default installs into that anvi'o installation's own package data
  directory. If you don't have write access there, or want a shared
  location, pass `--scgs-taxonomy-data-dir /shared/path` — you must then
  pass that **same** flag to every downstream anvi'o SCG command too
  (not something this pipeline currently threads through, so in practice
  stick with the default location unless you have a specific reason not to).
- `--gtdb-release <N>` pins a specific GTDB release instead of the latest.
- **Verify it worked**: count the diamond databases —
  ```bash
  find $(python3 -c "import anvio,os; print(os.path.dirname(anvio.__file__))")/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES -iname "*.dmnd" | wc -l
  # should print 22
  ```
  (or wherever `--scgs-taxonomy-data-dir` pointed, if you set it).
- **If skipped**: `anvi-run-scg-taxonomy` fails immediately with `Config
  Error: ... missing 22 of 22 databases required` — a clear, unambiguous
  signal to come back and run this step.

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
  templates using `{output_dir}`/`{assembly_group}` placeholders), or
  `extra_binner_dirs_by_group` (a plain dict keyed by the group in *this*
  run, mapping to already-resolved absolute paths — for when group naming
  doesn't match a template, e.g. a different binner run with its own
  naming convention), e.g.:
  ```yaml
  extra_binner_dirs:
    - "{output_dir}/binnerY/{assembly_group}/bins"
    - "/absolute/path/binnerZ_results/{assembly_group}"
  extra_binner_dirs_by_group:
    spa_co1:
      - "/absolute/path/some_other_binner_run/spa_co1_output"
  ```
  Binette itself doesn't care about filenames, only extensions
  (`.fa`/`.fasta`/`.fna`, optionally gzipped) — bins named however you like
  (`binXXXX.fa`, `binYYY.fa`, ...) are all picked up.

  **Hard constraint**: every bin passed to Binette must reference contigs
  from *this group's own* assembly FASTA (the same file passed as `-c`) —
  Binette errors out (`ValueError: N contigs from the input bins were not
  found in the contigs file`) otherwise. This only ever holds for another
  binner run on the *same* assembly as this group — never for bins from a
  genuinely different assembly (see below).

#### Comparing bins across separate pipeline runs (e.g. co-assembly vs single-sample per participant)

If a participant was binned two different ways in two separate runs — say
a co-assembly run and a single-sample run — those two runs produced
**different assemblies** (different contigs, even from the same reads), so
their bins can never be fed into one Binette call together (see the hard
constraint above). The right place to compare them is one level
downstream, at CheckM/dRep/GTDB-Tk: those tools work on standalone genome
FASTA files with no shared-contig requirement, and dRep specifically
compares genomes by whole-sequence ANI — exactly the right tool for "is
this co-assembly bin the same organism as that single-sample bin."

Use `extra_bin_source_dirs`: a flat list of other runs' bin-source
directories (their own `binette_renamed_bins/` if they used Binette, or
`renamed_bins/` if not) to pool in wholesale alongside this run's own bins,
before CheckM/dRep/GTDB-Tk run:
```yaml
extra_bin_source_dirs:
  - "/abs/path/results_singlesample/renamed_bins"
  - "/abs/path/results_singlesample/binette_renamed_bins"
```
No group/sample matching is needed here — dRep's ANI dereplication safely
collapses genomes that really are the same organism across runs/approaches
and ignores everything else, so you can pool in an *entire* other run
(all its participants) without needing to cherry-pick which groups
correspond to which. Bins are symlinked (not copied) into
`{output_dir}/pooled_bins/`, which then becomes what CheckM/dRep/GTDB-Tk
actually scan.

Note: Binette's own CheckM2-based quality report and the separate CheckM1 run
above both estimate completeness/contamination — by default CheckM1 always
runs regardless of `run_binette`, so genome quality is computed twice when
Binette is enabled. Set `use_binette_checkm2_for_drep: true` to skip the
redundant CheckM1 pass and build dRep's `genomeInfo` directly from Binette's
own `final_bins_quality_reports.tsv` (per group) instead — only takes effect
when `run_binette: true` (nothing to reuse otherwise), off by default so
existing configs keep today's behavior unchanged.

##### Scoped per-participant comparison: `per_group_drep_comparisons`

`extra_bin_source_dirs` pools bins **globally** — every bin from the other
run joins every bin from this one before a single dataset-wide dRep call.
That's the right tool when you want one dereplicated set across the whole
dataset, but it's the wrong tool when what you actually want is: "for this
one participant, compare their co-assembly bins against their single-sample
bins, without mixing in anyone else's bins." For that, use
`per_group_drep_comparisons` instead — a per-group dict that runs its own
scoped CheckM + dRep + GTDB-Tk, alongside (not instead of) the global ones:

```yaml
per_group_drep_comparisons:
  spa_co1:                                    # this run's assembly_group to scope the comparison to
    output_dir:     "results_singlesample"    # the other run's output_dir
    assemblies_tsv: "assemblies_singlesample.tsv"   # the other run's assemblies_tsv
```

The other run's matching groups are auto-detected the same way as
`extra_binner_dirs_by_group` — by reading the other run's `assemblies_tsv`
and matching `fq1` paths against this group's constituent samples (via
`mappings_tsv`), so `spa_co1` (a co-assembly of several samples) correctly
pulls in every single-sample group from the other run built from any of
those same fastqs, regardless of how the two runs name their groups.

For each configured group, this produces its own subtree:
```
{output_dir}/per_group_drep/{group}/
├── pooled_bins/            # symlinks: this group's own bins + matched groups' bins from the other run
├── checkm1/checkm.txt
├── dRep/dereplicated_genomes/
└── gtdbtk_classify/
```
Pooling prefers each side's Binette-refined bins (`binette_renamed_bins/`)
and only falls back to raw per-binner bins (`renamed_bins/`) for a group that
has no Binette output — mirroring how this run's own bin source is chosen.
Because that preference is resolved by scanning the other run's directory at
pooling time (not by a Snakemake-tracked dependency across the two separate
workflow invocations), if the other run gets a fresh Binette pass *after*
`pooled_bins/` was already built here, delete
`{output_dir}/per_group_drep/{group}/pooled_bins/` and re-run to pick it up.

This is fully additive: omit `per_group_drep_comparisons` (the default) and
nothing changes — only the global CheckM/dRep/GTDB-Tk run. Both mechanisms
can be used together, and both can be combined with `extra_bin_source_dirs`
in the same run.

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

Six test scenarios are provided, each with its own config pair and output directory.

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

#### `build_anvio_profile_db` — anvi'o profile-db without CONCOCT

`run_concoct: true` above builds a full CONCOCT clustering track, which
requires `anvi-merge`'d profile-dbs and therefore >=2 bowtie2-mapped
samples per group. If you just want a profile-db (for `anvi-summarize`
gene coverage/detection export, GFF3/protein-FASTA export via
`anvi-get-sequences-for-gene-calls`, DESMAN's `anvi-gen-variability-profile`
SNV step, `anvi-interactive`, etc.) without committing to CONCOCT
clustering, set `build_anvio_profile_db: true` instead (independent of
`run_concoct` — both can be on at once with no conflict, since a group
already CONCOCT-clustered already has its profile-db and this is a
no-op for it):

```yaml
build_anvio_profile_db: true
```

For every binnable group with at least one bowtie2-mapped sample, this
targets `concoct/{group}/.hmms.done` (HMM annotation on the contigs-db)
and a final profile-db:
- **>=2 bowtie2 samples**: `concoct/{group}/MERGED/PROFILE.db` (via the
  same `anvi_merge` rule CONCOCT uses).
- **Exactly 1 bowtie2 sample** (the common case for single-sample-assembly
  projects, where that one sample is the group's own self-sample):
  `concoct/{group}/PROFILE_{sample}/` directly — `anvi-merge` itself
  refuses to merge a single profile, so there's nothing to merge.

Contigs-db itself (`concoct/{group}/{group}.db`) is built in
`metagenome_assemble.smk`, shared with the optional SCG-taxonomy/contig-stats
step there (`anvi_taxonomy_and_stats: true`) — if that was already set for
the assemble run, the contigs-db already exists for every group regardless
of this flag.

**Exporting a GFF3 + protein FASTA from the resulting contigs-db** (e.g.
as pre-computed input to a downstream annotation pipeline like EBI's
mobilome-annotation-pipeline, whose samplesheet accepts optional
`proteins_gff`/`proteins_faa` columns to skip its own gene-calling and
keep gene IDs consistent with this pipeline's own database):

```bash
anvi-get-sequences-for-gene-calls -c concoct/{group}/{group}.db --export-gff3 -o {group}.gff3
anvi-get-sequences-for-gene-calls -c concoct/{group}/{group}.db --get-aa-sequences -o {group}.faa \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}"
```

Two separate calls -- `--export-gff3` and `--get-aa-sequences` are
mutually exclusive in one invocation (confirmed: anvi'o errors with "AA
sequences can only be reported in FASTA format, please remove the
--export-gff3 flag"). The GFF3's attributes column is bare (just a gene
ID, e.g. `ID=spaS222___0`) -- no function/product annotation baked in
even if HMM/COG/KOfam annotation was run on the contigs-db, so this only
carries gene *structure* (coordinates + sequence), not descriptions.

**`--defline-format` on the second call is required, not optional** -- without it, anvi'o
defaults the FASTA defline to bare `{gene_caller_id}` (e.g. `>0`), while `--export-gff3`
always writes `ID={contigs_db_project_name}___{gene_caller_id}` (e.g. `ID=spaS222___0`) --
same gene, two different ID strings, with no way to reconcile them downstream except by
matching row order. A downstream pipeline that joins by ID string (not row position) will
silently get wrong/empty results, not an error -- hit this for real feeding an early export
into EBI's mobilome-annotation-pipeline and nf-core/funcscan, both of which validate file
existence/format on their schemas but never checked ID correspondence. `--export-gff3`
itself rejects `--defline-format` outright ("not compatible with the GFF3 output mode") --
so the fix only goes on the `--get-aa-sequences` call, matching what `--export-gff3` already
writes.

#### PlasMAAG — plasmid + MAG binning, run as its own pipeline (`dev_plasmaag` branch)

[PlasMAAG](https://github.com/RasmussenLab/PlasMAAG) ([Nature Biotechnology
2026](https://www.nature.com/articles/s41587-026-03005-7)) recovers plasmids
*and* cellular genomes together from assembly-alignment graphs + contrastive
learning. It's a full multi-stage Snakemake pipeline in its own right (BLAST
all-vs-all across samples, node2vec, a VAMB-based VAE, geNomad classification)
and bundles VAMB itself as a component — not designed for DAS_Tool/Binette-
style ensembling with MetaBAT2/CONCOCT/SemiBin2 anyway. Given that, it's
**not** wired into `metagenome_binning.smk` as a rule: wrapping a whole
independent Snakemake pipeline inside one rule of this one would flatten its
own internal parallelism into a single opaque job, and it doesn't need to
share a DAG with the other binners since it never ensembles with them. Same
pattern as MAP/funcscan/anvi'o-KEGG elsewhere in this project — standalone
pipelines that consume this pipeline's outputs, not rules inside it.

It does its **own** internal read mapping (strobealign — the same mapper
this pipeline already uses elsewhere — but a separate pass from our own
BAMs; its CLI has no option to reuse a pre-computed BAM instead).

**This pipeline's only role**: retain the 3 files SPAdes normally discards
that PlasMAAG needs. One assemble-stage config flag:
```yaml
# metagenome_assemble.smk config -- SPAdes groups only, no effect on MEGAHIT
keep_plasmaag_files: true
```
retains SPAdes' own `contigs.fasta` (pre-scaffolding — **not** the same file
as this pipeline's renamed `final.contigs.fa`), `assembly_graph_after_simplification.gfa`,
and `contigs.paths` under `assembly/spades/{group}/plasmaag_input/` — all
three mutually consistent, all three normally deleted with the rest of the
SPAdes tmpdir. Default off (same "don't keep what's not needed" policy as
everywhere else — this is real extra disk per group). MEGAHIT groups need no
extra retention; PlasMAAG's `--reads_and_contigs` mode works from
`final.contigs.fa` directly.

**Fully independent of `anvi_reformat`** — `keep_plasmaag_files` copies
straight out of `rule spades_assemble`'s own SPAdes tmpdir, before
`reformat_contigs` (a separate, downstream rule operating on the different
file `final.contigs.fa` → `final.contigs.reformatted.fa`) ever runs.
`anvi_reformat: true` is this pipeline's normal default for real runs (see
Scenario B above) and works alongside `keep_plasmaag_files: true` with no
interaction at all — don't disable it for PlasMAAG's sake. (The validation
run this feature was proven against happened to use `anvi_reformat: false`
purely to keep that particular test minimal, not because it's required —
confirmed directly by tracing the rule dependency graph, not just assumed.)

**Two sharp edges, both confirmed the hard way against real `spa_single_all` data:**

1. **Cannot retroactively recover files for already-built assemblies.**
   `keep_plasmaag_files` only captures files while `spades_assemble` is
   actually running. Any group whose assembly came from a prebuilt path
   (`assembly_path` set in `assemblies_tsv`, common for assemblies ingested
   from before this pipeline existed) never runs that rule at all — checked
   `spa_single_all`'s own prebuilt directories directly: only the final
   scaffolds file and a report remain, the original SPAdes working
   directory (and its graph/paths files) is long gone. For an existing
   project like this one, using PlasMAAG on already-built groups means a
   genuine full re-assembly from raw reads, not a cheap incremental add.
2. **Turning this on retroactively re-triggers assembly for *every* already-built
   group**, even ones that don't care about PlasMAAG. `plasmaag_input` is a
   new declared output of `spades_assemble` — Snakemake reruns a rule if
   *any* of its declared outputs are missing, and this one never existed
   before, so every previously-built group gets re-assembled from scratch
   the first time this flag is flipped on. Confirmed directly via dry-run
   (`spades_assemble` scheduled to rerun for a group whose `final.contigs.fa`
   already existed, solely because `plasmaag_input` didn't).

**Running PlasMAAG itself**, once `keep_plasmaag_files` has produced real
input (or against MEGAHIT's `final.contigs.fa` directly): install it
separately (`git clone` + its own `conda env create`, not installable from a
plain conda dependency list — see its README), then build its samplesheet
with `scripts/write_plasmaag_samplesheet.py` and invoke it directly:
```bash
conda activate plasmaag
python3 scripts/write_plasmaag_samplesheet.py \
    --r1 <sample1_R1> [<sample2_R1> ...] --r2 <sample1_R2> [<sample2_R2> ...] \
    --assembly-path results/assembly/spades/{group}/plasmaag_input \
    --third-column assembly_dir --out samplesheet.tsv
PlasMAAG --reads_and_assembly_dir samplesheet.tsv --output plasmaag_out/{group} --threads 16
```
(`--third-column contigs` + point `--assembly-path` at `final.contigs.fa` for
MEGAHIT groups, per PlasMAAG's own two input modes.) For real cluster runs,
PlasMAAG's own README documents submitting via its Snakemake directly with
`--executor cluster-generic` instead of the CLI wrapper.

**Mapping strategy / one samplesheet row per sample**: multiple rows (a
co-assembly group with several contributing samples, from `ASM_READS_R1`/
`ASM_READS_R2` in `assemblies_tsv`) is how PlasMAAG gets its main
cross-sample alignment-graph advantage; a single-assembly group typically has
just one row, which still runs correctly but forgoes that advantage (the
same single- vs. multi-sample tradeoff VAMB itself has always had).

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

### Scenario E — `extra_bin_source_dirs`: comparing bins across separate runs

**Config files:** `config_assemble_toy_crossrun.yaml` + `config_binning_toy_crossrun.yaml`
**Output:** `results_toy_crossrun/`

A minimal, fast companion to Scenario D that validates `extra_bin_source_dirs`
(see [Comparing bins across separate pipeline runs](#comparing-bins-across-separate-pipeline-runs-eg-co-assembly-vs-single-sample-per-participant)):
`single_144`/`single_178` — deliberately named differently from Scenario D's
`spa_co1` — are single-sample MEGAHIT assemblies of 2 of `spa_co1`'s 5
constituent samples (the exact same fastq files), binned with MetaBAT2 +
VAMB only (kept small on purpose; this scenario exists to exercise the
cross-run pooling mechanism, not to be biologically meaningful on its own).

```bash
snakemake -s ../metagenome_assemble.smk --configfile config_assemble_toy_crossrun.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS --rerun-incomplete --latency-wait 60

snakemake -s ../metagenome_binning.smk --configfile config_binning_toy_crossrun.yaml \
  --profile ../cluster/ $SNAKEMAKE_FLAGS --rerun-incomplete --latency-wait 60
```

Then point a `summarise_mags.smk` config for the Scenario D run at it:
```yaml
extra_bin_source_dirs:
  - "results_toy_crossrun/renamed_bins"
```
Validated end-to-end: pooling correctly combined Scenario D's own 26
Binette-refined bins with 35 raw bins from this separate run (19 from
`single_144`, 16 from `single_178`), all scored together by CheckM with
real, varied completeness/contamination values — confirming bins from a
genuinely different assembly can be compared downstream of Binette (where
whole-genome ANI applies) even though they can never be fed into Binette
itself together (different contigs).

**`per_group_drep_comparisons` validation:** after also running Binette on
Scenario E itself (`config_summarise_toy_crossrun.yaml`, `run_binette: true`),
a `per_group_drep_comparisons` entry scoping `spa_co1` to this run's
`single_144`/`single_178` was validated end-to-end:
```yaml
per_group_drep_comparisons:
  spa_co1:
    output_dir:     "results_toy_crossrun"
    assemblies_tsv: "assemblies_toy_crossrun.tsv"
```
`fq1` auto-detection correctly matched only `single_144` + `single_178` (not
Scenario D's other participants). Pooling correctly preferred each side's
Binette-refined bins over raw ones: 15 (`spa_co1`) + 6 (`single_144`) + 8
(`single_178`) = 29 bins, no duplication. dRep dereplicated 29 → 24 genomes,
correctly collapsing species recovered by both assembly strategies for the
same participant (e.g. *Faecalibacterium longum*, *Roseburia rectalis*) while
keeping species unique to the co-assembly (*Prevotella copri*, *Akkermansia*,
*Bifidobacterium infantis*) — exactly the "which assembly strategy recovered
this species" comparison the feature is for. GTDB-Tk classified all 29 pooled
bins to species level. The global (non-scoped) CheckM/dRep/GTDB-Tk for
Scenario D ran unaffected alongside it in the same invocation.

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

**Capping concurrent assemblies specifically:** `--jobs`/`jobs:` caps the
total number of SLURM jobs in flight across *every* rule combined — it
doesn't let you cap just the heavy, long-running ones (`spades_assemble`,
`megahit_assemble`, each 16 threads/64GB/up to 24h) separately from cheap,
fast rules (`reformat_contigs`, `use_prebuilt_assembly`, etc.). Both
assembly rules declare `resources: assembly_slots = 1`; `cluster/config.yaml`
sets a default pool of 10 (`resources: ["assembly_slots=10"]`), so no more
than 10 SPAdes/MEGAHIT assemblies run concurrently regardless of the overall
`--jobs` value — useful when assembling many single-sample groups at once on
a shared cluster. Override per-invocation with `--resources assembly_slots=N`.

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
