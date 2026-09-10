import os
import glob
import re
import shutil
import pandas as pd
from collections import defaultdict

# Setup base configurations
OUT = config.get("output_dir", "results")
ASSEMBLER = config.get("assembler", "megahit")  # "spades" or "megahit"
ASM_TSV = config.get("assemblies_tsv", "assemblies.tsv")
MAP_TSV = config.get("mappings_tsv", None)

# Binner / refinement selection
BINNERS       = config.get("binners", ["metabat2"])   # list: metabat2, semibin2, vamb
RUN_BINETTE   = config.get("run_binette", False)       # merge binner results with Binette
ANVI_REFORMAT = config.get("anvi_reformat", False)     # filter+rename contigs via anvi-script-reformat-fasta
RUN_CONCOCT   = config.get("run_concoct", False)        # CONCOCT binning via Anvi'o (metagenome_binning.smk)
RUN_PLASMAAG  = config.get("run_plasmaag", False)       # PlasMAAG plasmid+MAG binning, invoked as its own Snakemake pipeline (metagenome_binning.smk)
ANVI_TAXONOMY_AND_STATS = config.get("anvi_taxonomy_and_stats", False)  # anvi-run-scg-taxonomy + anvi-display-contigs-stats (metagenome_assemble.smk)

# Optional pre-installed conda env paths. When set, Snakemake uses the existing
# directory instead of creating a new env from the yaml file.
_CONDA_DIRS = config.get("conda_dirs", {})

def conda_env(name, yaml_path):
    """Return pre-existing env path if configured, otherwise fall back to yaml."""
    p = _CONDA_DIRS.get(name)
    if p and str(p).strip().lower() not in ("", "none", "null"):
        return str(p).strip()
    return yaml_path

# Ensure runtime directories exist before any cluster submission.
for d in [OUT, f"{OUT}/logs", f"{OUT}/logs/per_rule"]:
    os.makedirs(d, exist_ok=True)

# Parse Sourcing Sheet
asm_df = pd.read_csv(ASM_TSV, sep="\t").fillna("None")
asm_df.columns = [c.strip() for c in asm_df.columns]

PREBUILT_ASM = {}
BUILT_ASM_GROUPS = set()
ASM_READS_R1 = defaultdict(list)
ASM_READS_R2 = defaultdict(list)

# Per-group assembler override: an optional "assembler" column in
# assemblies_tsv (same convention as mapping_tool in mappings_tsv) lets
# different assembly groups use different assemblers in the same run.
# Falls back to the global config's assembler for any row/TSV without it.
HAS_ASSEMBLER_COL = "assembler" in asm_df.columns

# Short labels used by auto-naming below (e.g. "spa" for spades, matching
# this dataset's existing "spaS144"-style convention). Override/extend via
# the assembler_abbrev config key.
ASSEMBLER_ABBREV = config.get("assembler_abbrev", {"spades": "spa", "megahit": "mh"})

# Optional auto-naming: an optional "group_key" column in assemblies_tsv
# groups rows that should be assembled together (a shared key = co-assembly;
# a key used by only one row = single-sample assembly) independently of the
# final assembly_group name. Leave assembly_group blank ("None"/empty/"auto")
# on those rows to auto-generate it as "{assembler_abbrev}_co{n}" for
# co-assemblies (n increments per assembler) or "{assembler_abbrev}_{group_key}"
# for single-sample assemblies. Give assembly_group explicitly to opt out.
if "group_key" in asm_df.columns:
    _group_key_counts = asm_df["group_key"].value_counts().to_dict()
    _co_counters = defaultdict(int)
    _group_key_name_cache = {}

    def _assembler_for_row(row):
        row_assembler = str(row["assembler"]).strip() if HAS_ASSEMBLER_COL else ""
        return row_assembler if row_assembler.lower() not in ("", "none") else ASSEMBLER

    def _auto_name_for_row(row):
        gk = row["group_key"]
        if gk not in _group_key_name_cache:
            abbrev = ASSEMBLER_ABBREV.get(_assembler_for_row(row), _assembler_for_row(row))
            if _group_key_counts[gk] > 1:
                _co_counters[abbrev] += 1
                _group_key_name_cache[gk] = f"{abbrev}_co{_co_counters[abbrev]}"
            else:
                _group_key_name_cache[gk] = f"{abbrev}_{gk}"
        return _group_key_name_cache[gk]

    def _resolve_assembly_group(row):
        val = str(row["assembly_group"]).strip()
        return _auto_name_for_row(row) if val.lower() in ("", "none", "auto") else val

    asm_df["assembly_group"] = asm_df.apply(_resolve_assembly_group, axis=1)

ASSEMBLER_FOR = {}

for _, row in asm_df.iterrows():
    grp = row["assembly_group"]
    path = str(row["assembly_path"]).strip()

    row_assembler = str(row["assembler"]).strip() if HAS_ASSEMBLER_COL else ""
    ASSEMBLER_FOR[grp] = row_assembler if row_assembler.lower() not in ("", "none") else ASSEMBLER

    if path and path.lower() != "none" and path != "":
        PREBUILT_ASM[grp] = path
    else:
        BUILT_ASM_GROUPS.add(grp)
        if str(row["fq1"]).lower() != "none":
            ASM_READS_R1[grp].append(row["fq1"])
            ASM_READS_R2[grp].append(row["fq2"])

# Parse Mapping Instructions (optional — not needed for assemble-only workflow)
MAP_TSV = config.get("mappings_tsv", None)

MAPPING_READS = {}
MAPPING_TOOL_FOR = {}
ASM_TO_MAPPED_SAMPLES = defaultdict(list)

if MAP_TSV and os.path.exists(MAP_TSV):
    map_df = pd.read_csv(MAP_TSV, sep="\t").fillna("None")
    map_df.columns = [c.strip() for c in map_df.columns]

    MAPPING_READS = {row["sample_id"]: (row["fq1"], row["fq2"]) for _, row in map_df.iterrows()}
    # Keyed by (assembly_group, sample_id): a sample can use different tools for different groups.
    MAPPING_TOOL_FOR = {(row["assembly_group"], row["sample_id"]): row["mapping_tool"] for _, row in map_df.iterrows()}

    for _, row in map_df.iterrows():
        ASM_TO_MAPPED_SAMPLES[row["assembly_group"]].append(row["sample_id"])

# Wildcard constraints protect against duplicate generation conflicts
prebuilt_constraint = "|".join(map(re.escape, PREBUILT_ASM.keys())) if PREBUILT_ASM else "none_placeholder"
built_constraint = "|".join(map(re.escape, BUILT_ASM_GROUPS)) if BUILT_ASM_GROUPS else "none_placeholder"

# Helper list of all assembly groups
all_groups = list(PREBUILT_ASM.keys()) + list(BUILT_ASM_GROUPS)

# Assembly groups that actually have mapped samples (i.e. can be binned).
# Shared by metagenome_binning.smk and summarise_mags.smk.
binnable_groups = [g for g in all_groups if g in ASM_TO_MAPPED_SAMPLES]

# anvi-merge refuses to merge a single profile, so CONCOCT (which needs a
# merged profile) is only usable for groups with >=2 bowtie2-mapped samples.
# Shared by metagenome_assemble.smk (contigs DB creation) and
# metagenome_binning.smk (the rest of the CONCOCT track).
concoct_capable_groups = [
    g for g in binnable_groups
    if sum(1 for s in ASM_TO_MAPPED_SAMPLES[g] if MAPPING_TOOL_FOR[(g, s)] == "bowtie2") >= 2
]

# snakemake -s metagenome_binning.smk --configfile config_binning.yaml --profile cluster --use-conda
# snakemake -s metagenome_assemble.smk --configfile config_assemble.yaml --profile cluster --use-conda
# snakemake -s summarise_mags.smk --configfile config_summarise.yaml --profile cluster --use-conda