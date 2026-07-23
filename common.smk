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
for d in [OUT, "logs", "logs/per_rule"]:
    os.makedirs(d, exist_ok=True)

# Parse Sourcing Sheet
asm_df = pd.read_csv(ASM_TSV, sep="\t").fillna("None")
asm_df.columns = [c.strip() for c in asm_df.columns]

PREBUILT_ASM = {}
BUILT_ASM_GROUPS = set()
ASM_READS_R1 = defaultdict(list)
ASM_READS_R2 = defaultdict(list)

for _, row in asm_df.iterrows():
    grp = row["assembly_group"]
    path = str(row["assembly_path"]).strip()
    
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

# snakemake -s metagenome_binning.smk --configfile config_binning.yaml --profile cluster --use-conda
# snakemake -s metagenome_assemble.smk --configfile config_assemble.yaml --profile cluster --use-conda
# snakemake -s summarise_mags.smk --configfile config_summarise.yaml --profile cluster --use-conda