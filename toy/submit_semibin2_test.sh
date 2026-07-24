#!/bin/bash
#SBATCH --job-name=semibin2_spaS144_test
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/submit_semibin2_test.%j.log

# Runs the snakemake driver itself as a real SLURM job, fully detached from
# any interactive session, so it (and the semibin2_bin job it submits)
# survives session disconnects. --rerun-triggers mtime avoids re-running
# already-correct MetaBAT2/VAMB steps; only semibin2_bin is actually missing.

set -e
set -o pipefail
cd "$SLURM_SUBMIT_DIR"

source /home/ljc444/.conda/envs/snakemake/bin/activate

snakemake -s ../metagenome_binning.smk \
  --configfile config_binning_toy_full.yaml \
  --profile ../cluster/ \
  --jobs 4 \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --latency-wait 60
