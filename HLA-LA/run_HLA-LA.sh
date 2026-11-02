#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Directory of the HLA*LA workflow (defaults to the folder this script is in).
# Set WORKFLOW_DIR if you copy this script into a run folder.
WORKFLOW_DIR="${WORKFLOW_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

RUN_NAME="hlala_ont.ont_open_data"
SAMPLESHEET="samplesheet.csv"
PROFILE="docker"                 # docker | singularity | conda, optionally ",slurm"

# References - edit to match your setup
REFERENCE="/data/irfan/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"   # what the BAMs are aligned to
HLA_REFERENCE="/data/irfan/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"     # decoy the graph was built on
HLA_REFERENCE_MMI="/data/irfan/ref/GRCh38_full_analysis_set_plus_decoy_hla.mmi" # optional, "" to build per sample
GRAPH_DIR="/data/irfan/ref/graph/PRG_MHC_GRCh38_withIMGT"                      # prepared with prepare_hla_graph.sh
HLA_CONTAINER="hla-la:1.0.3"     # image built from docker/Dockerfile

# Resources - must fit the machine (check with: free -g, nproc)
THREADS=8
HLA_MEMORY="60 GB"               # HLA*LA task
ALIGN_MEMORY="32 GB"             # minimap2 realignment

[[ -f "${WORKFLOW_DIR}/main.nf" ]] || { echo "main.nf not found in ${WORKFLOW_DIR} - set WORKFLOW_DIR" >&2; exit 1; }
[[ -f "${SAMPLESHEET}" ]] || { echo "Samplesheet not found: ${SAMPLESHEET}" >&2; exit 1; }

MMI_ARG=()
[[ -n "${HLA_REFERENCE_MMI}" ]] && MMI_ARG=(--hla_reference_mmi "${HLA_REFERENCE_MMI}")

echo "Launching Nextflow in the background..."
NXF_VER=26.04.4 nextflow run "${WORKFLOW_DIR}/main.nf" \
    --input "${SAMPLESHEET}" \
    -profile "${PROFILE}" \
    --reference "${REFERENCE}" \
    --hla_reference "${HLA_REFERENCE}" \
    "${MMI_ARG[@]}" \
    --graph_dir "${GRAPH_DIR}" \
    --hla_container "${HLA_CONTAINER}" \
    --threads "${THREADS}" \
    --hla_memory "${HLA_MEMORY}" \
    --align_memory "${ALIGN_MEMORY}" \
    --min_coverage 20 \
    --min_q1 0.95 \
    --outdir "results/${RUN_NAME}" \
    -resume > "nextflow_run.${RUN_NAME}.log" 2>&1 &

echo "PID: $!  |  log: nextflow_run.${RUN_NAME}.log"
