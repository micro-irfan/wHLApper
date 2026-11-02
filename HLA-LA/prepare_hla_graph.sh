#!/usr/bin/env bash
#
# Index (prepare) an HLA*LA PRG graph with the HLA*LA Docker image:
#   HLA-LA --action prepareGraph --PRG_graph_dir <graph_dir>
#
# The graph folder is modified in place (serializedGRAPH etc. are written into it),
# and can then be passed to the pipeline with --graph_dir.
#
# Example:
#   ./prepare_hla_graph.sh -g /ref/graph/PRG_MHC_GRCh38_withIMGT
#
# Requires: docker, and the image built from docker/Dockerfile (in this folder)

set -euo pipefail

IMAGE="hla-la:1.0.3"
GRAPH_DIR=""
FORCE=0

usage() {
    cat <<EOF
Usage: $(basename "$0") -g <PRG_graph_dir> [-i image] [-f]

  -g  HLA*LA PRG graph folder to index, e.g. PRG_MHC_GRCh38_withIMGT (required)
  -i  Docker image with HLA*LA [$IMAGE]
  -f  Re-index even if the graph already looks prepared
EOF
    exit "${1:-0}"
}

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

while getopts "g:i:fh" opt; do
    case $opt in
        g) GRAPH_DIR=$OPTARG ;;
        i) IMAGE=$OPTARG ;;
        f) FORCE=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
[[ -n $GRAPH_DIR ]] || usage 1
[[ -d $GRAPH_DIR ]] || die "Graph folder not found: $GRAPH_DIR"
[[ -w $GRAPH_DIR ]] || die "Graph folder is not writable: $GRAPH_DIR"
command -v docker >/dev/null || die "docker not found on PATH"
docker image inspect "$IMAGE" >/dev/null 2>&1 || die "Docker image not found: $IMAGE (build it with: docker build -t $IMAGE $(dirname "$0")/docker)"

GRAPH_DIR=$(realpath "$GRAPH_DIR")

if [[ -e $GRAPH_DIR/serializedGRAPH && $FORCE -eq 0 ]]; then
    log "$GRAPH_DIR already contains serializedGRAPH - looks prepared. Use -f to re-index."
    exit 0
fi

# Mounted at the same path so HLA*LA's messages refer to the real location
log "Preparing graph $GRAPH_DIR with $IMAGE (this takes a while and needs a lot of memory)"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$GRAPH_DIR:$GRAPH_DIR" \
    "$IMAGE" \
    HLA-LA --action prepareGraph --PRG_graph_dir "$GRAPH_DIR"

[[ -e $GRAPH_DIR/serializedGRAPH ]] || die "prepareGraph finished but $GRAPH_DIR/serializedGRAPH was not created"
log "Done. Use it with: --graph_dir $GRAPH_DIR"
