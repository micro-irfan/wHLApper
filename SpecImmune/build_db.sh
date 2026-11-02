#!/usr/bin/env bash
#
# build_db.sh - build the SpecImmune database with Docker (scripts/make_db.py)
#
# Families: HLA, KIR, IG_TR, CYP. Each is built with:
#   python scripts/make_db.py -o <db> -i <family>
# CYP also needs the PharmVar database (https://www.pharmvar.org/download, unzipped):
# all its *.fasta files are merged and passed as --CYP_fa.
#
# Example:
#   ./build_db.sh -o ./db                                  # HLA only
#   ./build_db.sh -o ./db -g HLA,KIR,IG_TR
#   ./build_db.sh -o ./db -g CYP -p ./pharmvar-6.2
#
# The db folder is then passed to run_specimmune.sh with -d.

set -euo pipefail

IMAGE="specimmune:v1.0.0"
DB_DIR=""
FAMILIES="HLA"
PHARMVAR=""
FORCE=false

usage() {
  cat <<EOF
Usage: $(basename "$0") -o <db_dir> [options]

Required:
  -o  Database folder to create/fill (later passed to run_specimmune.sh -d)

Options:
  -g  Families to build, comma-separated: HLA,KIR,IG_TR,CYP  [default: ${FAMILIES}]
  -p  Unzipped PharmVar download folder (required for CYP)
  -i  Docker image                                           [default: ${IMAGE}]
  -f  Rebuild families whose folder already exists in the db
  -h  Show this help
EOF
  exit "${1:-0}"
}

while getopts ":o:g:p:i:fh" opt; do
  case "$opt" in
    o) DB_DIR="$OPTARG" ;;
    g) FAMILIES="$OPTARG" ;;
    p) PHARMVAR="$OPTARG" ;;
    i) IMAGE="$OPTARG" ;;
    f) FORCE=true ;;
    h) usage 0 ;;
    :) echo "Option -$OPTARG needs a value" >&2; usage 1 ;;
    *) echo "Unknown option -$OPTARG" >&2; usage 1 ;;
  esac
done

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "$DB_DIR" ]] || usage 1
command -v docker >/dev/null || die "docker is not installed or not on PATH"
docker image inspect "$IMAGE" >/dev/null 2>&1 || die "Docker image not found: $IMAGE (see README: Build the image)"

IFS=',' read -r -a FAMILY_LIST <<< "$FAMILIES"
for fam in "${FAMILY_LIST[@]}"; do
  case "$fam" in
    HLA|KIR|IG_TR) ;;
    CYP)
      [[ -n "$PHARMVAR" ]] || die "CYP needs the unzipped PharmVar folder (-p), from https://www.pharmvar.org/download"
      [[ -d "$PHARMVAR" ]] || die "PharmVar folder not found: $PHARMVAR"
      PHARMVAR=$(realpath "$PHARMVAR") ;;
    *) die "Unknown family '$fam' (use HLA, KIR, IG_TR, CYP)" ;;
  esac
done

mkdir -p "$DB_DIR"
DB_DIR=$(realpath "$DB_DIR")
OWNER="$(id -u):$(id -g)"

failed=()
for fam in "${FAMILY_LIST[@]}"; do
  if [[ -d "$DB_DIR/$fam" ]] && [[ -n "$(ls -A "$DB_DIR/$fam" 2>/dev/null)" ]] && ! $FORCE; then
    log "$fam | $DB_DIR/$fam already exists, skipping (use -f to rebuild)"
    continue
  fi

  mounts=(-v "$DB_DIR:/db")
  cyp_cmd=""
  cyp_arg=""
  if [[ "$fam" == "CYP" ]]; then
    mounts+=(-v "$PHARMVAR:/pharmvar:ro")
    # Merge every PharmVar FASTA into one file, as in the SpecImmune README
    cyp_cmd="find /pharmvar -type f -name '*.fasta' -exec cat {} + > /tmp/CYP.all.fasta
             [ -s /tmp/CYP.all.fasta ] || { echo 'No *.fasta files found in the PharmVar folder' >&2; exit 2; }
             echo \"Merged \$(grep -c '^>' /tmp/CYP.all.fasta) PharmVar sequences\""
    cyp_arg="--CYP_fa /tmp/CYP.all.fasta"
  fi

  # make_db.py may write next to the SpecImmune code, so run as root and hand
  # the database back to the current user afterwards
  cmd="
    set -euo pipefail
    cd /SpecImmune
    $cyp_cmd
    python3 scripts/make_db.py -o /db -i $fam $cyp_arg
    chown -R $OWNER /db
  "

  log "$fam | building into $DB_DIR"
  if docker run --rm "${mounts[@]}" "$IMAGE" bash -c "$cmd"; then
    log "$fam | done"
  else
    log "$fam | FAILED"
    failed+=("$fam")
  fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
  die "Failed: ${failed[*]}"
fi
log "Database ready: $DB_DIR  (use: run_specimmune.sh -d $DB_DIR)"
