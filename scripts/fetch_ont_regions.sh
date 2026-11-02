#!/usr/bin/env bash
#
# Download ONT open data CRAMs and extract per-target BAMs for typing:
#   <sample_id>.HLA.bam     chr6 (+ chr6 alt / HLA-* contigs) + unmapped reads
#   <sample_id>.KIR.bam     LRC / KIR region (chr19)
#   <sample_id>.IGTCR.bam   IGH, IGK, IGL, TRA/TRD, TRB, TRG
#   <sample_id>.CYP2D6.bam  CYP2D6 / CYP2D7 / CYP2D8P
# The CRAM is deleted once all of its BAMs are written (unless -k). If extraction fails the
# CRAM is kept in <work_dir>/<sample_id>/ and a rerun picks up from there without re-downloading.
#
# Samplesheet (CSV, header optional): sample_id,flowcell_id,s3_path
#   HG001_1,PAW79146,s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/hac/HG001/PAW79146/output/
# 
# nohup bash -c './fetch_ont_regions.sh -i samplesheet.csv -o giab_open_data.2025 -t 8 \
#   -r /data/irfan/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' > Download.test.log 2>&1 &
# 
# Requires: aws CLI, and docker (samtools from ontresearch/wf-common) or samtools on PATH with -L

set -euo pipefail

CRAM_NAME="SAMPLE.haplotagged.cram"
SAMTOOLS_IMAGE="ontresearch/wf-common:sha72f3517dd994984e0e2da0b97cb3f23f8540be4b"

# GRCh38 regions per target. Add '*' to a target to include unmapped reads.
declare -A TARGETS=(
    [KIR]="chr19:54000000-56000000"
    [IGTCR]="chr14:105586437-106879844 chr2:88857361-90235368 chr22:22026076-22922913 chr14:21621904-22552132 chr7:142299011-142813287 chr7:38240024-38368055"
    [CYP2D6]="chr22:42000000-42300000"
)
# HLA regions are built from the CRAM header (chr6, chr6_*_alt, HLA-*) plus '*' - see hla_regions()
TARGET_NAMES=(HLA KIR IGTCR CYP2D6)

usage() {
    cat <<EOF
Usage: $(basename "$0") -i samplesheet.csv -o out_dir [-r ref.fa] [-t threads] [-w work_dir] [-k] [-I image | -L]

  -i  Samplesheet CSV: sample_id,flowcell_id,s3_path
  -o  Output folder for <sample_id>.<target>.bam
  -r  Reference FASTA used to decode the CRAMs (default: fetch from EBI by MD5, cached in work_dir)
  -t  Threads for samtools [4]
  -w  Download folder [out_dir/work]
  -k  Keep downloaded CRAMs
  -I  Docker image providing samtools [$SAMTOOLS_IMAGE]
  -L  Use samtools from PATH instead of Docker
EOF
    exit "${1:-0}"
}

log() { echo "[$(date '+%F %T')] $*" >&2; }

SHEET="" OUTDIR="" REF="" THREADS=4 WORKDIR="" KEEP=0 LOCAL_SAMTOOLS=0
while getopts "i:o:r:t:w:kI:Lh" opt; do
    case $opt in
        i) SHEET=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
        r) REF=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        w) WORKDIR=$OPTARG ;;
        k) KEEP=1 ;;
        I) SAMTOOLS_IMAGE=$OPTARG ;;
        L) LOCAL_SAMTOOLS=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
[[ -n $SHEET && -n $OUTDIR ]] || usage 1
[[ -f $SHEET ]] || { log "Samplesheet not found: $SHEET"; exit 1; }
required=(aws)
if [[ $LOCAL_SAMTOOLS -eq 1 ]]; then required+=(samtools); else required+=(docker); fi
for tool in "${required[@]}"; do
    command -v "$tool" >/dev/null || { log "$tool not found on PATH"; exit 1; }
done

WORKDIR=${WORKDIR:-$OUTDIR/work}
mkdir -p "$OUTDIR" "$WORKDIR"
# Absolute paths so they are identical inside the samtools container
OUTDIR=$(realpath "$OUTDIR")
WORKDIR=$(realpath "$WORKDIR")

REF_ARG=()
if [[ -n $REF ]]; then
    [[ -f $REF ]] || { log "Reference not found: $REF"; exit 1; }
    REF=$(realpath "$REF")
    REF_ARG=(-T "$REF")
else
    # Reference blocks are downloaded once by MD5 and reused for every sample
    export REF_CACHE="$WORKDIR/ref_cache/%2s/%2s/%s"
fi

# samtools from the wf-common image, with the work, output and reference folders
# mounted at the same paths so file arguments work unchanged. -i passes piped stdin through.
if [[ $LOCAL_SAMTOOLS -eq 0 ]]; then
    DOCKER_ARGS=(--rm -i --user "$(id -u):$(id -g)" -v "$WORKDIR:$WORKDIR" -v "$OUTDIR:$OUTDIR")
    [[ -n $REF ]] && DOCKER_ARGS+=(-v "$(dirname "$REF"):$(dirname "$REF"):ro")
    [[ -n ${REF_CACHE:-} ]] && DOCKER_ARGS+=(-e "REF_CACHE=$REF_CACHE")
    log "Using samtools from $SAMTOOLS_IMAGE"
    docker image inspect "$SAMTOOLS_IMAGE" >/dev/null 2>&1 || docker pull "$SAMTOOLS_IMAGE" >&2
    samtools() { docker run "${DOCKER_ARGS[@]}" "$SAMTOOLS_IMAGE" samtools "$@"; }
fi
samtools --version | head -n1 >&2

# chr6, chr6 alt and HLA-* contigs present in the header, plus unmapped reads
hla_regions() {
    local contigs
    contigs=$(samtools view -H "$1" \
        | awk -F'\t' '$1=="@SQ" { for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) { sn = substr($i, 4);
              if (sn == "chr6" || sn == "6" || sn ~ /^chr6_.*_alt$/ || sn ~ /^HLA-/) print sn } }' \
        | tr '\n' ' ')
    [[ -n $contigs ]] || { log "No chr6/HLA contigs found in $1 header"; return 1; }
    echo "$contigs *"
}

done_already() {
    local id=$1 t
    for t in "${TARGET_NAMES[@]}"; do
        [[ -s $OUTDIR/$id.$t.bam && -s $OUTDIR/$id.$t.bam.bai ]] || return 1
    done
}

extract() {
    local cram=$1 out=$2 regions
    shift 2
    regions=("$@")
    # Sort because regions may span several contigs in any order; write to .tmp so reruns never see partial files
    samtools view -@ "$THREADS" -u ${REF_ARG[@]+"${REF_ARG[@]}"} "$cram" "${regions[@]}" \
        | samtools sort -@ "$THREADS" -o "$out.tmp.bam" -
    samtools index -@ "$THREADS" "$out.tmp.bam"
    mv "$out.tmp.bam" "$out"
    mv "$out.tmp.bam.bai" "$out.bai"
}

process_sample() {
    local id=$1 flowcell=$2 path=$3
    local dir="$WORKDIR/$id"
    local cram="$dir/$CRAM_NAME"

    while [[ $path == */ ]]; do path=${path%/}; done

    # Download the CRAM unless it is already in the work folder. It is fetched to .part and
    # renamed when complete, so an interrupted download is never mistaken for a finished one.
    mkdir -p "$dir"
    if [[ -s $cram ]]; then
        log "$id: $CRAM_NAME already downloaded, reusing"
    else
        log "$id ($flowcell): downloading $path/$CRAM_NAME"
        aws s3 cp --no-sign-request --only-show-errors "$path/$CRAM_NAME" "$cram.part"
        mv "$cram.part" "$cram"
    fi

    # The published .crai can't be used, so index the CRAM ourselves. The marker tells our
    # index apart from a .crai downloaded by an older version of this script.
    local idx_marker="$dir/.crai_generated"
    if [[ ! -s $cram.crai || ! -e $idx_marker ]]; then
        log "$id: indexing $CRAM_NAME"
        rm -f "$cram.crai" "$idx_marker"
        samtools index -@ "$THREADS" "$cram" "$cram.crai.part"
        mv "$cram.crai.part" "$cram.crai"
        touch "$idx_marker"
    else
        log "$id: $CRAM_NAME.crai already built, reusing"
    fi

    local t regions
    for t in "${TARGET_NAMES[@]}"; do
        if [[ $t == HLA ]]; then
            regions=$(hla_regions "$cram")
        else
            regions=${TARGETS[$t]}
        fi
        log "$id: extracting $t"
        # shellcheck disable=SC2086  # word-split the region list; '*' must not glob
        (set -f; extract "$cram" "$OUTDIR/$id.$t.bam" $regions)
    done

    # Only reached when every target succeeded (set -e); on failure the CRAM stays for a rerun
    if [[ $KEEP -eq 0 ]]; then
        rm -rf "$dir"
    fi
    log "$id: done"
}

shopt -s extglob
failed=()
while IFS=, read -r id flowcell path _ || [[ -n ${id:-} ]]; do
    id=${id//$'\r'/} flowcell=${flowcell//$'\r'/} path=${path//$'\r'/}
    id=$(echo "$id" | xargs) flowcell=$(echo "$flowcell" | xargs) path=$(echo "$path" | xargs)
    [[ -z $id || $id == sample_id || $id == \#* ]] && continue
    if [[ -z $path ]]; then
        log "$id: missing s3 path, skipping"
        failed+=("$id")
        continue
    fi
    if done_already "$id"; then
        log "$id: outputs exist, skipping"
        continue
    fi
    # Not run as an `if` condition so that set -e still applies inside; stdin kept off the samplesheet
    set +e
    (set -e; process_sample "$id" "$flowcell" "$path") < /dev/null
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
        log "$id: FAILED (exit $rc) - any downloaded files are kept in $WORKDIR/$id, rerun to retry"
        failed+=("$id")
    fi
done < "$SHEET"

if [[ ${#failed[@]} -gt 0 ]]; then
    log "Failed samples: ${failed[*]}"
    exit 1
fi
log "All samples done -> $OUTDIR"
