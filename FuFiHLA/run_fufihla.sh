#!/usr/bin/env bash
#
# run_fufihla.sh - batch FuFiHLA typing (HLA-A, -B, -C, -DQA1, -DQB1, -DRB1)
#                  from long-read BAMs, using Docker
#
# Sample sheet: CSV or TSV, one sample per line, with or without a header:
#   sample_id,bam
#   S001,/data/ont/S001.hg38.bam
#   S002,/data/ont/S002.hg38.bam
#
# For each sample it:
#   1. pulls reads overlapping the six HLA genes out of the BAM (hg38 regions
#      from the FuFiHLA README), dropping secondary/supplementary alignments
#      - or, with -A, uses every read in the BAM (targeted/unaligned data)
#   2. runs fufihla on them
# and finally merges all calls into one table.
#
# Example:
#   ./run_fufihla.sh -s samples.csv -o fufihla_results            # ONT (default)
#   ./run_fufihla.sh -s samples.csv -o fufihla_results -p hifi    # PacBio HiFi

set -euo pipefail

# ---------------------------------------------------------------- defaults
IMAGE="fufihla:v0.2.4"
PLATFORM="ont"         # ont | hifi
THREADS=4              # used by samtools for extraction
REFDIR=""              # optional custom IMGT reference (from fufihla-ref-prep)
HOM_THRESHOLD=""       # optional --hom-threshold
CRAM_REF=""            # reference FASTA for CRAM inputs
ALL_READS=false        # true = skip region extraction, use all reads
DEBUG=false
SAMPLESHEET=""
OUTDIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") -s <samplesheet> -o <outdir> [options]

Required:
  -s  Sample sheet (CSV/TSV): sample_id, path to BAM/CRAM
  -o  Output directory

Options:
  -p  Read type: ont or hifi                       [default: ${PLATFORM}]
  -i  Docker image                                 [default: ${IMAGE}]
  -t  Threads for read extraction                  [default: ${THREADS}]
  -R  FuFiHLA reference folder (from fufihla-ref-prep; default: bundled IMGT set)
  -H  Homozygosity threshold (--hom-threshold)
  -c  Reference FASTA for CRAM inputs
  -A  Use all reads in the BAM (targeted/amplicon or unaligned data)
  -D  Keep intermediate files (--debug)
  -h  Show this help
EOF
  exit "${1:-0}"
}

while getopts ":s:o:p:i:t:R:H:c:ADh" opt; do
  case "$opt" in
    s) SAMPLESHEET="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    p) PLATFORM="$OPTARG" ;;
    i) IMAGE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    R) REFDIR="$OPTARG" ;;
    H) HOM_THRESHOLD="$OPTARG" ;;
    c) CRAM_REF="$OPTARG" ;;
    A) ALL_READS=true ;;
    D) DEBUG=true ;;
    h) usage 0 ;;
    :) echo "Option -$OPTARG needs a value" >&2; usage 1 ;;
    *) echo "Unknown option -$OPTARG" >&2; usage 1 ;;
  esac
done

# ---------------------------------------------------------------- checks
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "$SAMPLESHEET" && -n "$OUTDIR" ]] || usage 1
[[ -f "$SAMPLESHEET" ]] || die "Sample sheet not found: $SAMPLESHEET"
case "$PLATFORM" in ont|hifi) ;; *) die "-p must be 'ont' or 'hifi'" ;; esac
command -v docker >/dev/null || die "docker is not installed or not on PATH"
docker image inspect "$IMAGE" >/dev/null 2>&1 \
  || die "Docker image not found: $IMAGE (build it with: docker build -f Dockerfile.fufihla -t $IMAGE .)"

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
if [[ -n "$REFDIR" ]]; then
  [[ -d "$REFDIR" ]] || die "Reference folder not found: $REFDIR"
  REFDIR=$(realpath "$REFDIR")
fi
if [[ -n "$CRAM_REF" ]]; then
  [[ -f "$CRAM_REF" ]] || die "CRAM reference not found: $CRAM_REF"
  CRAM_REF=$(realpath "$CRAM_REF")
fi

DOCKER_USER="$(id -u):$(id -g)"

# hg38 regions of the six genes FuFiHLA types (from the FuFiHLA README)
BED="$OUTDIR/hla_6genes.hg38.bed"
cat > "$BED" <<'EOF'
chr6	29942254	29945755	HLA-A
chr6	31268254	31272571	HLA-C
chr6	31353362	31357442	HLA-B
chr6	32578769	32589848	HLA-DRB1
chr6	32636717	32643200	HLA-DQA1
chr6	32660031	32667132	HLA-DQB1
EOF

# Optional fufihla arguments
FUFI_ARGS="--$PLATFORM"
[[ -n "$REFDIR" ]]        && FUFI_ARGS+=" --refdir /refdir"
[[ -n "$HOM_THRESHOLD" ]] && FUFI_ARGS+=" --hom-threshold $HOM_THRESHOLD"
$DEBUG                    && FUFI_ARGS+=" --debug"

SUMMARY="$OUTDIR/run_summary.tsv"
printf "sample_id\tstatus\treads_used\tlog\n" > "$SUMMARY"

# ---------------------------------------------------------------- per-sample
run_one() {
  local sample="$1" bam="$2"
  local sdir="$OUTDIR/$sample"
  local logf="$sdir/${sample}.log"
  mkdir -p "$sdir"

  local bam_name q_in
  bam_name=$(basename "$bam")
  q_in=$(printf '%q' "/in/$bam_name")

  local mounts=(
    -v "$(dirname "$bam"):/in:ro"
    -v "$OUTDIR:/out"
  )
  [[ -n "$REFDIR" ]] && mounts+=(-v "$REFDIR:/refdir:ro")
  local ref_opt=""
  if [[ -n "$CRAM_REF" ]]; then
    mounts+=(-v "$(dirname "$CRAM_REF"):/cramref:ro")
    ref_opt="--reference $(printf '%q' "/cramref/$(basename "$CRAM_REF")")"
  fi

  # Step 1: reads -> gzipped FASTA. -F 0x900 drops secondary and supplementary
  # alignments so each read is only passed once.
  local extract
  if $ALL_READS; then
    extract="samtools fasta -@ $THREADS -F 0x900 $ref_opt $q_in"
  else
    extract="samtools view -@ $THREADS -bh -F 0x900 $ref_opt $q_in --region-file /out/hla_6genes.hg38.bed \
               | samtools fasta -"
  fi

  # Step 2: FuFiHLA writes <out> (folder) and <out>.out (allele calls) side by side
  local cmd="
    set -euo pipefail
    cd /out/$sample
    rm -rf $sample $sample.out
    $extract | gzip -c > ${sample}.hla_reads.fa.gz
    n=\$(zcat ${sample}.hla_reads.fa.gz | grep -c '^>' || true)
    echo \"reads_used=\$n\"
    if [ \"\$n\" -eq 0 ]; then
      echo 'No HLA reads extracted - is the BAM aligned to hg38 (chr6 naming) and indexed?' >&2
      exit 2
    fi
    fufihla --fa ${sample}.hla_reads.fa.gz --out $sample $FUFI_ARGS
  "

  log "$sample | starting"
  if docker run --rm --user "$DOCKER_USER" "${mounts[@]}" "$IMAGE" \
       bash -c "$cmd" > "$logf" 2>&1 && [[ -s "$sdir/$sample.out" ]]; then
    local nreads
    nreads=$(grep -m1 '^reads_used=' "$logf" | cut -d= -f2)
    log "$sample | done ($nreads reads)"
    printf "%s\tOK\t%s\t%s\n" "$sample" "$nreads" "$logf" >> "$SUMMARY"
  else
    local nreads
    nreads=$(grep -m1 '^reads_used=' "$logf" 2>/dev/null | cut -d= -f2 || true)
    log "$sample | FAILED (see $logf)"
    printf "%s\tFAILED\t%s\t%s\n" "$sample" "${nreads:--}" "$logf" >> "$SUMMARY"
  fi
}

# ---------------------------------------------------------------- main loop
n=0
while IFS=$',\t' read -r sample bam _rest || [[ -n "${sample:-}" ]]; do
  sample="${sample//$'\r'/}"; bam="${bam//$'\r'/}"
  sample="$(echo "$sample" | xargs)"; bam="$(echo "$bam" | xargs)"
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  lc=$(echo "$sample" | tr '[:upper:]' '[:lower:]')
  [[ "$lc" == "sample_id" || "$lc" == "sample" ]] && continue

  if [[ ! -f "$bam" ]]; then
    log "$sample | SKIPPED: file not found: $bam"
    printf "%s\tMISSING_BAM\t-\t-\n" "$sample" >> "$SUMMARY"
    continue
  fi
  bam=$(realpath "$bam")

  if ! $ALL_READS && [[ ! -f "$bam.bai" && ! -f "${bam%.*}.bai" && ! -f "$bam.crai" && ! -f "$bam.csi" ]]; then
    log "$sample | WARNING: no index next to $bam - region extraction needs one (samtools index)"
  fi

  run_one "$sample" "$bam"
  n=$((n + 1))
done < "$SAMPLESHEET"

# ---------------------------------------------------------------- merge calls
# Column 1 of <sample>.out = allele call, column 2 = consensus ID, last column = cs:Z tag.
# A cs tag that is only a match length (cs:Z::3503) means an exact known allele;
# anything else means the consensus differs from the database (possible novel allele).
ALL="$OUTDIR/all_samples_hla_calls.tsv"
printf "sample_id\tallele\tconsensus_id\texact_match\tcs_tag\n" > "$ALL"
for f in "$OUTDIR"/*/*.out; do
  [[ -s "$f" ]] || continue
  s=$(basename "$(dirname "$f")")
  awk -v s="$s" 'BEGIN{FS=OFS="\t"} NF>0 && $1 !~ /^#/ {
      cs=$NF; exact = (cs ~ /^cs:Z::[0-9]+$/) ? "yes" : "no"
      print s, $1, $2, exact, cs }' "$f" >> "$ALL"
done

log "Finished $n sample(s)."
log "Run summary: $SUMMARY"
log "All calls:   $ALL"
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
