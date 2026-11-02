#!/usr/bin/env bash
#
# run_specimmune.sh - batch SpecImmune typing of Nanopore BAMs with Docker
#
# Sample sheet: CSV or TSV, one sample per line, with or without a header:
#   sample_id,bam
#   NA12878,/data/ont/NA12878.hg38.bam
#   HG002,/data/ont/HG002.hg38.bam
#
# For each sample it first checks whether the BAM is aligned (@SQ lines in the header).
# Unaligned BAMs are aligned to the no-alt hg38 reference with minimap2 (map-ont) and
# written next to the original as <name>.hg38.bam. Then for each gene family it:
#   1. extracts gene-region reads from the BAM (scripts/ExtractReads.sh)
#      - or, with --targeted, converts the whole BAM to FASTQ (amplicon data)
#   2. runs scripts/main.py with Nanopore settings
#
# Example:
#   ./run_specimmune.sh -s samples.csv -o results -d ./db \
#       -g HLA,KIR,CYP -r ./ref/hg38_no_alt.fa -t 16
#
#   nohup bash -c './run_specimmune.sh -s samplesheet.csv -o results -d ./db \
#       -g HLA -r ./ref/hg38_no_alt.fa -t 16' > SpecImmune.log 2>&1 &

set -euo pipefail

# ---------------------------------------------------------------- defaults
IMAGE="specimmune:v1.0.0"
GENES="HLA"
THREADS=8
DB_DIR=""
HG38_NOALT=""      # needed for CYP and IG_TR, and to align unaligned BAMs
HG38_MMI=""        # optional prebuilt minimap2 index of HG38_NOALT
CRAM_REF=""        # only needed if inputs are CRAM
TARGETED=false     # true = amplicon/targeted BAM, skip region extraction
EXTRA_ARGS=""      # passed straight to main.py
SAMPLESHEET=""
OUTDIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") -s <samplesheet> -o <outdir> -d <db_dir> [options]

Required:
  -s  Sample sheet (CSV/TSV): sample_id, path to BAM/CRAM
  -o  Output directory
  -d  SpecImmune database folder (built with scripts/make_db.py)

Options:
  -g  Gene families, comma-separated: HLA,KIR,CYP,IG_TR   [default: ${GENES}]
  -r  No-alt hg38 FASTA (required for CYP and IG_TR, and for unaligned BAMs)
  -m  Prebuilt minimap2 index (.mmi) of the -r FASTA, used when aligning unaligned
      BAMs (otherwise minimap2 indexes the FASTA for every sample)
  -c  Reference FASTA used for CRAM inputs
  -t  Threads per sample                                     [default: ${THREADS}]
  -i  Docker image                                           [default: ${IMAGE}]
  -T  Targeted/amplicon data: skip read extraction and use all reads
  -x  Extra arguments for main.py, quoted (e.g. "--min_identity 0.8")
  -h  Show this help
EOF
  exit "${1:-0}"
}

while getopts ":s:o:d:g:r:m:c:t:i:x:Th" opt; do
  case "$opt" in
    s) SAMPLESHEET="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    d) DB_DIR="$OPTARG" ;;
    g) GENES="$OPTARG" ;;
    r) HG38_NOALT="$OPTARG" ;;
    m) HG38_MMI="$OPTARG" ;;
    c) CRAM_REF="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    i) IMAGE="$OPTARG" ;;
    x) EXTRA_ARGS="$OPTARG" ;;
    T) TARGETED=true ;;
    h) usage 0 ;;
    :) echo "Option -$OPTARG needs a value" >&2; usage 1 ;;
    *) echo "Unknown option -$OPTARG" >&2; usage 1 ;;
  esac
done

# ---------------------------------------------------------------- checks
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "$SAMPLESHEET" && -n "$OUTDIR" && -n "$DB_DIR" ]] || usage 1
[[ -f "$SAMPLESHEET" ]] || die "Sample sheet not found: $SAMPLESHEET"
[[ -d "$DB_DIR" ]]      || die "Database folder not found: $DB_DIR"
command -v docker >/dev/null || die "docker is not installed or not on PATH"
docker image inspect "$IMAGE" >/dev/null 2>&1 || die "Docker image not found: $IMAGE"

IFS=',' read -r -a GENE_LIST <<< "$GENES"
for g in "${GENE_LIST[@]}"; do
  case "$g" in
    HLA|KIR) ;;
    CYP|IG_TR)
      [[ -n "$HG38_NOALT" ]] || die "$g typing needs the no-alt hg38 reference (-r)"
      [[ -f "$HG38_NOALT" ]] || die "Reference not found: $HG38_NOALT" ;;
    *) die "Unknown gene family '$g' (use HLA, KIR, CYP, IG_TR)" ;;
  esac
done

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
DB_DIR=$(realpath "$DB_DIR")
[[ -n "$HG38_NOALT" ]] && HG38_NOALT=$(realpath "$HG38_NOALT")
[[ -n "$CRAM_REF" ]]   && CRAM_REF=$(realpath "$CRAM_REF")
if [[ -n "$HG38_MMI" ]]; then
  [[ -f "$HG38_MMI" ]] || die "minimap2 index not found: $HG38_MMI"
  HG38_MMI=$(realpath "$HG38_MMI")
fi

# Run containers as the current user so output files aren't owned by root
DOCKER_USER="$(id -u):$(id -g)"

# ---------------------------------------------------------------- per-sample
SUMMARY="$OUTDIR/run_summary.tsv"
printf "sample_id\tgene\tstatus\tlog\n" > "$SUMMARY"

# Align an unaligned BAM to the no-alt hg38 reference, writing
# <bam name>.hg38.bam next to the original. Sets ALIGNED_BAM to the BAM to use
# downstream (the input itself if it is already aligned).
align_if_needed() {
  local sample="$1" bam="$2"
  local bam_dir bam_name
  bam_dir=$(dirname "$bam")
  bam_name=$(basename "$bam")
  local out_name="${bam_name%.*}.hg38.bam"
  local out="$bam_dir/$out_name"
  local logf="$OUTDIR/$sample/logs/${sample}.align.log"
  ALIGNED_BAM="$bam"

  local mounts=(-v "$bam_dir:/in:ro")
  local cram_opt=""
  if [[ -n "$CRAM_REF" ]]; then
    mounts+=(-v "$(dirname "$CRAM_REF"):/cramref:ro")
    cram_opt="--reference /cramref/$(basename "$CRAM_REF")"
  fi

  # Step 1: is the BAM already aligned? Aligned BAMs have @SQ lines in the header
  local n_sq
  n_sq=$(docker run --rm --user "$DOCKER_USER" "${mounts[@]}" "$IMAGE" \
           bash -c "samtools view -H $cram_opt $(printf '%q' "/in/$bam_name") | grep -c '^@SQ' || true")
  if [[ "${n_sq:-0}" -gt 0 ]]; then
    log "$sample | align | BAM is already aligned ($n_sq contigs), skipping minimap2"
    return 0
  fi

  # Step 2: reuse a previous alignment
  if [[ -s "$out" && -s "$out.bai" ]]; then
    log "$sample | align | BAM is unaligned, reusing $out"
    ALIGNED_BAM="$out"
    return 0
  fi

  # Step 3: align with minimap2, output next to the original BAM
  [[ -n "$HG38_NOALT" ]] || { log "$sample | align | unaligned BAM needs the no-alt hg38 reference (-r)"; return 1; }
  [[ -w "$bam_dir" ]] || { log "$sample | align | cannot write to $bam_dir"; return 1; }
  mkdir -p "$OUTDIR/$sample/logs"
  mounts=(-v "$bam_dir:/in" "${mounts[@]:2}" -v "$(dirname "$HG38_NOALT"):/ref:ro")
  local index="/ref/$(basename "$HG38_NOALT")"
  if [[ -n "$HG38_MMI" ]]; then
    mounts+=(-v "$(dirname "$HG38_MMI"):/mmi:ro")
    index="/mmi/$(basename "$HG38_MMI")"
  fi

  local sort_threads=$(( THREADS / 4 > 0 ? (THREADS / 4 < 4 ? THREADS / 4 : 4) : 1 ))
  local map_threads=$(( THREADS - sort_threads > 0 ? THREADS - sort_threads : 1 ))
  local q_out q_tmp
  q_out=$(printf '%q' "/in/$out_name")
  q_tmp=$(printf '%q' "/in/${out_name%.bam}.tmp.bam")
  # -T '*' carries the BAM tags (e.g. MM/ML) into the FASTQ comment, -y copies them back.
  # Written to a .tmp name first so a failed run is never mistaken for a finished one.
  local cmd="
    set -euo pipefail
    samtools fastq -@ 2 -T '*' $cram_opt $(printf '%q' "/in/$bam_name") \
      | minimap2 -t $map_threads -ayL -x map-ont $(printf '%q' "$index") - \
      | samtools sort -@ $sort_threads -o $q_tmp -
    samtools index -@ $THREADS $q_tmp
    mv $q_tmp $q_out
    mv $q_tmp.bai $q_out.bai
  "

  log "$sample | align | BAM is unaligned, aligning to $(basename "$HG38_NOALT") with minimap2 -> $out"
  if docker run --rm --user "$DOCKER_USER" "${mounts[@]}" "$IMAGE" \
       bash -c "$cmd" > "$logf" 2>&1; then
    log "$sample | align | done"
    ALIGNED_BAM="$out"
  else
    log "$sample | align | FAILED (see $logf)"
    return 1
  fi
}

run_one() {
  local sample="$1" bam="$2" gene="$3"
  local sdir="$OUTDIR/$sample"
  local logf="$sdir/logs/${sample}.${gene}.log"
  mkdir -p "$sdir/logs" "$sdir/reads_${gene}"

  local bam_dir bam_name
  bam_dir=$(dirname "$bam")
  bam_name=$(basename "$bam")

  # Mounts: BAM folder read-only, outputs, db, and references if given
  local mounts=(
    -v "$bam_dir:/in:ro"
    -v "$OUTDIR:/out"
    -v "$DB_DIR:/db:ro"
  )
  local ref_arg=""
  if [[ -n "$HG38_NOALT" ]]; then
    mounts+=(-v "$(dirname "$HG38_NOALT"):/ref:ro")
    ref_arg="--hg38 /ref/$(basename "$HG38_NOALT")"
  fi
  local cram_arg=""
  if [[ -n "$CRAM_REF" ]]; then
    mounts+=(-v "$(dirname "$CRAM_REF"):/cramref:ro")
    cram_arg="-r /cramref/$(basename "$CRAM_REF")"
  fi

  local reads_dir="/out/$sample/reads_${gene}"
  # Shell-quote paths so filenames with spaces/special characters survive
  local q_in q_reads
  q_in=$(printf '%q' "/in/$bam_name")
  q_reads=$(printf '%q' "$reads_dir")

  # Step 1: get the reads as FASTQ
  local fq_cmd
  if $TARGETED; then
    fq_cmd="samtools fastq -@ $THREADS ${CRAM_REF:+--reference /cramref/$(basename "$CRAM_REF")} \
              $q_in | gzip > $q_reads/${sample}.${gene}.fq.gz"
  else
    fq_cmd="bash /SpecImmune/scripts/ExtractReads.sh -s $sample -i $q_in \
              -g $gene -o $q_reads $cram_arg"
  fi

  # Step 2: type with Nanopore settings. ExtractReads.sh names its output
  # itself, so pick up whichever FASTQ it wrote.
  local cmd="
    set -euo pipefail
    cd /SpecImmune
    $fq_cmd
    fq=\$(ls -t $q_reads/*.f*q* 2>/dev/null | head -n1 || true)
    if [ -z \"\$fq\" ] || [ ! -s \"\$fq\" ]; then
      echo 'No reads extracted for $gene - check the BAM is aligned to hg38 with chr-style names and is indexed' >&2
      exit 2
    fi
    echo \"Typing with reads: \$fq\"
    python3 scripts/main.py \
      -r \"\$fq\" \
      -n $sample \
      -o /out/$sample/typing \
      -i $gene \
      -y nanopore \
      -j $THREADS \
      --db /db \
      --align_method_1 minimap2 \
      $ref_arg $EXTRA_ARGS
  "

  log "$sample | $gene | starting"
  if docker run --rm --user "$DOCKER_USER" "${mounts[@]}" "$IMAGE" \
       bash -c "$cmd" > "$logf" 2>&1; then
    log "$sample | $gene | done"
    printf "%s\t%s\tOK\t%s\n" "$sample" "$gene" "$logf" >> "$SUMMARY"
  else
    log "$sample | $gene | FAILED (see $logf)"
    printf "%s\t%s\tFAILED\t%s\n" "$sample" "$gene" "$logf" >> "$SUMMARY"
  fi
}

# ---------------------------------------------------------------- main loop
n=0
while IFS=$',\t' read -r sample bam _rest || [[ -n "${sample:-}" ]]; do
  sample="${sample//$'\r'/}"; bam="${bam//$'\r'/}"         # handle Windows line endings
  sample="$(echo "$sample" | xargs)"; bam="$(echo "$bam" | xargs)"
  [[ -z "$sample" || "$sample" == \#* ]] && continue      # skip blanks/comments
  lc_sample=$(echo "$sample" | tr '[:upper:]' '[:lower:]')   # works on macOS bash 3.2 too
  [[ "$lc_sample" == "sample_id" || "$lc_sample" == "sample" ]] && continue   # header

  if [[ ! -f "$bam" ]]; then
    log "$sample | SKIPPED: file not found: $bam"
    printf "%s\t-\tMISSING_BAM\t-\n" "$sample" >> "$SUMMARY"
    continue
  fi
  bam=$(realpath "$bam")

  # Region extraction needs aligned reads; targeted mode uses all reads as-is
  if ! $TARGETED; then
    if ! align_if_needed "$sample" "$bam"; then
      printf "%s\t-\tALIGN_FAILED\t%s\n" "$sample" "$OUTDIR/$sample/logs/${sample}.align.log" >> "$SUMMARY"
      continue
    fi
    bam="$ALIGNED_BAM"
  fi

  if ! $TARGETED && [[ ! -f "$bam.bai" && ! -f "${bam%.bam}.bai" && ! -f "$bam.crai" ]]; then
    log "$sample | WARNING: no index found next to $bam - region extraction needs one (samtools index)"
  fi

  for gene in "${GENE_LIST[@]}"; do
    run_one "$sample" "$bam" "$gene"
  done
  n=$((n + 1))
done < "$SAMPLESHEET"

log "Finished $n sample(s). Summary: $SUMMARY"
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
