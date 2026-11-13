#!/usr/bin/env bash
#
# specimmune_wrapper.sh - batch CYP2D6 genotyping with SpecImmune, using Docker
#
# Runs SpecImmune/run_specimmune.sh with the CYP gene family for every sample in
# the sample sheet (same sheet format as chinook_wrapper.sh / run_specimmune.sh):
#   sample_id,bam
#   S001,/data/ont/S001.hg38.bam
#   S002,/data/ont/S002.hg38.bam
# then pulls the CYP2D6 rows of each sample's CYP result into one table,
# <outdir>/cyp2d6_calls.tsv (one row per haplotype).
#
# Unaligned BAMs are aligned to -r first (by run_specimmune.sh). The SpecImmune
# database (-d) must include CYP: ../../SpecImmune/build_db.sh -o ./db -g CYP -p <pharmvar>
#
# Example:
#   ./specimmune_wrapper.sh -s samples.csv -o cyp2d6_specimmune -d ./db -r ./ref/hg38_no_alt.fa -t 16
#   nohup ./specimmune_wrapper.sh -s samples.csv -o cyp2d6_specimmune -d ./db \
#       -r ./ref/hg38_no_alt.fa -t 16 > specimmune_cyp2d6.log 2>&1 &


set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENE="CYP2D6"

# ---------------------------------------------------------------- defaults
RUNNER="$SCRIPT_DIR/../../SpecImmune/run_specimmune.sh"
SAMPLESHEET=""
OUTDIR=""
DB_DIR=""
HG38_NOALT=""
PASS_ARGS=()       # forwarded to run_specimmune.sh

usage() {
  cat <<EOF
Usage: $(basename "$0") -s <samplesheet> -o <outdir> -d <db_dir> -r <hg38_no_alt.fa> [options]

Required:
  -s  Sample sheet (CSV/TSV): sample_id, path to BAM/CRAM
  -o  Output directory
  -d  SpecImmune database folder (must include CYP)
  -r  No-alt hg38 FASTA (SpecImmune needs it for CYP, and to align unaligned BAMs)

Options:
  -R  Path to run_specimmune.sh                    [default: ${RUNNER}]

Options (passed to run_specimmune.sh):
  -m  Prebuilt minimap2 index (.mmi) of the -r FASTA
  -c  Reference FASTA used for CRAM inputs
  -t  Threads per sample
  -i  Docker image
  -T  Targeted/amplicon data: skip read extraction and use all reads
  -x  Extra arguments for SpecImmune main.py, quoted
  -h  Show this help
EOF
  exit "${1:-0}"
}

while getopts ":s:o:d:r:R:m:c:t:i:x:Th" opt; do
  case "$opt" in
    R) RUNNER="$OPTARG" ;;
    s) SAMPLESHEET="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    d) DB_DIR="$OPTARG" ;;
    r) HG38_NOALT="$OPTARG" ;;
    m|c|t|i|x) PASS_ARGS+=("-$opt" "$OPTARG") ;;
    T) PASS_ARGS+=(-T) ;;
    h) usage 0 ;;
    :) echo "Option -$OPTARG needs a value" >&2; usage 1 ;;
    *) echo "Unknown option -$OPTARG" >&2; usage 1 ;;
  esac
done

# ---------------------------------------------------------------- checks
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[[ -n "$SAMPLESHEET" && -n "$OUTDIR" && -n "$DB_DIR" && -n "$HG38_NOALT" ]] || usage 1
[[ -f "$RUNNER" ]] || die "run_specimmune.sh not found: $RUNNER (set it with -R)"

# ---------------------------------------------------------------- typing
# run_specimmune.sh loops over the sample sheet and does its own checks
bash "$RUNNER" -s "$SAMPLESHEET" -o "$OUTDIR" -d "$DB_DIR" -r "$HG38_NOALT" -g CYP ${PASS_ARGS[@]+"${PASS_ARGS[@]}"}
OUTDIR=$(realpath "$OUTDIR")

# ---------------------------------------------------------------- collect CYP2D6
# CYP results follow the HLA layout:
#   <outdir>/<sample>/typing/<sample>/<sample>.CYP.final.type.result.formatted.txt
# with columns Locus, Chromosome, Genotype, ..., Reads_num, ..., One_guess.
# Allele = One_guess, or the first Genotype candidate when that is NA.
CALLS="$OUTDIR/cyp2d6_calls.tsv"
printf "sample_id\thaplotype\tallele\tcandidates\treads_num\tresult\n" > "$CALLS"

n=0
while IFS=$'\t' read -r sample gene status _log; do
  [[ "$sample" == "sample_id" || "$gene" != "CYP" || "$status" != "OK" ]] && continue
  rdir="$OUTDIR/$sample/typing/$sample"
  res="$rdir/$sample.CYP.final.type.result.formatted.txt"
  [[ -f "$res" ]] || res=$(ls "$rdir"/*CYP*final*formatted*.txt 2>/dev/null | head -n1 || true)
  if [[ -z "$res" || ! -f "$res" ]]; then
    log "$sample | WARNING: no CYP result file in $rdir"
    continue
  fi

  rows=$(awk -v s="$sample" -v g="$GENE" -v f="$res" 'BEGIN{FS=OFS="\t"}
    /^#/ || NF == 0 { next }
    !hdr { for (i = 1; i <= NF; i++) col[$i] = i; hdr = 1; next }
    $col["Locus"] == g {
      geno  = $col["Genotype"]; split(geno, cand, ";")
      guess = ("One_guess" in col) ? $col["One_guess"] : "NA"
      allele = (guess != "" && guess != "NA" && guess != "-") ? guess : (cand[1] == "NA" ? "" : cand[1])
      reads = ("Reads_num" in col) ? $col["Reads_num"] : ""
      print s, $col["Chromosome"], allele, geno, reads, f
    }' "$res")
  if [[ -z "$rows" ]]; then
    log "$sample | WARNING: no $GENE rows in $res"
    continue
  fi
  printf "%s\n" "$rows" >> "$CALLS"
  n=$((n + 1))
done < "$OUTDIR/run_summary.tsv"

log "$GENE calls for $n sample(s): $CALLS"
column -t -s $'\t' <(cut -f1-3,5 "$CALLS") 2>/dev/null || cut -f1-3,5 "$CALLS"
