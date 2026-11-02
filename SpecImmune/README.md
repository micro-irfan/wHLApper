# SpecImmune

Batch runner for [SpecImmune](https://github.com/deepomicslab/SpecImmune) typing of Nanopore BAMs (HLA, KIR, CYP, IG/TR) in Docker. It is separate from the Nextflow pipeline.

## Build the image

The Dockerfile copies the SpecImmune source code into the image. Build it from a SpecImmune clone that has this `docker/` folder copied in:

```bash
git clone https://github.com/deepomicslab/SpecImmune.git
cp -r docker/ SpecImmune/docker/
cd SpecImmune
docker build -f docker/Dockerfile -t specimmune:v1.0.0 .
```

The database (`-d`) is built separately with SpecImmune's `scripts/make_db.py`. See the SpecImmune README.

## Usage

```bash
./run_specimmune.sh -s samplesheet.csv -o results -d ./db \
    -g HLA,KIR,CYP -r ./ref/hg38_no_alt.fa -m ./ref/hg38_no_alt.mmi -t 16

# Long runs
nohup bash -c './run_specimmune.sh -s samplesheet.csv -o results -d ./db -g HLA \
    -r ./ref/hg38_no_alt.fa -t 16' > SpecImmune.log 2>&1 &
```

| option | |
|---|---|
| `-s` | samplesheet (required) |
| `-o` | output folder (required) |
| `-d` | SpecImmune database folder (required) |
| `-g` | gene families, comma-separated: `HLA`, `KIR`, `CYP`, `IG_TR` (default `HLA`) |
| `-r` | no-alt hg38 FASTA. Required for `CYP`, `IG_TR` and any unaligned BAM |
| `-m` | prebuilt minimap2 index of `-r` (`minimap2 -x map-ont -d hg38_no_alt.mmi hg38_no_alt.fa`) |
| `-c` | reference FASTA for CRAM inputs |
| `-t` | threads per sample (default 8) |
| `-i` | Docker image (default `specimmune:v1.0.0`) |
| `-T` | targeted/amplicon data: use all reads, skip alignment and region extraction |
| `-x` | extra arguments for `main.py`, quoted |

## Samplesheet

A CSV or TSV with `sample_id,bam`. The header row is optional:

```
sample_id,bam
HG00276,/data/pgx/cohort_1/HG00276/barcode21.bam
```

## Steps

For each sample:

1. **Align if needed.** If the BAM header has no `@SQ` lines, the reads are aligned to `-r` with `minimap2 -x map-ont`. The result is written next to the original as `<name>.hg38.bam` (+ `.bai`), so that folder must be writable. If that file already exists, it is reused. BAMs that are already aligned skip this step.
2. **Extract reads** for each gene family with `scripts/ExtractReads.sh`. With `-T`, the whole BAM is converted to FASTQ instead.
3. **Type** with `scripts/main.py -y nanopore`.

## Outputs

```
results/
  run_summary.tsv              sample, gene, OK / FAILED / ALIGN_FAILED / MISSING_BAM, log
  <sample>/reads_<gene>/       extracted FASTQ
  <sample>/typing/             SpecImmune results
  <sample>/logs/               <sample>.align.log, <sample>.<gene>.log
```
