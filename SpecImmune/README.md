# SpecImmune

Batch runner for [SpecImmune](https://github.com/deepomicslab/SpecImmune) typing of Nanopore BAMs (HLA, KIR, CYP, IG/TR) in Docker. It is separate from the Nextflow pipeline.

## Build the image

The Dockerfile copies the SpecImmune source code into the image. Build it from a SpecImmune clone that has this `docker/` folder copied in:

```bash
git clone https://github.com/micro-irfan/wHLApper.git
git clone https://github.com/deepomicslab/SpecImmune.git

# Copy files over from the docker directory into SpecImmune repository 
cp -r wHLApper/SpecImmune/docker/* SpecImmune/docker/

cd SpecImmune
docker build -f docker/Dockerfile -t specimmune:v1.0.0 .
```

## Build the database

`run_specimmune.sh -d` needs a SpecImmune database for every gene family you type. `build_db.sh` (Not tested) builds it inside the image with SpecImmune's `scripts/make_db.py`, running `python scripts/make_db.py -o <db> -i <family>` once per family 

```bash
## Documented
mkdir -p db
docker run --rm   -v "$(pwd)/db:/SpecImmune/db"   -w /SpecImmune   --user "$(id -u):$(id -g)"   specimmune:v1.0.0   python scripts/make_db.py -o ./db -i HLA

# it will lead to an error ~missing hla_gen.fasta: 
cd db/HLA 
unzip hla_gen.fasta.zip 
cd ../../

# Rerun the above command - it should complete with no error afterwards
docker run --rm   -v "$(pwd)/db:/SpecImmune/db"   -w /SpecImmune   --user "$(id -u):$(id -g)"   specimmune:v1.0.0   python scripts/make_db.py -o ./db -i HLA

```

CYP also needs the [PharmVar](https://www.pharmvar.org/download) database. Download the complete set, unzip it, and pass the folder with `-p`. All its `*.fasta` files are merged and passed to `make_db.py --CYP_fa`:

```bash
unzip pharmvar-6.2.zip                       # use your PharmVar version
./build_db.sh -o ./db -g CYP -p ./pharmvar-6.2
```

| option | |
|---|---|
| `-o` | database folder (required). Pass the same folder to `run_specimmune.sh -d` |
| `-g` | families, comma-separated: `HLA`, `KIR`, `IG_TR`, `CYP` (default `HLA`) |
| `-p` | unzipped PharmVar folder (required for `CYP`) |
| `-i` | Docker image (default `specimmune:v1.0.0`) |
| `-f` | rebuild a family whose folder already exists in the database |

Families whose `<db>/<family>` folder already exists are skipped, so you can add a family to an existing database later. The container runs as root, because `make_db.py` may write next to the SpecImmune code, and then hands the database back to your user.

Typing CYP and IG_TR also needs the no-alt hg38 reference at run time (`run_specimmune.sh -r`), not when building the database. SpecImmune provides one at https://doi.org/10.5281/zenodo.14722796.

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

The HLA calls for each sample are in `results/<sample>/typing/<sample>/<sample>.HLA.final.type.result.formatted.txt`. See `sample_output/` for an example.

## Summarise the HLA calls

`parse_specimmune.py` compiles the HLA result files into one CSV with two rows per locus (`<locus>_1`, `<locus>_2`) and one column per sample. It uses only the Python standard library.

```bash
# HLA-A and HLA-B (default)
python3 parse_specimmune.py results -o hla_calls.csv

# Add loci ("C" or "HLA-C"; non-HLA loci like TAP1 or MICA work too), or include every locus
python3 parse_specimmune.py results -g C DQB1 DRB1 -o hla_calls.csv
python3 parse_specimmune.py results --all-genes -o hla_calls.csv

# 2-field resolution, e.g. HLA-A*03:01
python3 parse_specimmune.py results --fields 2 -o hla_calls.csv
```

```
allele,HG00276
HLA-A_1,HLA-A*03:01:01:01
HLA-A_2,HLA-A*03:01:01:01
HLA-B_1,HLA-B*35:01:01:05
HLA-B_2,HLA-B*51:01:01:04
```

| option | |
|---|---|
| `inputs` | results folder(s) from `run_specimmune.sh`, or result files |
| `-o` | output CSV (default: print to screen) |
| `-g` | extra loci besides HLA-A and HLA-B |
| `--all-genes` | every locus in the result files |
| `--call` | `one_guess` (default), `top` or `all`. See below |
| `--fields N` | shorten alleles to N fields (default: keep all) |

SpecImmune lists several tied candidate alleles per haplotype in `Genotype`, best first. It also gives a single best guess in `One_guess`, which isn't always the first candidate (e.g. HLA-DPB1 in the sample).
- `one_guess` reports `One_guess`. If that is `NA`, as for loci typed without it (e.g. HLA-E, TAP1), it falls back to the first `Genotype` candidate.
- `top` always reports the first `Genotype` candidate.
- `all` reports every candidate, `;`-separated. Combined with `--fields 2`, duplicates are merged, which shows how ambiguous a call is at 2-field resolution.

Loci that weren't typed (`NA`) get empty cells.
