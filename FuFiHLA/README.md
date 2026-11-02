# FuFiHLA

Batch runner for FuFiHLA typing of long-read BAMs (ONT or PacBio HiFi) in Docker. It types six genes: HLA-A, -B, -C, -DQA1, -DQB1 and -DRB1. It is separate from the Nextflow pipeline.

```
FuFiHLA/
  docker/Dockerfile.fufihla   image: fufihla + samtools from bioconda
  run_fufihla.sh              runs FuFiHLA for every sample in a samplesheet
  parse_fufihla.py            summarises the calls into a CSV (alleles x samples)
  sample_output/              example <sample>.out
```

## 1. Build the image

```bash
cd FuFiHLA/docker
docker build -f Dockerfile.fufihla -t fufihla:v0.2.4 .
```

The Dockerfile installs the latest `fufihla` from bioconda, so the `v0.2.4` tag is only a label. To pin the version, change `fufihla` to `fufihla=0.2.4` in the `mamba install` line.

Check it works:

```bash
docker run --rm fufihla:v0.2.4 fufihla --help
```

## 2. Prepare the samplesheet

A CSV or TSV with `sample_id,bam`. The header row is optional:

```
sample_id,bam
HG00276,/data/pgx/cohort_1/HG00276/HG00276.hg38.bam
HG01190,/data/pgx/cohort_1/HG01190/HG01190.hg38.bam
```

- By default, BAMs must be **aligned to hg38** with `chr`-style names, and **indexed** (`.bai`/`.crai` next to the file). The runner only pulls reads from the six gene regions.
- For **unaligned** BAMs, or targeted/amplicon data, use `-A` to give FuFiHLA every read in the BAM.
- CRAMs need the reference passed with `-c`.

## 3. Run FuFiHLA

Run these from the `FuFiHLA/` folder:

```bash
chmod +x run_fufihla.sh parse_fufihla.py    # needed after copying from Windows

./run_fufihla.sh -s samplesheet.csv -o fufihla_results            # ONT (default)
./run_fufihla.sh -s samplesheet.csv -o fufihla_results -p hifi    # PacBio HiFi
./run_fufihla.sh -s samplesheet.csv -o fufihla_results -A         # unaligned / targeted BAMs

# Long runs
nohup ./run_fufihla.sh -s samplesheet.csv -o fufihla_results > fufihla.log 2>&1 &
```

| option | |
|---|---|
| `-s` | samplesheet (required) |
| `-o` | output folder (required) |
| `-p` | read type: `ont` or `hifi` (default `ont`) |
| `-i` | Docker image (default `fufihla:v0.2.4`) |
| `-t` | threads for read extraction (default 4) |
| `-R` | custom FuFiHLA reference folder from `fufihla-ref-prep` (default: the bundled IMGT set) |
| `-H` | homozygosity threshold (`--hom-threshold`) |
| `-c` | reference FASTA for CRAM inputs |
| `-A` | use all reads in the BAM (unaligned or targeted data) |
| `-D` | keep FuFiHLA's intermediate files (`--debug`) |

For each sample, the runner:

1. **Extracts reads.** It takes reads overlapping the six genes (hg38 regions from the FuFiHLA README) and drops secondary and supplementary alignments, so each read is used once. The reads are written as `<sample>.hla_reads.fa.gz`. With `-A`, it takes every read instead.
2. **Types** them with `fufihla --fa <reads> --ont` (or `--hifi`).

A sample with no extracted reads fails with a hint to check the alignment and index. Failed samples are recorded in `run_summary.tsv`, and the runner moves on to the next sample.

### Outputs

```
fufihla_results/
  run_summary.tsv                 sample, OK / FAILED / MISSING_BAM, reads used, log
  all_samples_hla_calls.tsv       every call in long format: sample, allele, consensus, exact_match, cs tag
  hla_6genes.hg38.bed             regions used for extraction
  <sample>/<sample>.out           FuFiHLA allele calls (see sample_output/)
  <sample>/<sample>/              FuFiHLA working files
  <sample>/<sample>.hla_reads.fa.gz
  <sample>/<sample>.log
```

`exact_match` in `all_samples_hla_calls.tsv` is `yes` when the consensus matches a known allele exactly (the cs tag is a plain match length, e.g. `cs:Z::2943`). `no` means it differs from the database, which may point to a novel allele.

## 4. Summarise the calls

`parse_fufihla.py` turns the `.out` files into one CSV with two rows per gene (`<gene>_1`, `<gene>_2`) and one column per sample. It uses only the Python standard library.

```bash
# HLA-A and HLA-B (default)
python3 parse_fufihla.py fufihla_results -o hla_calls.csv

# Add more genes ("C" or "HLA-C"), or include every gene found
python3 parse_fufihla.py fufihla_results -g C DQA1 DQB1 DRB1 -o hla_calls.csv
python3 parse_fufihla.py fufihla_results --all-genes -o hla_calls.csv

# 2-field resolution, e.g. HLA-A*03:01
python3 parse_fufihla.py fufihla_results --fields 2 -o hla_calls.csv
```

```
allele,HG00276
HLA-A_1,HLA-A*03:01:01:59
HLA-A_2,HLA-A*03:01:01:59
HLA-B_1,HLA-B*35:01:01:05
HLA-B_2,HLA-B*51:01:01:04
```

| option | |
|---|---|
| `inputs` | results folder(s) from `run_fufihla.sh`, or `.out` files |
| `-o` | output CSV (default: print to screen) |
| `-g` | extra genes besides HLA-A and HLA-B |
| `--all-genes` | every gene found in the `.out` files |
| `--fields N` | shorten alleles to N fields (default: keep all) |
| `--no-fill-homozygous` | leave `<gene>_2` empty when only one allele is reported |

FuFiHLA reports a single line for a gene it calls homozygous, as with HLA-A in `sample_output/HG00276.out`. By default that allele fills both `_1` and `_2`. Use `--no-fill-homozygous` to tell homozygous calls apart from missed alleles. A gene with no call gets empty cells.
