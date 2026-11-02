# wHLApper

Wrappers for running and comparing HLA (and other immune/PGx gene) typing tools on Oxford Nanopore long reads. Each tool lives in its own folder with its own Docker image, runner and README.

| folder | tool | what it types | how it runs |
|---|---|---|---|
| [`HLA-LA/`](HLA-LA/README.md) | [HLA*LA](https://github.com/DiltheyLab/HLA-LA) | HLA | Nextflow pipeline (Docker, Singularity or conda; local or Slurm) |
| [`SpecImmune/`](SpecImmune/README.md) | [SpecImmune](https://github.com/deepomicslab/SpecImmune) | HLA, KIR, CYP, IG/TR | bash + Docker (`run_specimmune.sh`) |
| [`FuFiHLA/`](FuFiHLA/) | FuFiHLA (bioconda `fufihla`) | HLA-A, -B, -C, -DQA1, -DQB1, -DRB1 | bash + Docker (`run_fufihla.sh`) |
| [`scripts/`](scripts/README.md) | helpers | | `fetch_ont_regions.sh` builds test data from ONT open data |

## Quick start

```bash
# HLA*LA (Nextflow)
nextflow run HLA-LA/main.nf -profile docker --input samplesheet.csv \
    --reference GRCh38.fa --hla_reference GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --graph_dir /ref/graph/PRG_MHC_GRCh38_withIMGT

# SpecImmune
SpecImmune/run_specimmune.sh -s samplesheet.csv -o results -d ./db -g HLA,KIR -r hg38_no_alt.fa

# FuFiHLA
FuFiHLA/run_fufihla.sh -s samplesheet.csv -o fufihla_results
```

All three take a samplesheet whose first two columns are `sample_id,bam`, so the same sheet can be reused across tools. See each folder's README for its image build, extra columns and options.

## Test data

[`scripts/fetch_ont_regions.sh`](scripts/README.md) downloads GIAB CRAMs from the ONT open data bucket and writes small per-sample BAMs for HLA, KIR, IG/TCR and CYP2D6. These can go straight into any of the tools above.

## Notes

- Inputs are expected to be on GRCh38/hg38 with `chr`-style contig names.
- Run `chmod +x` on the scripts after copying the repo from Windows.
