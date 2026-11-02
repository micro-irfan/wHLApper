# HLA*LA

A Nextflow pipeline that runs [HLA*LA](https://github.com/DiltheyLab/HLA-LA) on aligned ONT BAMs/CRAMs. It QCs each call and compiles cohort matrices, with optional comparison against expected (truth) types.

## Layout

```
HLA-LA/
  main.nf                 pipeline: input parsing + workflow
  nextflow.config         params, resources, profiles
  module/                 processes, one file per tool
    samtools.nf             MERGE_BAMS, SAMTOOLS_INDEX, EXTRACT_HLA_REGION
    minimap2.nf             MINIMAP2_ALIGN (realign to the HLA decoy reference)
    hlala.nf                HLALA, HLALA_QC
    compile.nf              COMPILE
  bin/                    hlala_qc.py, compile_hla.py (on PATH for processes)
  assets/                 example samplesheet, NO_TRUTH placeholder
  envs/                   conda environments for the conda profile
  docker/Dockerfile       HLA*LA 1.0.3 image (samtools, bwa, picard, minimap2)
  prepare_hla_graph.sh    indexes the HLA*LA PRG graph
  run_HLA-LA.sh           launcher: runs the pipeline in the background with your settings
  reference/              the earlier PGx pipeline this was adapted from (not used)
```

## Setup

All commands below are run from this folder.

```bash
git clone https://github.com/micro-irfan/wHLApper.git
cd wHLApper

# 1. Build the image (slow the first time: Boost is compiled from source)
docker build -t hla-la:1.0.3 docker

# 2. Index the PRG graph once. Files are written into the graph folder
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz

# If the above does not work: 
wget -O PRG_MHC_GRCh38_withIMGT.tar.gz   "https://zenodo.org/records/19336310/files/PRG_MHC_GRCh38_withIMGT.tar.gz?download=1"

md5sum PRG_MHC_GRCh38_withIMGT.tar.gz
# md5sum for PRG_MHC_GRCh38_withIMGT.tar.gz is 525a8aa0c7f357bf29fe2c75ef1d477d

tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz

# This should take about 20mins on a GridION 
./prepare_hla_graph.sh -g /ref/graph/PRG_MHC_GRCh38_withIMGT

# 3. Optional: build a minimap2 index of the decoy reference, so it isn't rebuilt per sample
BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
wget -c $BASE/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c $BASE/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

minimap2 -x map-ont -d GRCh38_full_analysis_set_plus_decoy_hla.mmi GRCh38_full_analysis_set_plus_decoy_hla.fa
```

`prepare_hla_graph.sh` options: `-g` graph folder (required), `-i` image (default `hla-la:1.0.3`), `-f` re-index even if `serializedGRAPH` already exists. Indexing takes a long time and needs a lot of memory.

The default `--hla_container` in `nextflow.config` points at the team's ECR image. Pass `--hla_container hla-la:1.0.3` to use the one you built.

## Usage

```bash
# Single BAM/CRAM
nextflow run main.nf -profile docker \
    --input sample.bam \
    --reference /ref/GRCh38.fa \
    --hla_reference /ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --hla_reference_mmi /ref/GRCh38_full_analysis_set_plus_decoy_hla.mmi \
    --graph_dir /ref/graph/PRG_MHC_GRCh38_withIMGT

# Samplesheet
nextflow run main.nf -profile singularity \
    --input assets/samplesheet_example.csv --reference ... --hla_reference ... --graph_dir ...
```

### Running with `run_HLA-LA.sh`

For routine runs, use [`run_HLA-LA.sh`](run_HLA-LA.sh). It starts Nextflow in the background with `-resume`, so you can log out while it runs.

1. Make a run folder and copy in the script and your samplesheet:

   ```bash
   mkdir -p /data/runs/hlala_2026_11 && cd /data/runs/hlala_2026_11
   cp /path/to/wHLApper/HLA-LA/run_HLA-LA.sh .
   cp /path/to/samplesheet.csv .
   ```

2. Edit the settings at the top of `run_HLA-LA.sh`:

   | variable | set it to |
   |---|---|
   | `WORKFLOW_DIR` | `/path/to/wHLApper/HLA-LA`. Required once the script is copied out of the repo. You can also pass it as `WORKFLOW_DIR=... ./run_HLA-LA.sh` |
   | `RUN_NAME` | a name for this run. Results go to `results/<RUN_NAME>` |
   | `SAMPLESHEET` | your samplesheet (default `samplesheet.csv` in the run folder) |
   | `PROFILE` | `docker`, `singularity` or `conda`, optionally with `,slurm` |
   | `REFERENCE` | FASTA your BAMs/CRAMs were aligned to |
   | `HLA_REFERENCE`, `HLA_REFERENCE_MMI` | decoy FASTA and its minimap2 index (set the `.mmi` to `""` if you don't have one) |
   | `GRAPH_DIR` | graph folder prepared with `prepare_hla_graph.sh` |
   | `HLA_CONTAINER` | `hla-la:1.0.3` for the image you built, or the ECR image |

   Threads and QC thresholds (`--threads`, `--min_coverage`, `--min_q1`) are set as flags in the `nextflow run` command further down.

3. Launch it:

   ```bash
   bash run_HLA-LA.sh
   # Launching Nextflow in the background...
   # PID: 12345  |  log: nextflow_run.<RUN_NAME>.log
   ```

4. Follow progress, and stop the run if needed:

   ```bash
   tail -f nextflow_run.<RUN_NAME>.log
   kill <PID>        # stops Nextflow; running jobs are cleaned up
   ```

To resume after a failure or a `kill`, run `bash run_HLA-LA.sh` again from the same folder. `-resume` reuses finished tasks from the `work/` folder. When results are final, you can delete `work/` and `.nextflow*`.

Profiles: `docker`, `singularity` or `conda`. Any of them can be combined with `slurm`, e.g. `-profile singularity,slurm --queue long`. Whichever profile you use, the graph must already be prepared. The conda package doesn't include one.

| parameter | |
|---|---|
| `--input` | aligned `.bam`/`.cram`, or a samplesheet `.csv`/`.tsv` (required) |
| `--reference` | FASTA the inputs were aligned to, `.fai` alongside (required) |
| `--hla_reference` | decoy FASTA the graph was built against, e.g. `GRCh38_full_analysis_set_plus_decoy_hla.fa`, `.fai` alongside (required) |
| `--graph_dir` | prepared HLA*LA graph folder, e.g. `PRG_MHC_GRCh38_withIMGT` (required) |
| `--hla_reference_mmi` | prebuilt minimap2 index of `--hla_reference` |
| `--sample_id` | sample name in single-BAM mode (default: BAM file name) |
| `--outdir` | output folder (default `results`) |
| `--min_coverage`, `--min_q1` | QC thresholds (default `20`, `0.95`) |
| `--skip_extract` | send the whole BAM to realignment instead of only the chr6/HLA reads (much slower) |
| `--threads`, `--hla_memory`, `--align_memory` | resources for HLA*LA and minimap2 |
| `--hla_container`, `--python_container` | container images |

## Samplesheet (CSV or TSV)

| column | required | notes |
|---|---|---|
| `sample_id` | yes | |
| `bam` | yes | one BAM/CRAM, or a folder whose `*.bam`/`*.cram` files are merged (they must be coordinate-sorted) |
| `HLA-A_1`, `HLA-A_2`, `HLA-B_1`, `HLA-B_2` | no | expected types, e.g. `A*11:01`. Other genes also work, e.g. `HLA-C_1` |

## Steps

1. Merge (for folders) and index the alignments.
2. Extract the chr6, chr6 alt and `HLA-*` contigs plus unmapped reads. Use `--skip_extract` to turn this off.
3. Realign those reads with `minimap2 -x map-ont` to the HLA decoy reference (`--hla_reference`).
4. Run `HLA-LA.pl --longReads ont2d` against the decoy reference.
5. QC each call. It passes if `AverageCoverage > 20` and `Q1 > 0.95` (set with `--min_coverage` and `--min_q1`).
6. Compile the cohort matrices.

## Outputs

```
results/
  <sample>/hla/                    raw HLA*LA output (R1_bestguess_G.txt, ...)
  <sample>/<sample>.hla_qc.tsv     per-call QC
  summary/hla_matrix.tsv           loci x samples, 3-field G calls; failed calls tagged "(FAIL)"
  summary/hla_qc_matrix.tsv        PASS / FAIL:<reason>
  summary/hla_truth_matrix.tsv     MATCH / MISMATCH / NA + concordance (only when truth columns are given)
  summary/hla_truth_details.tsv    expected vs called, in long format
  pipeline_info/                   trace, timeline, report
```

Calls are compared with the truth at 2 fields. The order of the two alleles in a gene doesn't matter: `A*01/A*11` matches `A*11/A*01`. A call that failed QC is still compared but tagged, e.g. `MATCH (FAIL)`.
