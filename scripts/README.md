# scripts

Helper scripts shared by all the typing tools in this repo. Tool-specific scripts live in each tool's folder; for example, the HLA*LA graph indexer is [`HLA-LA/prepare_hla_graph.sh`](../HLA-LA/README.md#setup).

## fetch_ont_regions.sh

Downloads `SAMPLE.haplotagged.cram` files from the [ONT open data](https://labs.epi2me.io/dataindex/) S3 bucket and cuts out small per-target BAMs for typing tests. The published `.crai` index can't be used, so each CRAM is re-indexed with `samtools index` after download. Each CRAM is deleted once all four of its BAMs are written. If extraction fails, the CRAM stays in `<work>/<sample_id>/` and a rerun extracts from it without downloading again.

Requires the `aws` CLI and Docker. The bucket is public, so no AWS credentials are needed. samtools runs from `ontresearch/wf-common:sha72f3517dd994984e0e2da0b97cb3f23f8540be4b` (change it with `-I`, or use `-L` for samtools on your PATH).

```bash
./fetch_ont_regions.sh -i samplesheet.csv -o giab_open_data.2025 -t 8 \
    -r /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Long runs
nohup bash -c './fetch_ont_regions.sh -i samplesheet.csv -o giab_open_data.2025 -t 8 -r ...' > download.log 2>&1 &
```

| option | |
|---|---|
| `-i` | samplesheet (required) |
| `-o` | output folder (required) |
| `-r` | reference FASTA to decode the CRAMs. Without it, samtools fetches the reference from EBI and caches it in the work folder (slow the first time) |
| `-t` | samtools threads (default 4) |
| `-w` | download folder (default `<out>/work`) |
| `-k` | keep the downloaded CRAMs |
| `-I` | Docker image that provides samtools (default `ontresearch/wf-common:sha72f3517dd994984e0e2da0b97cb3f23f8540be4b`) |
| `-L` | use samtools from your PATH instead of Docker |

### Samplesheet

A CSV with `sample_id,flowcell_id,s3_path`. The header row is optional. `s3_path` is the folder that contains `SAMPLE.haplotagged.cram`:

```
sample_id,flowcell_id,s3_path
HG001_1,PAW79146,s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/hac/HG001/PAW79146/output/
```

### Outputs

One sorted and indexed BAM per target, using GRCh38 regions:

| file | regions |
|---|---|
| `<sample_id>.HLA.bam` | `chr6`, any `chr6_*_alt` / `HLA-*` contigs, and unmapped reads |
| `<sample_id>.KIR.bam` | `chr19:54000000-56000000` |
| `<sample_id>.IGTCR.bam` | IGH, IGK, IGL, TRA/TRD, TRB, TRG |
| `<sample_id>.CYP2D6.bam` | `chr22:42000000-42300000` (CYP2D6, CYP2D7, CYP2D8P) |

To change the regions, edit `TARGETS` at the top of the script. Add `*` to a target to include unmapped reads.

Samples whose four BAMs already exist are skipped, so you can rerun after a failure. Failed samples are listed at the end, and the script exits with code 1.
