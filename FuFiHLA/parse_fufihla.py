#!/usr/bin/env python3
"""
Summarise FuFiHLA allele calls into one CSV: one row per allele (two per gene),
one column per sample.

Input is the run_fufihla.sh results folder (<results_dir>/<sample>/<sample>.out),
or .out files given directly. The sample name is the .out file name without .out.

FuFiHLA writes a single line for a gene when it calls it homozygous, so that
allele is reported for both <gene>_1 and <gene>_2 (use --no-fill-homozygous to
leave _2 empty instead).

If run_fufihla.sh was run with -D (keep intermediates), read counts, coverage and
mean depth of each assembly BAM (<results_dir>/<sample>/<sample>/bam/*_asm*.bam)
are written to a second CSV via `samtools coverage` (default: <output>.coverage.csv).
Without -D those BAMs don't exist and this step is skipped.

Examples:
    parse_fufihla.py fufihla_results -o hla_calls.csv
    parse_fufihla.py fufihla_results -g C DRB1 --fields 2 -o hla_calls.csv
    parse_fufihla.py sample_output/HG00276.out
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

DEFAULT_GENES = ["HLA-A", "HLA-B"]
COVERAGE_COLUMNS = ["sample_id", "bam", "contig", "gene", "haplotype", "allele",
                    "numreads", "coverage", "meandepth"]


def gene_name(gene):
    """'c' / 'HLA-C' / 'hla-c' -> 'HLA-C'."""
    gene = gene.strip().upper()
    return gene if gene.startswith("HLA-") else f"HLA-{gene}"


def truncate(allele, fields):
    """HLA-A*03:01:01:59 -> HLA-A*03:01 with fields=2. fields=0 keeps it as is."""
    if not fields or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    return f"{gene}*{':'.join(rest.split(':')[:fields])}"


def parse_out(path):
    """Return {gene: [allele, ...]} in file order.

    Call lines start with the allele (e.g. 'HLA-B*35:01:01:05'), followed by a
    tab or a space; everything else (log lines like '[FuFiHLA] ...', 'Step 1: ...')
    is skipped.
    """
    calls = OrderedDict()
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or not fields[0].startswith("HLA-") or "*" not in fields[0]:
                continue
            allele = fields[0]
            calls.setdefault(allele.split("*", 1)[0], []).append(allele)
    return calls


def debug_bams(out_file):
    """Assembly BAMs FuFiHLA keeps with --debug: <out>/bam/*_asm*.bam next to <out>.out."""
    out_dir = out_file.with_name(out_file.name[: -len(".out")] if out_file.name.endswith(".out") else out_file.stem)
    return sorted((out_dir / "bam").glob("*_asm*.bam"))


def gene_of(contig, bam_name):
    """'HLA-A*03:01:01' / 'A*03:01' / 'HLA-A_asm1' -> 'HLA-A' ('' if no gene is recognisable)."""
    for text in (contig, bam_name):
        m = (re.search(r"HLA-([A-Z]+[0-9]*)", text, re.I)
             or re.search(r"(?<![A-Za-z0-9])([A-Z]+[0-9]*)\*[0-9]", text))
        if m:
            return f"HLA-{m.group(1).upper()}"
    return ""


def haplotypes_of(bam_name):
    """'HLA-A_asm1' -> [1], 'HLA-A_asm2' -> [2]; no assembly number -> both haplotypes [1, 2]."""
    m = re.search(r"_asm(\d+)$", bam_name)
    return [int(m.group(1))] if m else [1, 2]


def contig_allele(contig, fields):
    """Assembly contig 'HLA-A*03_01_01_01' -> allele 'HLA-A*03:01:01:01' (then --fields)."""
    if "*" not in contig:
        return ""
    gene, rest = contig.split("*", 1)
    return truncate(f"{gene}*{rest.replace('_', ':')}", fields)


def samtools_coverage(bam):
    """Rows of `samtools coverage <bam>` as dicts keyed by its header (rname, numreads, ...)."""
    res = subprocess.run(["samtools", "coverage", str(bam)], capture_output=True, text=True)
    if res.returncode != 0:
        print(f"WARNING: samtools coverage failed for {bam}: {res.stderr.strip()}", file=sys.stderr)
        return []
    lines = res.stdout.splitlines()
    if not lines:
        return []
    header = lines[0].lstrip("#").split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:] if l.strip()]


def write_coverage(path, bams_by_sample, genes, fields):
    """
    One row per haplotype of each BAM contig of the requested genes (genes=None
    keeps every contig). The allele is the database allele the assembly contig
    was aligned to.
    """
    rows, unmatched = [], set()
    for sample, bams in bams_by_sample.items():
        for bam in bams:
            name = bam.name[: -len(".bam")]
            for r in samtools_coverage(bam):
                gene = gene_of(r["rname"], name)
                if genes is not None and gene not in genes:
                    if not gene:
                        unmatched.add(r["rname"])
                    continue
                for hap in haplotypes_of(name):
                    rows.append([sample, name, r["rname"], gene, hap, contig_allele(r["rname"], fields),
                                 r["numreads"], r["coverage"], r["meandepth"]])
    if not rows and unmatched:
        print("WARNING: no coverage rows matched the genes; contigs seen: "
              + ", ".join(sorted(unmatched)[:5]), file=sys.stderr)

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COVERAGE_COLUMNS)
        w.writerows(rows)
    print(f"Wrote {len(rows)} coverage rows to {path}", file=sys.stderr)


def find_out_files(inputs):
    """Expand folders to their <sample>/<sample>.out (or *.out) files."""
    files = []
    for p in map(Path, inputs):
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            found = sorted(p.glob("*/*.out")) or sorted(p.glob("*.out"))
            if not found:
                print(f"WARNING: no .out files found in {p}", file=sys.stderr)
            files.extend(found)
        else:
            sys.exit(f"ERROR: not found: {p}")
    return files


def main():
    ap = argparse.ArgumentParser(
        description="Summarise FuFiHLA .out files into a CSV (rows: alleles, columns: samples).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:" + __doc__.split("Examples:", 1)[1],
    )
    ap.add_argument("inputs", nargs="+", help="results folder(s) from run_fufihla.sh, or .out files")
    ap.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    ap.add_argument("-g", "--genes", nargs="+", default=[], metavar="GENE",
                    help="additional genes besides HLA-A and HLA-B, e.g. C DQA1 DQB1 DRB1")
    ap.add_argument("--all-genes", action="store_true",
                    help="include every gene found in the .out files")
    ap.add_argument("--fields", type=int, default=0,
                    help="truncate alleles to this many fields, e.g. 2 -> A*03:01 (default: keep all)")
    ap.add_argument("--no-fill-homozygous", action="store_true",
                    help="leave <gene>_2 empty when FuFiHLA reports a single allele")
    ap.add_argument("--coverage", metavar="CSV",
                    help="coverage CSV from the -D assembly BAMs (default: <output>.coverage.csv)")
    ap.add_argument("--no-coverage", action="store_true", help="skip the samtools coverage step")
    args = ap.parse_args()

    files = find_out_files(args.inputs)
    if not files:
        sys.exit("ERROR: no .out files to parse")

    samples, sources = OrderedDict(), OrderedDict()
    for f in files:
        sample = f.name[: -len(".out")] if f.name.endswith(".out") else f.stem
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_out(f)
        sources[sample] = f

    genes = DEFAULT_GENES + [g for g in map(gene_name, args.genes) if g not in DEFAULT_GENES]
    if args.all_genes:
        found = sorted({g for calls in samples.values() for g in calls})
        genes += [g for g in found if g not in genes]
    genes = list(OrderedDict.fromkeys(genes))

    # rows: HLA-A_1, HLA-A_2, ... ; columns: samples
    table = OrderedDict((f"{g}_{i}", {}) for g in genes for i in (1, 2))
    for sample, calls in samples.items():
        for g in genes:
            alleles = calls.get(g, [])
            if len(alleles) > 2:
                print(f"WARNING: {sample} has {len(alleles)} {g} calls, keeping the first two",
                      file=sys.stderr)
            if len(alleles) == 1 and not args.no_fill_homozygous:
                alleles = alleles * 2
            for i in (1, 2):
                table[f"{g}_{i}"][sample] = truncate(alleles[i - 1], args.fields) if len(alleles) >= i else ""

    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        w = csv.writer(out)
        w.writerow(["allele"] + list(samples))
        for row, values in table.items():
            w.writerow([row] + [values[s] for s in samples])
    finally:
        if out is not sys.stdout:
            out.close()

    if args.output != "-":
        print(f"Wrote {len(table)} rows x {len(samples)} samples to {args.output}", file=sys.stderr)

    # Coverage of the assembly BAMs, only present when run_fufihla.sh used -D
    if args.no_coverage:
        return
    cov_out = args.coverage or (None if args.output == "-" else f"{Path(args.output).with_suffix('')}.coverage.csv")
    bams = OrderedDict((s, debug_bams(f)) for s, f in sources.items())
    if not any(bams.values()):
        print("NOTE: no FuFiHLA debug BAMs (<sample>/<sample>/bam/*_asm*.bam) - "
              "rerun run_fufihla.sh with -D for coverage; skipped", file=sys.stderr)
        return
    if not cov_out:
        print("NOTE: coverage skipped - give --coverage FILE when the calls go to stdout", file=sys.stderr)
        return
    if not shutil.which("samtools"):
        print("WARNING: samtools not on PATH - coverage skipped", file=sys.stderr)
        return
    for s in [s for s, b in bams.items() if not b]:
        print(f"WARNING: no debug BAMs for {s}", file=sys.stderr)
    write_coverage(cov_out, bams, None if args.all_genes else set(genes), args.fields)


if __name__ == "__main__":
    main()
