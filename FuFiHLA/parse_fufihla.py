#!/usr/bin/env python3
"""
Summarise FuFiHLA allele calls into one CSV: one row per allele (two per gene),
one column per sample.

Input is the run_fufihla.sh results folder (<results_dir>/<sample>/<sample>.out),
or .out files given directly. The sample name is the .out file name without .out.

FuFiHLA writes a single line for a gene when it calls it homozygous, so that
allele is reported for both <gene>_1 and <gene>_2 (use --no-fill-homozygous to
leave _2 empty instead).

Examples:
    parse_fufihla.py fufihla_results -o hla_calls.csv
    parse_fufihla.py fufihla_results -g C DRB1 --fields 2 -o hla_calls.csv
    parse_fufihla.py sample_output/HG00276.out
"""

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

DEFAULT_GENES = ["HLA-A", "HLA-B"]


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
    args = ap.parse_args()

    files = find_out_files(args.inputs)
    if not files:
        sys.exit("ERROR: no .out files to parse")

    samples = OrderedDict()
    for f in files:
        sample = f.name[: -len(".out")] if f.name.endswith(".out") else f.stem
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_out(f)

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


if __name__ == "__main__":
    main()
