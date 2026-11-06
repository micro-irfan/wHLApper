#!/usr/bin/env python3
"""
Gather HLA*LA calls, coverage and QC into one CSV: one row per allele (two per
gene), in the same long layout as the FuFiHLA / SpecImmune coverage CSVs.

Input is the wHLApper results folder (<outdir>/<sample>/hla/R1_bestguess_G.txt),
or R1_bestguess_G.txt files given directly. The sample name is the folder that
holds hla/.

HLA*LA reports coverage per locus, so both haplotypes of a gene share the same
numbers. Columns taken from R1_bestguess_G.txt:
  Q1, Q2                  HLA*LA call quality
  AverageCoverage         mean read depth over the locus (cf. samtools meandepth)
  CoverageFirstDecile     10th-percentile depth
  MinimumCoverage         lowest depth
  proportionkMersCovered  fraction of the allele's k-mers seen in reads
QC is PASS when AverageCoverage > --min-cov and Q1 > --min-q1, as in the pipeline.

Examples:
    parse_hlala.py results -o hlala_coverage.csv
    parse_hlala.py results -g C DRB1 --fields 2 -o hlala_coverage.csv
    parse_hlala.py results/HG00276/hla/R1_bestguess_G.txt
"""

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

RESULT_FILE = "R1_bestguess_G.txt"
DEFAULT_GENES = ["HLA-A", "HLA-B"]
METRICS = ["Q1", "Q2", "AverageCoverage", "CoverageFirstDecile", "MinimumCoverage", "proportionkMersCovered"]
COLUMNS = ["sample_id", "gene", "haplotype", "allele"] + METRICS + ["QC", "Reason"]


def gene_name(gene):
    """'a' / 'A' / 'HLA-A' -> 'HLA-A'."""
    gene = gene.strip().upper()
    return gene if gene.startswith("HLA-") else f"HLA-{gene}"


def truncate(allele, fields):
    """HLA-A*03:01:01G -> HLA-A*03:01 with fields=2. fields=0 keeps it as is."""
    if not fields or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    return f"{gene}*{':'.join(rest.split(':')[:fields])}"


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def qc(row, min_cov, min_q1):
    """(PASS/FAIL, reason) with the same thresholds as the pipeline's hlala_qc.py."""
    reasons = []
    cov, q1 = to_float(row.get("AverageCoverage")), to_float(row.get("Q1"))
    if cov is None or cov <= min_cov:
        reasons.append(f"cov={row.get('AverageCoverage')}<={min_cov:g}")
    if q1 is None or q1 <= min_q1:
        reasons.append(f"Q1={row.get('Q1')}<={min_q1:g}")
    return ("FAIL", ";".join(reasons)) if reasons else ("PASS", ".")


def parse_bestguess(path):
    """Return {gene: {haplotype: row}} from R1_bestguess_G.txt (Locus 'A' -> 'HLA-A')."""
    calls = OrderedDict()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get("Locus") or row.get("Chromosome") not in ("1", "2"):
                continue
            calls.setdefault(gene_name(row["Locus"]), {})[int(row["Chromosome"])] = row
    return calls


def find_result_files(inputs):
    """Expand folders to their <sample>/hla/R1_bestguess_G.txt files."""
    files = []
    for p in map(Path, inputs):
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            found = (sorted(p.glob(f"*/hla/{RESULT_FILE}"))
                     or sorted(p.glob(f"hla/{RESULT_FILE}"))
                     or sorted(p.rglob(RESULT_FILE)))
            if not found:
                print(f"WARNING: no {RESULT_FILE} found in {p}", file=sys.stderr)
            files.extend(found)
        else:
            sys.exit(f"ERROR: not found: {p}")
    return files


def sample_of(path):
    """<outdir>/<sample>/hla/R1_bestguess_G.txt -> <sample>."""
    parent = path.resolve().parent
    return parent.parent.name if parent.name == "hla" else parent.name


def main():
    ap = argparse.ArgumentParser(
        description="Gather HLA*LA calls, coverage and QC into one CSV (one row per allele).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="HLA*LA reports" + __doc__.split("HLA*LA reports", 1)[1],
    )
    ap.add_argument("inputs", nargs="+", help="wHLApper results folder(s), or R1_bestguess_G.txt files")
    ap.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    ap.add_argument("-g", "--genes", nargs="+", default=[], metavar="GENE",
                    help="additional genes besides HLA-A and HLA-B, e.g. C DQA1 DQB1 DRB1")
    ap.add_argument("--all-genes", action="store_true", help="include every gene in the result files")
    ap.add_argument("--fields", type=int, default=0,
                    help="truncate alleles to this many fields, e.g. 2 -> A*03:01 (default: keep all)")
    ap.add_argument("--min-cov", type=float, default=20.0,
                    help="QC: AverageCoverage must be greater than this (default: 20)")
    ap.add_argument("--min-q1", type=float, default=0.95,
                    help="QC: Q1 must be greater than this (default: 0.95)")
    args = ap.parse_args()

    files = find_result_files(args.inputs)
    if not files:
        sys.exit(f"ERROR: no {RESULT_FILE} files to parse")

    samples = OrderedDict()
    for f in files:
        sample = sample_of(f)
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_bestguess(f)

    genes = DEFAULT_GENES + [g for g in map(gene_name, args.genes) if g not in DEFAULT_GENES]
    if args.all_genes:
        found = OrderedDict((g, None) for calls in samples.values() for g in calls)
        genes += [g for g in found if g not in genes]
    genes = list(OrderedDict.fromkeys(genes))

    n = 0
    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        w = csv.writer(out)
        w.writerow(COLUMNS)
        for sample, calls in samples.items():
            for g in genes:
                for hap in (1, 2):
                    row = calls.get(g, {}).get(hap)
                    if row is None:
                        print(f"WARNING: {sample} has no {g} haplotype {hap}", file=sys.stderr)
                        continue
                    allele = truncate(gene_name(row["Allele"]), args.fields) if row.get("Allele") else ""
                    w.writerow([sample, g, hap, allele] + [row.get(m, "") for m in METRICS]
                               + list(qc(row, args.min_cov, args.min_q1)))
                    n += 1
    finally:
        if out is not sys.stdout:
            out.close()

    if args.output != "-":
        print(f"Wrote {n} rows ({len(samples)} samples) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
