#!/usr/bin/env python3
"""
Summarise SpecImmune HLA calls into one CSV: one row per allele (two per locus),
one column per sample.

Input is the run_specimmune.sh results folder
(<results_dir>/<sample>/typing/<sample>/<sample>.HLA.final.type.result.formatted.txt),
or the result files given directly. The sample name is the part of the file name
before ".HLA.final.type.result.formatted.txt".

Which allele is reported (--call):
  one_guess  SpecImmune's single best guess (One_guess column); when that is NA,
             the first candidate in Genotype is used instead (default)
  top        the first candidate in Genotype
  all        every candidate in Genotype, ';'-separated (deduplicated after --fields)

Examples:
    parse_specimmune.py results -o hla_calls.csv
    parse_specimmune.py results -g C DQB1 DRB1 --fields 2 -o hla_calls.csv
    parse_specimmune.py results --all-genes --call all -o hla_candidates.csv
    parse_specimmune.py sample_output/HG00276.HLA.final.type.result.formatted.txt
"""

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

SUFFIX = ".HLA.final.type.result.formatted.txt"
DEFAULT_GENES = ["HLA-A", "HLA-B"]


def truncate(allele, fields):
    """HLA-A*03:01:01:01 -> HLA-A*03:01 with fields=2. fields=0 keeps it as is."""
    if not fields or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    return f"{gene}*{':'.join(rest.split(':')[:fields])}"


def pick(row, call, fields):
    """The allele string to report for one Locus/Chromosome row ('' if not typed)."""
    candidates = [a for a in row["Genotype"].split(";") if a and a != "NA"]
    one_guess = row.get("One_guess", "NA")
    if call == "all":
        return ";".join(OrderedDict.fromkeys(truncate(a, fields) for a in candidates))
    if call == "one_guess" and one_guess not in ("", "NA", "-"):
        return truncate(one_guess, fields)
    return truncate(candidates[0], fields) if candidates else ""


def parse_result(path, call, fields):
    """Return {locus: {1: allele, 2: allele}} in file order."""
    calls = OrderedDict()
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        locus, hap = row["Locus"], row["Chromosome"]
        if hap not in ("1", "2"):
            continue
        calls.setdefault(locus, {})[int(hap)] = pick(row, call, fields)
    return calls


def find_result_files(inputs):
    """Expand folders to their <sample>/typing/<sample>/<sample>.HLA...txt files."""
    files = []
    for p in map(Path, inputs):
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            found = (sorted(p.glob(f"*/typing/*/*{SUFFIX}"))
                     or sorted(p.glob(f"*{SUFFIX}"))
                     or sorted(p.rglob(f"*{SUFFIX}")))
            if not found:
                print(f"WARNING: no *{SUFFIX} files found in {p}", file=sys.stderr)
            files.extend(found)
        else:
            sys.exit(f"ERROR: not found: {p}")
    return files


def resolve_gene(gene, known):
    """'c' -> 'HLA-C', 'tap1' -> 'TAP1', matching the loci found in the results."""
    g = gene.strip().upper()
    for cand in (g, f"HLA-{g}"):
        if cand in known:
            return cand
    return g if g.startswith("HLA-") else f"HLA-{g}"


def main():
    ap = argparse.ArgumentParser(
        description="Summarise SpecImmune HLA results into a CSV (rows: alleles, columns: samples).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Which allele" + __doc__.split("Which allele", 1)[1],
    )
    ap.add_argument("inputs", nargs="+", help="results folder(s) from run_specimmune.sh, or result files")
    ap.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    ap.add_argument("-g", "--genes", nargs="+", default=[], metavar="GENE",
                    help="additional loci besides HLA-A and HLA-B, e.g. C DQB1 DRB1 MICA")
    ap.add_argument("--all-genes", action="store_true",
                    help="include every locus in the result files")
    ap.add_argument("--call", choices=["one_guess", "top", "all"], default="one_guess",
                    help="which allele to report (default: one_guess)")
    ap.add_argument("--fields", type=int, default=0,
                    help="truncate alleles to this many fields, e.g. 2 -> A*03:01 (default: keep all)")
    args = ap.parse_args()

    files = find_result_files(args.inputs)
    if not files:
        sys.exit("ERROR: no SpecImmune result files to parse")

    samples = OrderedDict()
    for f in files:
        sample = f.name[: -len(SUFFIX)] if f.name.endswith(SUFFIX) else f.stem
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_result(f, args.call, args.fields)

    known = OrderedDict((locus, None) for calls in samples.values() for locus in calls)
    genes = DEFAULT_GENES + [resolve_gene(g, known) for g in args.genes]
    if args.all_genes:
        genes += sorted(known)
    genes = list(OrderedDict.fromkeys(genes))
    for g in genes:
        if g not in known:
            print(f"WARNING: locus {g} not found in any result file", file=sys.stderr)

    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        w = csv.writer(out)
        w.writerow(["allele"] + list(samples))
        for g in genes:
            for hap in (1, 2):
                w.writerow([f"{g}_{hap}"] + [samples[s].get(g, {}).get(hap, "") for s in samples])
    finally:
        if out is not sys.stdout:
            out.close()

    if args.output != "-":
        print(f"Wrote {2 * len(genes)} rows x {len(samples)} samples to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
