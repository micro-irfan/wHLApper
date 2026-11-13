#!/usr/bin/env python3
"""
Summarise SpecImmune CYP calls into one CSV: one row per sample and gene, in the
same layout as parse_chinook.py (diplotype *4/*5, allele_1, allele_2).

Input is the specimmune_wrapper.sh / run_specimmune.sh results folder
(<results_dir>/<sample>/typing/<sample>/<sample>.CYP.final.type.result.formatted.txt),
or the result files given directly. The file has the same columns as the HLA
result (Locus, Chromosome, Genotype, ..., Reads_num, ..., One_guess).

Which allele is reported (--call):
  one_guess  SpecImmune's single best guess (One_guess column); when that is NA,
             the first candidate in Genotype is used instead (default)
  top        the first candidate in Genotype
  all        every candidate in Genotype, ';'-separated

Alleles are written as star alleles without the gene prefix (CYP2D6*4 -> *4).
--no-suballele drops PharmVar suballele suffixes (*4.001 -> *4).

Examples:
    parse_specimmune_cyp.py cyp2d6_specimmune -o cyp2d6_specimmune_calls.csv
    parse_specimmune_cyp.py cyp2d6_specimmune -g CYP2C19 CYP2C9 -o cyp_calls.csv
    parse_specimmune_cyp.py results/HG00276/typing/HG00276/HG00276.CYP.final.type.result.formatted.txt
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

SUFFIX = ".CYP.final.type.result.formatted.txt"
DEFAULT_GENES = ["CYP2D6"]
COLUMNS = ["sample_id", "gene", "diplotype", "allele_1", "allele_2", "reads_num", "status"]


def star(allele, no_suballele):
    """'CYP2D6*4.001' -> '*4.001' (or '*4' with no_suballele); other strings unchanged."""
    if "*" in allele:
        allele = "*" + allele.split("*", 1)[1]
    if no_suballele:
        allele = re.sub(r"\.\d+$", "", allele)
    return allele


def pick(row, call, no_suballele):
    """The allele string to report for one Locus/Chromosome row ('' if not typed)."""
    candidates = [a for a in row["Genotype"].split(";") if a and a != "NA"]
    one_guess = row.get("One_guess", "NA")
    if call == "all":
        return ";".join(OrderedDict.fromkeys(star(a, no_suballele) for a in candidates))
    if call == "one_guess" and one_guess not in ("", "NA", "-"):
        return star(one_guess, no_suballele)
    return star(candidates[0], no_suballele) if candidates else ""


def parse_result(path, call, no_suballele):
    """Return {locus: {'alleles': {1: allele, 2: allele}, 'reads': Reads_num}} in file order."""
    calls = OrderedDict()
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        locus, hap = row["Locus"], row["Chromosome"]
        if hap not in ("1", "2"):
            continue
        entry = calls.setdefault(locus, {"alleles": {}, "reads": row.get("Reads_num", "")})
        entry["alleles"][int(hap)] = pick(row, call, no_suballele)
    return calls


def find_result_files(inputs):
    """Expand folders to their <sample>/typing/<sample>/<sample>.CYP...txt files."""
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


def main():
    ap = argparse.ArgumentParser(
        description="Summarise SpecImmune CYP results into a CSV (one row per sample and gene).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Which allele" + __doc__.split("Which allele", 1)[1],
    )
    ap.add_argument("inputs", nargs="+", help="results folder(s) from specimmune_wrapper.sh, or result files")
    ap.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    ap.add_argument("-g", "--genes", nargs="+", default=[], metavar="GENE",
                    help="additional genes besides CYP2D6, e.g. CYP2C19 CYP2C9")
    ap.add_argument("--all-genes", action="store_true", help="include every locus in the result files")
    ap.add_argument("--call", choices=["one_guess", "top", "all"], default="one_guess",
                    help="which allele to report (default: one_guess)")
    ap.add_argument("--no-suballele", action="store_true",
                    help="drop PharmVar suballele suffixes, e.g. *4.001 -> *4")
    args = ap.parse_args()

    files = find_result_files(args.inputs)
    if not files:
        sys.exit("ERROR: no SpecImmune CYP result files to parse")

    samples = OrderedDict()
    for f in files:
        sample = f.name[: -len(SUFFIX)] if f.name.endswith(SUFFIX) else f.stem
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_result(f, args.call, args.no_suballele)

    known = OrderedDict((locus, None) for calls in samples.values() for locus in calls)
    genes = DEFAULT_GENES + [g.strip().upper() for g in args.genes]
    if args.all_genes:
        genes += list(known)
    genes = list(OrderedDict.fromkeys(genes))
    for g in genes:
        if g not in known:
            print(f"WARNING: locus {g} not found in any result file", file=sys.stderr)

    n = 0
    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        w = csv.writer(out)
        w.writerow(COLUMNS)
        for sample, calls in samples.items():
            for g in genes:
                entry = calls.get(g, {"alleles": {}, "reads": ""})
                a1, a2 = entry["alleles"].get(1, ""), entry["alleles"].get(2, "")
                diplotype = f"{a1}/{a2}" if a1 or a2 else ""
                status = "genotyped" if a1 and a2 else "not genotyped"
                w.writerow([sample, g, diplotype, a1, a2, entry["reads"], status])
                n += 1
    finally:
        if out is not sys.stdout:
            out.close()

    if args.output != "-":
        print(f"Wrote {n} rows ({len(samples)} samples) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
