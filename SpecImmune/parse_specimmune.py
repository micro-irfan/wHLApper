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

Read counts, coverage and mean depth of SpecImmune's step-2 BAMs
(<results_dir>/<sample>/typing/<sample>/Genes_step2/*.bam) are written to a second
CSV via `samtools coverage` (default: <output>.coverage.csv); skipped if there are none.

Examples:
    parse_specimmune.py results -o hla_calls.csv
    parse_specimmune.py results -g C DQB1 DRB1 --fields 2 -o hla_calls.csv
    parse_specimmune.py results --all-genes --call all -o hla_candidates.csv
    parse_specimmune.py sample_output/HG00276.HLA.final.type.result.formatted.txt
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

SUFFIX = ".HLA.final.type.result.formatted.txt"
DEFAULT_GENES = ["HLA-A", "HLA-B"]
COVERAGE_COLUMNS = ["sample_id", "bam", "contig", "gene", "haplotype", "allele",
                    "numreads", "coverage", "meandepth"]


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


def gene_of(contig, bam_name):
    """'HLA-B_ref1' / 'HLA-B.0' -> 'HLA-B'; non-HLA BAMs ('HFE.0') fall back to the BAM name ('HFE')."""
    for text in (contig, bam_name):
        m = (re.search(r"HLA-([A-Z]+[0-9]*)", text, re.I)
             or re.search(r"(?<![A-Za-z0-9])([A-Z]+[0-9]*)\*[0-9]", text))
        if m:
            return f"HLA-{m.group(1).upper()}"
    return re.sub(r"\.\d+$", "", bam_name)


def haplotypes_of(bam_name):
    """'HLA-B.0' -> [1], 'HLA-B.1' -> [2]; a single BAM for the gene ('HLA-A') covers both -> [1, 2]."""
    m = re.search(r"\.(\d+)$", bam_name)
    return [int(m.group(1)) + 1] if m else [1, 2]


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


def write_coverage(path, bams_by_sample, genes, samples):
    """
    One row per haplotype of each BAM contig of the requested genes (genes=None
    keeps every contig). A gene SpecImmune aligned to a single reference (one BAM,
    e.g. homozygous) is reported for both haplotypes with the same numbers.
    The allele is the call from the result file for that gene and haplotype.
    """
    rows, unmatched = [], set()
    for sample, bams in bams_by_sample.items():
        for bam in bams:
            name = bam.name[: -len(".bam")]
            for r in samtools_coverage(bam):
                gene = gene_of(r["rname"], name)
                if genes is not None and gene not in genes:
                    unmatched.add(r["rname"])
                    continue
                haps = haplotypes_of(name)
                for hap in haps:
                    allele = samples[sample].get(gene, {}).get(hap, "")
                    # Split single-reference genes are labelled like the rest: HLA-A.0 / HLA-A_ref1
                    bam_label, contig = ((f"{name}.{hap - 1}", f"{r['rname']}_ref{hap}") if len(haps) > 1
                                         else (name, r["rname"]))
                    rows.append([sample, bam_label, contig, gene, hap, allele,
                                 r["numreads"], r["coverage"], r["meandepth"]])
    if not rows and unmatched:
        print("WARNING: no coverage rows matched the loci; contigs seen: "
              + ", ".join(sorted(unmatched)[:5]), file=sys.stderr)

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(COVERAGE_COLUMNS)
        w.writerows(rows)
    print(f"Wrote {len(rows)} coverage rows to {path}", file=sys.stderr)


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
    ap.add_argument("--coverage", metavar="CSV",
                    help="coverage CSV from the Genes_step2 BAMs (default: <output>.coverage.csv)")
    ap.add_argument("--no-coverage", action="store_true", help="skip the samtools coverage step")
    args = ap.parse_args()

    files = find_result_files(args.inputs)
    if not files:
        sys.exit("ERROR: no SpecImmune result files to parse")

    samples, sources = OrderedDict(), OrderedDict()
    for f in files:
        sample = f.name[: -len(SUFFIX)] if f.name.endswith(SUFFIX) else f.stem
        if sample in samples:
            print(f"WARNING: duplicate sample {sample} ({f}), keeping the first", file=sys.stderr)
            continue
        samples[sample] = parse_result(f, args.call, args.fields)
        sources[sample] = f

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

    # Coverage of SpecImmune's step-2 BAMs (Genes_step2/ next to each result file)
    if args.no_coverage:
        return
    cov_out = args.coverage or (None if args.output == "-" else f"{Path(args.output).with_suffix('')}.coverage.csv")
    bams = OrderedDict((s, sorted((f.parent / "Genes_step2").glob("*.bam"))) for s, f in sources.items())
    if not any(bams.values()):
        print("NOTE: no Genes_step2/*.bam next to the SpecImmune result files; coverage skipped", file=sys.stderr)
        return
    if not cov_out:
        print("NOTE: coverage skipped - give --coverage FILE when the calls go to stdout", file=sys.stderr)
        return
    if not shutil.which("samtools"):
        print("WARNING: samtools not on PATH - coverage skipped", file=sys.stderr)
        return
    for s in [s for s, b in bams.items() if not b]:
        print(f"WARNING: no BAMs for {s}", file=sys.stderr)
    write_coverage(cov_out, bams, None if args.all_genes else set(genes), samples)


if __name__ == "__main__":
    main()
