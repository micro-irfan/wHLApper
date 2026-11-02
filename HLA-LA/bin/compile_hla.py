#!/usr/bin/env python3

#####################################
## Compile HLA*LA Results (Cohort) ##
#####################################

## Inputs : per-sample <sample>.hla_qc.tsv (from hlala_qc.py), optional truth.tsv
## Outputs: hla_matrix.tsv, hla_qc_matrix.tsv
##          hla_truth_matrix.tsv, hla_truth_details.tsv (only with --truth)

import argparse
import csv
import os
import re
import sys

QC_SUFFIX = '.hla_qc.tsv'
TRUTH_FIELDS = 2
TRUTH_COLUMN = re.compile(r'^HLA-([A-Za-z0-9]+)_([12])$', re.IGNORECASE)


def read_qc(path):
    ''' Returns sample, ordered loci, {(locus, chromosome): row} '''
    name = os.path.basename(path)
    sample = name[:-len(QC_SUFFIX)] if name.endswith(QC_SUFFIX) else os.path.splitext(name)[0]
    loci, calls = [], {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['Locus'] not in loci:
                loci.append(row['Locus'])
            calls[(row['Locus'], row['Chromosome'])] = row
    return sample, loci, calls


def read_truth(path):
    ''' Returns ordered truth columns (HLA-A_1, ...), {sample: {column: allele}} '''
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        columns = []
        for col in reader.fieldnames[1:]:
            m = TRUTH_COLUMN.match(col.strip())
            if not m:
                sys.exit(f'Unrecognised truth column: {col} (expected e.g. HLA-A_1)')
            columns.append(f'HLA-{m.group(1).upper()}_{m.group(2)}')
        truth = {}
        id_col = reader.fieldnames[0]
        for row in reader:
            values = [row[c] for c in reader.fieldnames[1:]]
            truth[row[id_col].strip()] = {col: (v or '').strip() for col, v in zip(columns, values)}
    return columns, truth


def normalise(allele, gene=None, fields=TRUTH_FIELDS):
    ''' "HLA-A*11:01:01G" -> ("A", ["11", "01"]); None for empty/NA '''
    a = (allele or '').strip().upper()
    if a in ('', 'NA', 'N/A', '-', '.'):
        return None
    a = re.sub(r'^HLA-', '', a)
    if '*' not in a and gene:
        a = f'{gene}*{a}'
    g, _, rest = a.partition('*')
    parts = [re.sub(r'[A-Z]+$', '', p) for p in rest.split(':') if p]
    return g, parts[:fields]


def is_match(call_allele, truth_allele, gene):
    ''' Compare at 2 fields (or fewer if the truth has fewer). Ambiguous calls (a;b) match on any. '''
    t = normalise(truth_allele, gene)
    if t is None:
        return None
    for option in call_allele.split(';'):
        c = normalise(option, gene)
        if c and c[0] == t[0] and c[1][:len(t[1])] == t[1]:
            return True
    return False


def compare_gene(gene, calls, truths):
    '''
    calls : [row_chr1 or None, row_chr2 or None]
    truths: [truth_1 or '', truth_2 or '']
    Tries both pairings of calls to truths and keeps the one with most matches.
    Returns one (status, call_row) per truth slot.
    '''
    best, best_score = None, -1
    for order in ((0, 1), (1, 0)):
        res = []
        for t_idx, c_idx in enumerate(order):
            truth, call = truths[t_idx], calls[c_idx]
            if not truth or normalise(truth, gene) is None:
                res.append(('NA', call))
            elif call is None:
                res.append(('NA', None))
            else:
                res.append(('MATCH' if is_match(call['Allele'], truth, gene) else 'MISMATCH', call))
        score = sum(r[0] == 'MATCH' for r in res)
        if score > best_score:
            best, best_score = res, score
    return best


def fraction(n_match, n_total):
    return f'{n_match}/{n_total} ({100 * n_match / n_total:.1f}%)' if n_total else 'NA'


def write_matrix(path, row_names, samples, values, extra_col=None):
    with open(path, 'w') as f:
        header = ['Locus'] + samples + ([extra_col[0]] if extra_col else [])
        f.write('\t'.join(header) + '\n')
        for r in row_names:
            line = [r] + [values.get((r, s), 'NA') for s in samples]
            if extra_col:
                line.append(extra_col[1].get(r, 'NA'))
            f.write('\t'.join(line) + '\n')


def main(args):
    samples, loci, all_calls = [], [], {}
    for path in sorted(args.qc):
        sample, s_loci, calls = read_qc(path)
        samples.append(sample)
        all_calls[sample] = calls
        loci += [l for l in s_loci if l not in loci]

    rows = [(f'HLA-{l}_{c}', l, c) for l in loci for c in ('1', '2')]
    allele_values, qc_values = {}, {}
    for s in samples:
        for name, locus, chrom in rows:
            call = all_calls[s].get((locus, chrom))
            if call is None:
                continue
            failed = call['QC'] != 'PASS'
            allele_values[(name, s)] = call['Allele'] + (' (FAIL)' if failed else '')
            qc_values[(name, s)] = 'PASS' if not failed else f"FAIL:{call['Reason']}"

    row_names = [r[0] for r in rows]
    write_matrix(os.path.join(args.outdir, 'hla_matrix.tsv'), row_names, samples, allele_values)
    write_matrix(os.path.join(args.outdir, 'hla_qc_matrix.tsv'), row_names, samples, qc_values)

    if not args.truth:
        return

    truth_cols, truth = read_truth(args.truth)
    t_samples = samples + [s for s in truth if s not in samples]
    genes = []
    for col in truth_cols:
        g = TRUTH_COLUMN.match(col).group(1).upper()
        if g not in genes:
            genes.append(g)

    values, details = {}, []
    per_row = {c: [0, 0] for c in truth_cols}
    per_sample = {s: [0, 0] for s in t_samples}
    for s in t_samples:
        calls = all_calls.get(s, {})
        t = truth.get(s, {})
        for g in genes:
            slots = [f'HLA-{g}_1', f'HLA-{g}_2']
            result = compare_gene(g, [calls.get((g, '1')), calls.get((g, '2'))],
                                  [t.get(slots[0], ''), t.get(slots[1], '')])
            for slot, (status, call) in zip(slots, result):
                if slot not in truth_cols:
                    continue
                failed = call is not None and call['QC'] != 'PASS'
                values[(slot, s)] = status + (' (FAIL)' if failed and status != 'NA' else '')
                if status != 'NA':
                    per_row[slot][1] += 1
                    per_sample[s][1] += 1
                    if status == 'MATCH':
                        per_row[slot][0] += 1
                        per_sample[s][0] += 1
                details.append([s, slot, t.get(slot, '') or 'NA',
                                call['Allele'] if call else 'NA',
                                call['QC'] if call else 'NA', status])

    for s in t_samples:
        values[('Concordance', s)] = fraction(*per_sample[s])
    row_conc = {r: fraction(*per_row[r]) for r in truth_cols}
    row_conc['Concordance'] = fraction(sum(v[0] for v in per_row.values()),
                                       sum(v[1] for v in per_row.values()))
    write_matrix(os.path.join(args.outdir, 'hla_truth_matrix.tsv'), truth_cols + ['Concordance'],
                 t_samples, values, extra_col=('Concordance', row_conc))

    with open(os.path.join(args.outdir, 'hla_truth_details.tsv'), 'w') as f:
        f.write('\t'.join(['Sample', 'Locus', 'Expected', 'Called', 'QC', 'Status']) + '\n')
        for d in details:
            f.write('\t'.join(d) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compile HLA*LA results into cohort matrices")
    parser.add_argument("--qc", nargs='+', required=True, help="Per-sample *.hla_qc.tsv files")
    parser.add_argument("--truth", default=None, help="Optional truth.tsv (sample_id, HLA-A_1, ...)")
    parser.add_argument("--outdir", default='.', help="Output directory")
    main(parser.parse_args())
