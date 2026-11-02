#!/usr/bin/env python3

################################
## Process HLA Typing Results ##
################################

## Per-sample QC of HLA*LA R1_bestguess_G.txt
## A call passes when AverageCoverage > min_cov AND Q1 > min_q1

import argparse
import logging
import sys

__version__ = 'v1.1.0'

OUT_COLUMNS = ['Locus', 'Chromosome', 'Allele', 'Q1', 'AverageCoverage', 'MinimumCoverage',
               'proportionkMersCovered', 'QC', 'Reason']


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class Locus:

    def __init__(self, input_dict):
        for key, value in input_dict.items():
            setattr(self, key, value)
        self.q1 = to_float(input_dict.get('Q1'))
        self.coverage = to_float(input_dict.get('AverageCoverage'))

    def check_threshold(self, min_cov=20.0, min_q1=0.95):
        reasons = []
        if self.coverage is None or self.coverage <= min_cov:
            reasons.append(f'cov={self.AverageCoverage}<={min_cov:g}')
        if self.q1 is None or self.q1 <= min_q1:
            reasons.append(f'Q1={self.Q1}<={min_q1:g}')
        self.passed = not reasons
        self.QC = 'PASS' if self.passed else 'FAIL'
        self.Reason = ';'.join(reasons) if reasons else '.'


def process_file(args):
    '''
    Process Results from HLA-Typing
    Params args: args object containing thresholds and file from HLA-LA
    '''
    results = []

    with open(args.results, 'r') as f:
        for c, line in enumerate(f):
            line = line.rstrip('\n').split('\t')
            if c == 0:
                header = line
                continue
            if not any(line):
                continue

            tmp = Locus({k: v for k, v in zip(header, line)})
            tmp.check_threshold(args.min_cov, args.min_q1)
            results.append(tmp)

    return results


def write_results(sample, results, output=None):
    file_name = output or f'{sample}.hla_qc.tsv'
    with open(file_name, 'w') as write_file:
        write_file.write('\t'.join(OUT_COLUMNS) + '\n')
        for r in results:
            write_file.write('\t'.join(str(getattr(r, attr, 'NA')) for attr in OUT_COLUMNS) + '\n')


def main(args):
    try:
        results = process_file(args)
        write_results(args.sample, results, args.output)

    except Exception as Argument:
        logging.error(f"Error reading alleles from HLA-LA {__version__}: {Argument}")
        raise Exception("Error reading alleles from HLA-LA output")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="QC of HLA*LA typing results")

    parser.add_argument("-s", "--sample", type=str, required=True,
                        dest="sample", help="sample id for the run")

    parser.add_argument("-r", "--results", type=str, required=True,
                        dest="results", help="HLA-LA R1_bestguess_G.txt")

    parser.add_argument("--min-cov", type=float, default=20.0,
                        dest="min_cov", help="AverageCoverage must be greater than this (default: 20)")

    parser.add_argument("--min-q1", type=float, default=0.95,
                        dest="min_q1", help="Q1 must be greater than this (default: 0.95)")

    parser.add_argument("-o", "--output", type=str, default=None,
                        dest="output", help="Output file (default: <sample>.hla_qc.tsv)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    main(parser.parse_args())
