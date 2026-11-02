#!/usr/bin/env python3

################################
## Process HLA Typing Results ##
################################

## Brief QC - Check if Q1 > threshold Specified

import argparse
import sys
import logging

__version__ = 'v1.0.3'

LOCUS = ['A', 'B', 'C', 'DQA1', 'DQB1', 'DRB1']

class Locus:

    def __init__(self, input_dict):
        for key, value in input_dict.items():
            value = float(value) if key in ['Q1'] else value
            setattr(self, key, value)

    def check_threshold(self, threshold = 0.8):
        self.passed = self.Q1 > threshold
        fail_message = f'HLA-{self.Locus}-{self.Chromosome}-{str(self.Q1)}'
        self.message = 'Passed' if self.passed else fail_message

def process_file(args):
    ''' 
    Process Results from HLA-Typing 
    Params args: args object containintg threshold and file from HLA-LA
    '''
    results = []

    with open(args.results, 'r') as f:
        for c, line in enumerate(f):
            line = line.strip('\n').split('\t')
            if c == 0:
                header = line 
                continue

            tmp = Locus({k:v for k,v in zip(header, line)})
            if tmp.Locus not in LOCUS: continue
            tmp.check_threshold(args.threshold)
            results.append(tmp) 

    return results

def write_to_file(to_write, write_file, sep = '\t'):
    to_write = sep.join(to_write)
    write_file.write(f'{to_write}\n')

def write_wrapper(sample, results):
    attributes = ['Locus', 'Chromosome', 'Allele', 'Q1', 'AverageCoverage', 'MinimumCoverage', 'proportionkMersCovered', 'message']
    file_name = f'{sample}.hla_typing_results.txt'
    with open(file_name, 'w') as write_file:
        write_to_file(attributes, write_file)
        for r in results:
            to_write = [str(getattr(r, attr)) for attr in attributes]
            write_to_file(to_write, write_file)

def main(args):
    try:
        results = process_file(args)
        write_wrapper(args.sample, results)

    except Exception as Argument:
        logging.error(f"Error reading alleles from HLA-LA {__version__}: {Argument}")
        raise Exception("Error reading alleles from HLA-LA output")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process HLA Typing Results")

    parser.add_argument("-s", "--sample", type=str, required=True, 
                        dest="sample", help="sample id for the run")

    parser.add_argument("-r", "--results", type=str, required=True, 
                        dest="results", help="HLA-LA Results")
    
    parser.add_argument("-t", "--threshold", type=str, default=0.8, 
                        dest="threshold", help="Threshold Value")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    main(parser.parse_args())
