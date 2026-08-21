#! /usr/bin/env python3

"""
Script to read RNASeqQC output and write a JSON results file

Input files are written by the RNASeqQC WDL workflow:
https://github.com/oicr-gsi/RNASeqQC

Based on collate_metrics.py in the rnaseq-qc Niassa workflow:
https://github.com/oicr-gsi/rnaseq-qc/blob/6143b0cc2ab3f54731581f2396b1ae46efeaaac5/workflow-rnaseq-qc/workflow/bin/collate_metrics.py

"""

import argparse
import json
import os
import re
import sys


class collator:
    BAMQC_KEY = 'bamqc'
    PICARD_KEY = 'picard'
    CONTAM_KEY = 'rrna_contamination'
    UNIQUE_READS_KEY = 'unique_reads'
    NO_STRAND_SPECIFICITY_KEY = 'no_strand_specificity'

    ''' set these Picard metrics to null for non-strand-specific libraries'''
    STRAND_SPECIFIC_METRICS = [
        'CORRECT_STRAND_READS',
        'INCORRECT_STRAND_READS',
        'NUM_R1_TRANSCRIPT_STRAND_READS',
        'NUM_R2_TRANSCRIPT_STRAND_READS',
        'NUM_UNEXPLAINED_READS',
        'PCT_CORRECT_STRAND_READS',
        'PCT_R1_TRANSCRIPT_STRAND_READS',
        'PCT_R2_TRANSCRIPT_STRAND_READS',
    ]

    def __init__(self, config):
        for key in config.keys():
            val = config.get(key, None)
            if val is None:
                raise ValueError("Config value for %s is required" % key)
            elif key != self.NO_STRAND_SPECIFICITY_KEY and not os.path.exists(val):
                raise FileNotFoundError("Input file '%s' for %s does not exist" % (val, key))
        ''' need to set self.strand_specific before parsing Picard input '''
        self.strand_specific = not config[self.NO_STRAND_SPECIFICITY_KEY]
        self.bamqc_metrics = self.parse_bamqc(config[self.BAMQC_KEY])
        self.picard_metrics = self.parse_picard(config[self.PICARD_KEY])
        self.contam_metrics = self.parse_contam(config[self.CONTAM_KEY])
        self.unique_reads = self.parse_unique_reads(config[self.UNIQUE_READS_KEY])

    @staticmethod
    def parse_bamqc(input_path):
        return json.loads(open(input_path).read())

    @staticmethod
    def parse_contam(input_path):
        lines = open(input_path).readlines()
        output = {}
        simpleKeys = ["mapped", "properly paired", "singletons"]
        for line in lines:
            words = re.split(" ", line.strip())
            count = int(words.pop(0))
            if words[0] == '+' and words[1] == '0': words = words[2:]
            key = " ".join(words)
            ''' strip out additional chatter in the output, to have consistent keys'''
            for simpleKey in simpleKeys:
                if re.match(simpleKey, key):
                    key = simpleKey
                    break
            output[key] = count
        return output

    def parse_picard(self, input_path):
        lines = open(input_path).readlines()
        header = []
        in_header = False
        output = {}
        for i in range(len(lines)):
            line = lines[i].strip()
            if re.match('## htsjdk.samtools.metrics.StringHeader$', line):
                in_header = not in_header
            elif in_header:
                header.append(line)
            elif re.match(r'## METRICS CLASS\s+picard.analysis.RnaSeqMetrics$', line):
                keys = re.split(r'\t', lines[i + 1].strip())
                values = re.split(r'\t', lines[i + 2].strip())
                if len(keys) < len(values):
                    ''' 
                        may have len(keys) > len(values)
                        final fields ['SAMPLE', 'LIBRARY', 'INPUT_GROUP'] are optional
                    '''
                    raise ValueError("Too many values")
                metrics = {}
                for j in range(len(values)):
                    value = values[j]
                    if re.match(r"^[0-9]+\.[0-9]+$", value):
                        metrics[keys[j]] = float(value)
                    elif re.match(r"^[0-9]+$", value):
                        metrics[keys[j]] = int(value)
                    elif value == "":
                        metrics[keys[j]] = 0.0
                    else:
                        metrics[keys[j]] = value
                ''' If no strandedness parameters, set "correct strand" metrics to null'''
                if not self.strand_specific:
                    for key in self.STRAND_SPECIFIC_METRICS:
                        metrics[key] = None
                output['metrics'] = metrics
            elif re.match('## HISTOGRAM', line):
                ''' histogram of coverage'''
                tcu = {'headers': re.split(r'\s+', lines[i + 1].strip()), 'values': []}
                j = 2
                while re.match('[0-9]', lines[i + j]):
                    fields = re.split(r'\s+', lines[i + j].strip())
                    if len(fields) != 2:
                        raise ValueError("Incorrect number of fields")
                    tcu['values'].append([int(fields[0]), float(fields[1])])
                    j += 1
                output['normalized_coverage_histogram'] = tcu
            output['header'] = "\n".join(header)
        return output

    @staticmethod
    def parse_unique_reads(input_path):
        return int(open(input_path).read().strip())

    def write_json(self, out_path):
        if out_path == '-':
            out = sys.stdout
        else:
            out = open(out_path, 'w')
        output = {}
        output[self.BAMQC_KEY] = self.bamqc_metrics
        output[self.PICARD_KEY] = self.picard_metrics
        output[self.CONTAM_KEY] = self.contam_metrics
        output[self.UNIQUE_READS_KEY] = self.unique_reads
        out.write(json.dumps(output, sort_keys=True))
        if out_path != '-':
            out.close()


def main():
    parser = argparse.ArgumentParser(description='Write JSON summary of RNASeqQC output')
    parser.add_argument('--bamqc', help="Path to bamqc JSON input", metavar="PATH", required=True)
    parser.add_argument('--contam', help="Path to ribosome contamination text input", metavar="PATH", required=True)
    parser.add_argument('--picard', help="Path to Picard CollectRnaSeqMetrics text input",
                        metavar="PATH", required=True)
    parser.add_argument('--unique-reads', help="Path to unique reads text input", metavar="PATH", required=True)
    parser.add_argument('--no-strand-specificity', help='Replace Picard strand metrics with null values',
                        action="store_true")
    parser.add_argument('--out', help='Path for JSON output, or - for STDOUT', required=True)
    args = parser.parse_args()
    config = {
        collator.BAMQC_KEY: args.bamqc,
        collator.CONTAM_KEY: args.contam,
        collator.PICARD_KEY: args.picard,
        collator.NO_STRAND_SPECIFICITY_KEY: args.no_strand_specificity,
        collator.UNIQUE_READS_KEY: args.unique_reads
    }
    collator(config).write_json(args.out)


if __name__ == "__main__":
    main()
