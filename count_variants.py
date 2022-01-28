#!/usr/bin/env python3

# Count variants from mito_var Nextflow script to determine frequency/type of variants
# 1. Exclude control region (>15120bp)
# 2. Include only PASS variants

import sys
import gzip

# Point at which control region begins and from which we will ignore variants until end of sequence
control_begin = 15120
# Ignore alleles above this frequencyas presumably they were fixed prior to the experiment
fixed_freq = 0.9

def read_file(fn):
    '''Generator for reading files either zipped or not'''
    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            for line in f:
                yield line
    else:
        with open(fn) as f:
            for line in f:
                yield line

def read_vcf (vcf):
    ''' Generator to read each non-header line of a vcf in as a list of elements'''
    for x in read_file(vcf):
        if x.startswith('#'):
            continue
        else:
            x = x.rstrip()
            v = x.split('\t')
            yield(v)

def count_variants (vcf_file):
    c = 0
    snps = 0
    indels = 0
    for x in read_vcf(vcf_file):
        sample_field = x[-1].split(':')
        # There can be multiple frequency values where there are multiple alternative alleles
        freq_vals = sample_field[2].split(',')

        # Skip filtered variants
        if x[6] != 'PASS':
            continue
        # Exclude variants from control region
        elif int(x[1]) >= control_begin:
            continue
        # Exclude variants unless there is a low frequency one

        min_freq = min(freq_vals)
        if float(min_freq) > fixed_freq:
            continue
        
        # Assume this is kept in because it hasn't been skipped
        #print(x)
        c += 1

        # Is it a SNP or indel?
        if len(x[3]) > 1 or len(x[4]) > 1:
            indels += 1
        else:
            snps += 1
    return c, snps, indels

###############
##### MAIN

ss = sys.argv[1]

with open(ss, 'r') as s:
    for x in s.readlines():
        x = x.rstrip()
        v = x.split(',')
        if v[0].startswith('sample'):
            header = v
        else:
            #print(v[-1])
            c, snps, indels = count_variants(v[-1])
            print ("{} {} {} {} {} {} {} {}".format(v[0], v[1], v[2], v[3], v[4], c, snps, indels))
