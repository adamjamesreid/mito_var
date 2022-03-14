#!/usr/bin/env python3

# Count variants from mito_var Nextflow script to determine frequency/type of variants
# 1. Exclude control region (>15120bp)
# 2. Include only PASS variants

import sys
import gzip
import subprocess
import re

# Point at which control region begins and from which we will ignore variants until end of sequence
control_begin = 15120
# Ignore alleles above this frequencyas presumably they were fixed prior to the experiment
max_allele_freq = 0.9

# Minimum allele frequency
min_allele_freq = 1e-4

# Allow these AS_FilterStatus tags through (i.e. call is either PASS or these tags)
#['base_qual', 'contamination', 'map_qual', 'position', 'strand_bias', 'weak_evidence']
# n.b. SITE is equivalent to PASS
#allow_these_tags = ['SITE', 'strand_bias', 'base_qual', 'map_qual', 'weak_evidence']
allow_these_tags = ['SITE', 'strand_bias']


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

def read_vcf_header (vcf):
    ''' Generater to return header of a VCF file'''
    for x in read_file(vcf):
        if x.startswith('#'):
            yield(x)

def filter_vcf (vcf_file, control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
    ''' Generator to return lines of filtered VCF file suitable for rare variant calling'''
    for x in read_vcf(vcf_file):
        sample_field = x[-1].split(':')
        # There can be multiple frequency values where there are multiple alternative alleles
        freq_vals = sample_field[2].split(',')

        # Yield variant entry unless it should be filtered
        if call_filter(x, control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
            yield('\t'.join(x))
        

def call_filter (vcf_line, control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
    ''' Return True if variant passes filters '''
    # Exclude variants from control region
    if int(vcf_line[1]) >= control_begin:
        #print ("vcf_line[1] in control region")
        return False

    # Get frequency values for each allele
    # There can be multiple frequency values where there are multiple alternative alleles
    sample_field = vcf_line[-1].split(':')
    freq_vals = sample_field[2].split(',')

    # Look for a set of AS_filter tags which pass
    m = re.findall(r'AS_FilterStatus=([^;]+)', vcf_line[7]) 
    # Split status field per allele
    status_fields = m[0].split('|')
    # Loop through allele status and frequencies checking whether each passes
    for s, f in zip(status_fields, freq_vals):
        passed = True
        tags = s.split(',')
        # Loop through each tag to see if it is OK
        for t in tags:
            # set passed to false if t is not an allowed tag
            if t not in allow_these_tags:
                passed = False
                #print ("{} has tag {}".format(vcf_line[1], t))
        # Check frequency
        if passed:
            # IF frequnecy is in appropriate range, return True
            if float(f) >= min_allele_freq and float(f) <= max_allele_freq:
                return True
            else:
                pass
                #print ("{} has freq {}".format(vcf_line[1], f))
    # If filters are not passed, return False
    return False

def count_variants (vcf_file, control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
    ''' Generator to return counts of filtered snps and indels suitable for rare variant calling'''
    c = 0
    snps = 0
    indels = 0
    for x in read_vcf(vcf_file):
        sample_field = x[-1].split(':')
        # There can be multiple frequency values where there are multiple alternative alleles
        freq_vals = sample_field[2].split(',')

        # Check whether this variant passes filters
        if call_filter(x, control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
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
            c, snps, indels = count_variants(v[-1], control_begin, max_allele_freq, min_allele_freq, allow_these_tags)
            print ("{} {} {} {} {} {} {} {}".format(v[0], v[1], v[2], v[3], v[4], c, snps, indels))

            # Write out filtered VCF
            out_file = v[-1] + '_filtered.vcf'
            o = open(out_file, 'w')
            for i in read_vcf_header(v[-1]):
                o.write(i)
            for j in filter_vcf(v[-1], control_begin, max_allele_freq, min_allele_freq, allow_these_tags):
                o.write(j + '\n')
            o.close()
            zip_cmd = "bgzip {}".format(v[-1] + '_filtered.vcf')
            subprocess.call(zip_cmd, shell=True)
            tabix_cmd = "tabix -p vcf {}".format(v[-1] + '_filtered.vcf.gz')
            subprocess.call(tabix_cmd, shell=True)

