# take multiple mosdepth summary files and get results in suitable format

import sys

files = sys.argv[1:]

for fn in files:
    fns = fn.split('.')
    with open(fn) as f:
        for x in f.readlines():
            x = x.rstrip()
            v = x.split('\t')
            if v[0] == 'chrM':
                print(fns[0], v[3], sep='\t')
                # To make the format fit with multiqc output from sarek
                print("\n\n")
