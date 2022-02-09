# bind datasets together

import pandas as pd
import sys

counts_file = sys.argv[1]
cov_file = sys.argv[2]

cf = pd.read_csv(counts_file, sep='\t', index_col=1)
covf = pd.read_csv(cov_file, sep='\t', index_col=0)

combined_df = pd.concat([cf, covf], axis=1)

combined_df.to_csv("counts_coverage_combined.tsv", sep='\t')
