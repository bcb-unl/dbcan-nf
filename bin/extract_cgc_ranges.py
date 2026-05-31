#!/usr/bin/env python

import sys
import pandas as pd

# arguments: cgc_standard_out.tsv, substrate_prediction.tsv, output.tsv
cgc_tsv = sys.argv[1]
substrate_tsv = sys.argv[2]
out_tsv = sys.argv[3]

# read substrate_prediction.tsv and extract CGCID set
substrate = pd.read_csv(substrate_tsv, sep='\t')
substrate = substrate[substrate['PULID'].notna()]
cgcid_set = set(substrate['#cgcid'].astype(str))

# read cgc_standard_out.tsv
df = pd.read_csv(cgc_tsv, sep='\t')
df['CGCID'] = df['Contig ID'].astype(str) + '|' + df['CGC#'].astype(str)

#   only keep CGCIDs that are in the substrate set
df = df[df['CGCID'].isin(cgcid_set)]

# For each CGCID, get contig span from first/last gene (by Gene Start)
ranges = (
    df.sort_values(['CGCID', 'Gene Start'])
      .groupby('CGCID', as_index=False)
      .agg({
          'Contig ID': 'first',
          'CGC#': 'first',
          'Gene Start': 'first',
          'Gene Stop': 'last',
      })
)

ranges.to_csv(out_tsv, sep='\t', index=False, header=True)
