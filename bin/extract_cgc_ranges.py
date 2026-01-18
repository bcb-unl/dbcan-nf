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

# for each CGCID, get the range of gene start and stop positions
def get_cgc_range(group):
    group_sorted = group.sort_values('Gene Start')
    first_gene = group_sorted.iloc[0]
    last_gene = group_sorted.iloc[-1]
    return pd.Series({
        'Contig ID': first_gene['Contig ID'],
        'CGC#': first_gene['CGC#'],
        'Gene Start': first_gene['Gene Start'],
        'Gene Stop': last_gene['Gene Stop']
    })

ranges = (
    df.groupby('CGCID', group_keys=False)
      .apply(get_cgc_range)
      .reset_index()
)

ranges.to_csv(out_tsv, sep='\t', index=False, header=True)
