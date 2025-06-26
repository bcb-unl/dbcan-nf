#!/usr/bin/env python

import sys
import pandas as pd

# arguments: cgc_standard_out.tsv, substrate_prediction.tsv, output.tsv
cgc_tsv = sys.argv[1]
substrate_tsv = sys.argv[2]
out_tsv = sys.argv[3]

# read substrate_prediction.tsv and extract CGCID set
substrate = pd.read_csv(substrate_tsv, sep='\t')
#print(substrate)
cgcid_set = set(substrate['#cgcid'].astype(str))
#print(cgcid_set)

# read cgc_standard_out.tsv
df = pd.read_csv(cgc_tsv, sep='\t')
df['CGCID'] = df['Contig ID'].astype(str) + '|' + df['CGC#'].astype(str)
df = df.drop_duplicates(subset=['CGCID'])

#print(df['CGCID'])

# only keep rows where CGCID is in the substrate set
df = df[df['CGCID'].isin(cgcid_set)]

# group by CGCID and get the first Contig ID, min Gene Start, and max Gene Stop
ranges = df.groupby('CGCID').agg({
    'Contig ID': 'first',
    'CGC#': 'first',
    'Gene Start': 'min',
    'Gene Stop': 'max'
}).reset_index()

ranges.to_csv(out_tsv, sep='\t', index=False, header=True)
