import pandas as pd

df = pd.read_table("rpoC_cords",delimiter='\t', header= None)

df.columns =['seqid', 'cord1', 'cord2']

s = df['cord1'] > df['cord2']
df.loc[s, ['cord1','cord2']] = df.loc[s, ['cord2','cord1']].values
df['length'] = df.cord2 - df.cord1
df.to_csv('flip_cords.tsv', sep="\t")