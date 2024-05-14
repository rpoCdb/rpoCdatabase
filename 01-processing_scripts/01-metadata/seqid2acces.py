import pandas as pd
df1 = pd.read_table("flip_cords.tsv",delimiter='\t', header= None)
df1.columns =['seq', 'cords1','cord2', 'length']
print(df1.head())

df2 = pd.read_table("acces2seqid",delimiter='\t', header= None)
df2.columns =['GCA', 'seq']
print(df2.head())

combined = pd.merge(df1, df2, on ='seq')

print(combined['GCA'].unique())

combined.to_csv('cords.tsv', sep="\t", index=False)