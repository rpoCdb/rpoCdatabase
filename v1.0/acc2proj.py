import pandas as pd
df1 = pd.read_table("missing_taxa",delimiter='\t', header=None)
df1.columns =['GCA']

df2 = pd.read_table("acc2proj",delimiter='\t', header=None)
df2.columns =['GCA', 'proj']
print(df2.head())

df_all = pd.merge(df1, df2, on ='GCA')

df_all.to_csv('projs.tsv', sep="\t", index=False, header = False)