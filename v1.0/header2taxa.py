import pandas as pd
df1 = pd.read_table("header2taxid",delimiter='\t', header= None)
df1.columns =['header', 'taxid']
print(df1.head())

df2 = pd.read_table("ids_taxonomy.txt",delimiter='\t', header= None)
df2.columns =['taxid', 'taxa']
print(df2.head())

combined = pd.merge(df1, df2, on ='taxid')

print(combined['taxid'].unique())

combined.to_csv('header_taxa.tsv', sep="\t")