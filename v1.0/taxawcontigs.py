import pandas as pd
df1 = pd.read_table("length_taxa2",delimiter='\t', header= None)
df1.columns =['contig', 'gca', 'assembely_level', 'cord1', 'cord2', 'length', 'kingdom','phylum', 'class', 'order', 'family', 'genus', 'species']
print(df1.head())

df2 = pd.read_table("counts",delimiter='\t', header= None)
df2.columns =['gca','rpoC_gca', 'contig_number']
print(df2.head())

df3 = pd.read_table("contig_rpoC",delimiter='\t', header= None)
df3.columns =['contig', 'rpoC_num']

combined = pd.merge(df1, df2, on ='gca')
combined2 = pd.merge(combined, df3, on ='contig')

combined2.to_csv('metadata.tsv', sep="\t", index =False)