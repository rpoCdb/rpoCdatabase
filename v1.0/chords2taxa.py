import pandas as pd
df1 = pd.read_table("cords.tsv",delimiter='\t')

df2 = pd.read_table("accs2taxa",delimiter='\t')
df2.columns =['GCA', 'taxa']
print(df2.head())

df_all = pd.merge(df1, df2, on ='GCA', how='left', indicator=True)
needed_seq = df_all[df_all['_merge'] == 'left_only']

needed_seq.to_csv('get_seqs.tsv', sep="\t", index=False, header = False)

df_all.to_csv('taxa.tsv', sep="\t", index=False)