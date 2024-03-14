import pandas as pd
df1 = pd.read_table("seqids_needed",delimiter='\t', header= None)
df1.columns =['seq']
print(df1.head())

df2 = pd.read_table("annot_cords",delimiter='\t', header= None)
df2.columns =['seq']
print(df2.head())

df_all = pd.merge(df1, df2, on ='seq', how='left', indicator=True)
needed_seq = df_all[df_all['_merge'] == 'left_only']

needed_seq.to_csv('get_seqs.tsv', sep="\t", index=False, header = False)