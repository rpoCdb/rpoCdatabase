import pandas as pd

# Load the list of genera from genomes.txt
with open('pom_genera.txt', 'r') as file:
    filters = set(line.strip() for line in file if line.strip())

# Load the unique_species.tsv file into a DataFrame
df = pd.read_csv('unique_species.tsv', sep='\t')

# Filter the DataFrame based on the genera
filtered_df = df[df['genus'].isin(filters) | df['family'].isin(filters) | df['order'].isin(filters)]

# Output the filtered DataFrame to a new TSV file
filtered_df.to_csv('filtered_species.tsv', sep='\t', index=False)

print("Filtered data has been saved to 'filtered_species.tsv'.")
