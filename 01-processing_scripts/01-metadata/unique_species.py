import pandas as pd
import random

df = pd.read_csv("metadata.tsv", sep="\t")

# Select rows where rpoC_gca is 1
filtered_df = df[df['rpoC_gca'] == 1]
filtered_df2 = filtered_df[filtered_df['kingdom'] == "Bacteria"]

# Get stats
length_sd = filtered_df2['length'].std()
length_mean = filtered_df2['length'].mean()

# Define the bounds
lower_bound = length_mean - length_sd
upper_bound = length_mean + length_sd

filtered_df3 = filtered_df2[(filtered_df2['length'] >= lower_bound) & (filtered_df2['length'] <= upper_bound)]
filtered_df4 = filtered_df3[~filtered_df3['species'].str.contains('_sp[1234567890]')]

# Remove anything after the final underscore if it ends with a capital letter
def remove_suffix_if_capital(species):
    if species[-1].isupper():  # Check if it ends with a capital letter
        parts = species.rsplit('_', 1)  # Split at the last underscore
        return parts[0]  # Return the part before the last underscore
    return species

filtered_df4['species'] = filtered_df4['species'].apply(remove_suffix_if_capital)

# Create the edited species column
filtered_df5 = filtered_df4[~filtered_df4['species'].str.contains('_[A-Z]{2}')]
filtered_df5['species_edits'] = filtered_df5['species'].str.replace(r'_[A-Z]', '', regex=True)

# Resetting index for the species DataFrame to ensure proper alignment if needed
unique_samples = filtered_df5.groupby('species_edits').sample(n=1, random_state=1).reset_index(drop=True)

# Save the results to a new TSV file
unique_samples.to_csv('unique_species.tsv', index=False, sep="\t")
