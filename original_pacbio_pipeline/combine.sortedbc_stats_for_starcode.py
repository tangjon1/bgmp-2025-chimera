import pandas as pd
from glob import glob

# List of file names
file_names = [
    "part_00.noEnvZbc_stats_for_starcode.tsv",
    "part_01.noEnvZbc_stats_for_starcode.tsv",
    "part_02.noEnvZbc_stats_for_starcode.tsv",
    "part_03.noEnvZbc_stats_for_starcode.tsv",
    "part_04.noEnvZbc_stats_for_starcode.tsv"
]

# Initialize an empty dataframe
merged_df = pd.DataFrame()

# Iterate over each file
for file_name in file_names:
    df = pd.read_csv(file_name, sep='\t', header=None, names=['barcode', 'count'])
    merged_df = pd.concat([merged_df, df])

# Group by barcode and sum the counts
merged_df = merged_df.groupby('barcode').sum().reset_index()

# Save the merged dataframe to a new TSV file
merged_df.to_csv('merged.sortedbc_stats_for_starcode.tsv', sep='\t', index=False, header=False)

print("Files merged successfully!")
