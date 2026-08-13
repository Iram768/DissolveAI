# import pandas as pd

# files = [
#     "dataset-A.csv", "dataset-B.csv", "dataset-C.csv",
#     "dataset-D.csv", "dataset-E.csv", "dataset-F.csv",
#     "dataset-G.csv", "dataset-H.csv", "dataset-I.csv"
# ]

# # Read and merge all CSVs
# df_list = [pd.read_csv(file) for file in files]
# merged_df = pd.concat(df_list, ignore_index=True)

# # Show info
# print("Merged dataset shape:", merged_df.shape)
# print("Columns:", merged_df.columns)

# # Save merged dataset
# merged_df.to_csv("AqSolDB_merged.csv", index=False)
# print("✅ Merged file saved as AqSolDB_merged.csv")






import pandas as pd
import os

folder = "AqSolDB all files"

files = [
    "dataset-A.csv", "dataset-B.csv", "dataset-C.csv",
    "dataset-D.csv", "dataset-E.csv", "dataset-F.csv",
    "dataset-G.csv", "dataset-H.csv", "dataset-I.csv"
]

file_paths = [os.path.join(folder, f) for f in files]
df_list = [pd.read_csv(path) for path in file_paths]
merged_df = pd.concat(df_list, ignore_index=True)

print("Merged dataset shape:", merged_df.shape)
print("Columns:", merged_df.columns)

merged_df.to_csv("AqSolDB_merged.csv", index=False)

print("✅ Merged file saved as AqSolDB_merged.csv")





