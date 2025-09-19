



import pandas as pd

# Load BigSolDB raw file
bigsol = pd.read_csv("../data_raw/BigSolDBv2.0.csv")

# Step 1: Keep only important columns
useful_cols = ["SMILES_Solute", "Solvent", "Temperature_K", "Solubility(mol/L)"]
bigsol_clean = bigsol[useful_cols].copy()

# Step 2: Rename for consistency
bigsol_clean = bigsol_clean.rename(columns={
    "SMILES_Solute": "Solute",
    "Solvent": "Solvent",
    "Temperature_K": "Temperature_K",
    "Solubility(mol/L)": "Solubility"
})

# Step 3: Add a units column
bigsol_clean["Units"] = "mol/L"

# Step 4: Drop missing values
bigsol_clean = bigsol_clean.dropna()

# Step 5: Quick check
print("Clean dataset shape:", bigsol_clean.shape)
print(bigsol_clean.head())

# Step 6: Save cleaned file
bigsol_clean.to_csv("../data_processed/BigSolDB_clean.csv", index=False)
print("✅ Cleaned dataset saved as BigSolDB_clean.csv")





# import pandas as pd

# # Load the raw BigSolDB file
# bigsol = pd.read_csv("../data_raw/BigSolDBv2.0.csv")

# # Print column names
# print("Columns in BigSolDB:", bigsol.columns.tolist())

# # Show first few rows
# print(bigsol.head())







