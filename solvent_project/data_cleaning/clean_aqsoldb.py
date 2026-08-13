

import pandas as pd
import os

# Load your merged AqSolDB
aqsol = pd.read_csv(r"D:\AI_projects\DissolveAI\solvent_project\AqSolDB all files\AqSolDB_merged.csv")

# Show columns so you know what's inside
print("Columns in AqSolDB:", aqsol.columns)

# Standardize schema:
# AqSolDB gives solute SMILES and solubility in mol/L already.
# Solvent is always water, temperature usually 25 °C (298 K).
aqsol_clean = pd.DataFrame({
    "Solute": aqsol["SMILES"],  
    "Solvent": "water",  
    "Temperature_K": 298,  
    "Solubility": aqsol["Solubility"],  
    "Units": "mol/L"  
})

print("Clean dataset shape:", aqsol_clean.shape)
print(aqsol_clean.head())

# Save inside data_processed
os.makedirs("../data_processed", exist_ok=True)
aqsol_clean.to_csv("../data_processed/AqSolDB_clean.csv", index=False)

print("✅ Cleaned dataset saved as AqSolDB_clean.csv")















