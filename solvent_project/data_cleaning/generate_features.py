import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

# Ensure output folder exists
os.makedirs("../features", exist_ok=True)

# Helper function to compute descriptors
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol)
    }

# Process a dataset
def process_dataset(input_path, output_path):
    df = pd.read_csv(input_path)

    features = []
    for smi in df["Solute"]:
        desc = compute_descriptors(smi)
        if desc:
            features.append(desc)
        else:
            features.append({"MW": None, "LogP": None, "TPSA": None, "HBD": None, "HBA": None})

    feat_df = pd.DataFrame(features)
    result = pd.concat([df, feat_df], axis=1)

    result.to_csv(output_path, index=False)
    print(f"✅ Saved features to {output_path} with shape {result.shape}")

# Run for both datasets
process_dataset("data_processed/AqSolDB_clean.csv", "features/AqSolDB_features.csv")
process_dataset("data_processed/BigSolDB_clean.csv", "features/BigSolDB_features.csv")

