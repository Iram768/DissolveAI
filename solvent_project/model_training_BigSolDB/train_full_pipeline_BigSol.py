import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, DataStructs
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
import joblib

# ----------------------------
# Paths
# ----------------------------
DATA_PATH = "d:/AI_projects/DissolveAI/solvent_project/data_processed/BigSolDB_clean.csv"
MODEL_PATH = "d:/AI_projects/DissolveAI/solvent_project/models/solubility_model_bigsol.pkl"
ENCODER_PATH = "d:/AI_projects/DissolveAI/solvent_project/models/solvent_encoder_bigsol.pkl"
METADATA_PATH = "d:/AI_projects/DissolveAI/solvent_project/models/metadata.json"

# ----------------------------
# Feature engineering 
# ----------------------------
def compute_descriptors(smiles: str):
    """Compute RDKit descriptors for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "FracCSP3": Descriptors.FractionCSP3(mol)
    }

from rdkit.Chem import rdFingerprintGenerator

def compute_fingerprints(smiles, n_bits=1024, radius=2):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.GetFingerprint(mol)
    return list(fp)

# ----------------------------
# Load data
# ----------------------------
df = pd.read_csv(DATA_PATH)
print(f"Loaded {len(df)} rows from {DATA_PATH}")

# ----------------------------
# Compute features
# ----------------------------
desc_list, fps_list, solvents, y = [], [], [], []

for idx, row in df.iterrows():
    smi = str(row["Solute"])  
    desc = compute_descriptors(smi)
    if desc is None:
        continue
    fp = compute_fingerprints(smi, n_bits=1024)

    desc_list.append(desc)
    fps_list.append(fp)
    solvents.append(row["Solvent"])
    y.append(float(row["Solubility"]))

desc_df = pd.DataFrame(desc_list)
fps_arr = np.array(fps_list)
y = np.array(y)

# ----------------------------
# Encode solvents
# ----------------------------
encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
solvent_encoded = encoder.fit_transform(np.array(solvents).reshape(-1, 1))
solvent_cols = encoder.get_feature_names_out(["Solvent"])
solvent_df = pd.DataFrame(solvent_encoded, columns=solvent_cols)

# ----------------------------
# Combine all features
# ----------------------------
X = np.hstack([desc_df.values, fps_arr, solvent_encoded])

X = X.astype(np.float32)
y = y.astype(np.float32)

print(f"Feature matrix shape: {X.shape}, Target shape: {y.shape}")
print(f"X dtype: {X.dtype}, y dtype: {y.dtype}")
# ----------------------------
# Train/test split
# ----------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# ----------------------------
# Train model
# ----------------------------
model = RandomForestRegressor(
    n_estimators=500,
    n_jobs=-1,
    random_state=42
)
model.fit(X_train, y_train)

# ----------------------------
# Evaluate
# ----------------------------
y_pred = model.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print(f"Test RMSE: {rmse:.3f}")

cv_scores = -cross_val_score(
    model, X, y, cv=5,
    scoring="neg_mean_squared_error",
    n_jobs=1
)
print(f"5-fold CV RMSE: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

# ----------------------------
# Save model + encoder + metadata
# ----------------------------
os.makedirs("models", exist_ok=True)
joblib.dump(model, MODEL_PATH)
joblib.dump(encoder, ENCODER_PATH)

metadata = {
    "features": {
        "descriptors": list(desc_df.columns),
        "fingerprint_bits": fps_arr.shape[1],
        "solvent_columns": solvent_cols.tolist()
    },
    "target": "Solubility",
    "n_samples": len(y),
    "rmse": float(rmse)
}
with open(METADATA_PATH, "w") as f:
    json.dump(metadata, f, indent=2)

print("✅ Model, encoder, and metadata saved successfully!")
