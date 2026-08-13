import os
import json
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from difflib import get_close_matches
import pubchempy as pcp
import sys

# ----------------------------
# Path handling 
# ----------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODELS_DIR = os.path.join(BASE_DIR, "models")
MODEL_PATH = "d:/AI_projects/DissolveAI/solvent_project/models/solubility_model_bigsol.pkl"
ENCODER_PATH = "d:/AI_projects/DissolveAI/solvent_project/models/solvent_encoder_bigsol.pkl"
METADATA_PATH = os.path.join(MODELS_DIR, "metadata.json")

TRAIN_CSV_PATHS = [
    "d:/AI_projects/DissolveAI/solvent_project/data_processed/BigSolDB_clean.csv",
    "d:/AI_projects/DissolveAI/solvent_project/data_processed/AqSolDB_clean.csv"
]

# ----------------------------
# create metadata from CSVs if metadata.json missing
# ----------------------------
def build_and_save_metadata(paths, metadata_path):
    temp_min, temp_max = 1e9, -1e9
    solvents = set()
    found_any = False
    for p in paths:
        if os.path.exists(p):
            found_any = True
            try:
                df = pd.read_csv(p)
                if "Temperature_K" in df.columns:
                    tmin = pd.to_numeric(df["Temperature_K"], errors="coerce").min()
                    tmax = pd.to_numeric(df["Temperature_K"], errors="coerce").max()
                    if pd.notna(tmin):
                        temp_min = min(temp_min, float(tmin))
                    if pd.notna(tmax):
                        temp_max = max(temp_max, float(tmax))
                if "Solvent" in df.columns:
                    solvents.update(df["Solvent"].dropna().astype(str).unique().tolist())
            except Exception:
                continue

    if not found_any or temp_min == 1e9 or temp_max == -1e9:
        temp_min, temp_max = 273.0, 373.0

    metadata = {
        "temperature_min": float(temp_min),
        "temperature_max": float(temp_max),
        "solubility_lower_bound": -12.0,
        "solubility_upper_bound": 2.0,
        "solubility_threshold": -2.0,
        "total_solvents": len(solvents),
        "unique_solvents_example": list(list(solvents)[:10])
    }

    os.makedirs(os.path.dirname(metadata_path), exist_ok=True)
    with open(metadata_path, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)
    print(f"Created metadata at {metadata_path}")
    return metadata

# ----------------------------
# Load model / encoder / metadata
# ----------------------------
if not os.path.exists(MODEL_PATH) or not os.path.exists(ENCODER_PATH):
    print("Error: Trained model or encoder not found.")
    print(f"Expected model at: {MODEL_PATH}")
    print(f"Expected encoder at: {ENCODER_PATH}")
    print("Please run the training script (train_full_pipeline_bigsol.py) first.")
    sys.exit(1)

model = joblib.load(MODEL_PATH)
encoder = joblib.load(ENCODER_PATH)

if os.path.exists(METADATA_PATH):
    with open(METADATA_PATH, "r", encoding="utf-8") as fh:
        metadata = json.load(fh)
else:
    print("metadata.json not found — attempting to build metadata from CSV files...")
    metadata = build_and_save_metadata(TRAIN_CSV_PATHS, METADATA_PATH)

temp_min = metadata.get("temperature_min", 273.0)
temp_max = metadata.get("temperature_max", 373.0)
SOLUBILITY_LOWER_BOUND = metadata.get("solubility_lower_bound", -12.0)
SOLUBILITY_UPPER_BOUND = metadata.get("solubility_upper_bound", 2.0)
SOLUBILITY_THRESHOLD = metadata.get("solubility_threshold", -2.0)

# ----------------------------
# Descriptor helpers 
# ----------------------------
DESCRIPTOR_FUNCS = {
    "MW": Descriptors.MolWt,
    "LogP": Descriptors.MolLogP,
    "TPSA": Descriptors.TPSA,
}

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

def get_compound_name(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid compound"
    try:
        res = pcp.get_compounds(smiles, "smiles")
        if res and getattr(res[0], "iupac_name", None):
            return res[0].iupac_name
        if res and getattr(res[0], "synonyms", None) and res[0].synonyms:
            return res[0].synonyms[0]
    except Exception:
        pass
    return Chem.MolToSmiles(mol, canonical=True)

def find_encoder_match(solvent_input):
    if solvent_input is None:
        return None
    categories = encoder.categories_[0]
    norm_input = str(solvent_input).strip().lower()
    for s in categories:
        if s.lower() == norm_input:
            return s
    close = get_close_matches(norm_input, [s.lower() for s in categories], n=1, cutoff=0.6)
    if close:
        for s in categories:
            if s.lower() == close[0]:
                return s
    # last resort: try containment
    for s in categories:
        if norm_input in s.lower():
            return s
    return None

def build_feature_row(temperature, mw, logp, tpsa, solvent_name):
    solvent_encoded = encoder.transform(pd.DataFrame({"Solvent": [solvent_name]}))
    feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(["Solvent"]).tolist()
    row = [temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]
    return pd.DataFrame([row], columns=feature_names)

# ----------------------------
# CLI Prediction  
# ----------------------------
if __name__ == "__main__":
    print("Solubility Prediction Tool (robust paths + metadata fallback)")
    solute_smiles = input("Enter solute SMILES: ").strip()
    solvent_input = input("Enter solvent name: ").strip()
    temp_raw = input("Enter temperature in Kelvin [default=298]: ").strip()

    desc = compute_descriptors(solute_smiles)
    if desc is None:
        print("Invalid SMILES string.")
        sys.exit(1)
    mw, logp, tpsa = desc
    solute_name = get_compound_name(solute_smiles)

    try:
        temperature = float(temp_raw) if temp_raw else 298.0
    except Exception:
        temperature = 298.0

    clipped_note = None
    if temperature < temp_min:
        temperature, clipped_note = temp_min, f"(Note: clamped to {temp_min} K)"
    if temperature > temp_max:
        temperature, clipped_note = temp_max, f"(Note: clamped to {temp_max} K)"

    matched = find_encoder_match(solvent_input)
    if matched is None:
        print(f"No data for solvent '{solvent_input}'. Try one of: {list(encoder.categories_[0])[:12]}")
        sys.exit(1)

    X = build_feature_row(temperature, mw, logp, tpsa, matched)
    pred_raw = float(model.predict(X)[0])
    pred = max(min(pred_raw, SOLUBILITY_UPPER_BOUND), SOLUBILITY_LOWER_BOUND)

    print("\n================= Prediction =================")
    print(f"Compound: {solute_name}")
    print(f"Solvent : {matched}")
    print(f"Temp    : {temperature} K {clipped_note or ''}")
    print(f"Predicted solubility: {pred:.3f} log mol/L")

    if pred > SOLUBILITY_THRESHOLD:
        print("Good solubility predicted.")
    else:
        print("Poor solubility in this solvent. Suggesting alternatives...")
        suggestions = []
        for alt in encoder.categories_[0]:
            if alt == matched:
                continue
            try:
                X_alt = build_feature_row(temperature, mw, logp, tpsa, alt)
                alt_pred = float(model.predict(X_alt)[0])
                alt_pred = max(min(alt_pred, SOLUBILITY_UPPER_BOUND), SOLUBILITY_LOWER_BOUND)
                if alt_pred > SOLUBILITY_THRESHOLD:
                    suggestions.append((alt, round(alt_pred, 3)))
            except Exception:
                continue
        suggestions = sorted(suggestions, key=lambda x: -x[1])
        if suggestions:
            print("Alternative solvents with better solubility:")
            for s, val in suggestions[:8]:
                print(f"  - {s}: {val} log mol/L")
        else:
            print("No better solvents found.")
