from flask import Flask, request, jsonify
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from flask_cors import CORS
from difflib import get_close_matches
import pubchempy as pcp
import json
import os

app = Flask(__name__)
CORS(app)

# -------------------
# Paths
# -------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, ".."))

MODEL_PATH = os.path.join(PROJECT_ROOT, "models", "solubility_model_bigsol.pkl")
ENCODER_PATH = os.path.join(PROJECT_ROOT, "models", "solvent_encoder_bigsol.pkl")
METADATA_PATH = os.path.join(PROJECT_ROOT, "models", "metadata.json")

# -------------------
# Load model/encoder/metadata
# -------------------
model = joblib.load(MODEL_PATH)
encoder = joblib.load(ENCODER_PATH)

with open(METADATA_PATH, "r") as f:
    metadata = json.load(f)

temp_min = metadata["temperature_min"]
temp_max = metadata["temperature_max"]
SOLUBILITY_LOWER_BOUND = metadata["solubility_lower_bound"]
SOLUBILITY_UPPER_BOUND = metadata["solubility_upper_bound"]
SOLUBILITY_THRESHOLD = metadata["solubility_threshold"]

# -------------------
# Helpers
# -------------------
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
    except Exception:
        pass
    return Chem.MolToSmiles(mol, canonical=True)

def find_encoder_match(solvent_input):
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
    return None

def build_feature_row(temperature, mw, logp, tpsa, solvent_name):
    solvent_encoded = encoder.transform(pd.DataFrame({"Solvent": [solvent_name]}))
    feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(["Solvent"]).tolist()
    row = [temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]
    return pd.DataFrame([row], columns=feature_names)

# -------------------
# Main route
# -------------------
@app.route("/predict", methods=["POST"])
def predict():
    try:
        data = request.get_json() or {}
        solute_smiles = data.get("solute_smiles")
        solvent_input = data.get("solvent_name")
        temp_raw = data.get("temperature", None)

        if not solute_smiles:
            return jsonify({"error": "No solute SMILES provided."}), 400

        desc = compute_descriptors(solute_smiles)
        if desc is None:
            return jsonify({"error": "Invalid SMILES string provided."}), 400
        mw, logp, tpsa = desc

        solute_name = get_compound_name(solute_smiles)

        try:
            temperature = float(temp_raw) if temp_raw else 298.0
        except Exception:
            temperature = 298.0

        clipped_note = None
        if temperature < temp_min:
            temperature, clipped_note = temp_min, f"Temperature was clamped to {temp_min} K"
        if temperature > temp_max:
            temperature, clipped_note = temp_max, f"Temperature was clamped to {temp_max} K"

        matched = find_encoder_match(solvent_input)
        if matched is None:
            return jsonify({
                "solute_name": solute_name,
                "message": f"No data available for solvent '{solvent_input}'. Please try a valid solvent name."
            }), 200

        X = build_feature_row(temperature, mw, logp, tpsa, matched)
        pred_raw = float(model.predict(X)[0])
        pred = max(min(pred_raw, SOLUBILITY_UPPER_BOUND), SOLUBILITY_LOWER_BOUND)

        response = {
            "solute_name": solute_name,
            "solvent_name": matched,
            "temperature_used_K": temperature,
            "predicted_solubility": round(pred, 3)
        }
        if clipped_note:
            response["note"] = clipped_note

        if pred > SOLUBILITY_THRESHOLD:
            response["message"] = f"The solubility of {solute_name} in {matched} at {temperature} K is {pred:.3f} log mol/L."
            return jsonify(response), 200

        # Poor solubility → suggest better ones
        suggestions = []
        for alt in encoder.categories_[0]:
            if alt == matched:
                continue
            try:
                X_alt = build_feature_row(temperature, mw, logp, tpsa, alt)
                alt_pred = float(model.predict(X_alt)[0])
                alt_pred = max(min(alt_pred, SOLUBILITY_UPPER_BOUND), SOLUBILITY_LOWER_BOUND)
                if alt_pred > SOLUBILITY_THRESHOLD:
                    suggestions.append({"solvent": alt, "predicted_solubility": round(alt_pred, 3)})
            except Exception:
                continue
        suggestions = sorted(suggestions, key=lambda x: -x["predicted_solubility"])

        response["message"] = f"{solute_name} shows poor solubility in {matched} ({pred:.3f} log mol/L)."
        response["suggested_solvents"] = suggestions[:8]
        return jsonify(response), 200

    except Exception as e:
        return jsonify({"error": "Unexpected server error", "detail": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, use_reloader=False)
