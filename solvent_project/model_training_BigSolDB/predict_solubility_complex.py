# predict_solubility_complex.py
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# ----------------------------
# Load trained model + encoder
# ----------------------------
model = joblib.load("models/solubility_model_bigsol.pkl")
encoder = joblib.load("models/solvent_encoder_bigsol.pkl")

# ----------------------------
# Helper: compute solute descriptors
# ----------------------------
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("❌ Invalid SMILES string")
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

# ----------------------------
# Prediction function
# ----------------------------
def predict_solubility(model, encoder, solute_smiles, solvent_name, temperature):
    # Compute descriptors
    mw, logp, tpsa = compute_descriptors(solute_smiles)

    # Encode solvent
    solvent_encoded = encoder.transform(pd.DataFrame({'Solvent': [solvent_name]}))

    # Build feature row
    feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(['Solvent']).tolist()
    X_input = pd.DataFrame([[temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]], columns=feature_names)

    # Predict
    return model.predict(X_input)[0]

# ----------------------------
# Main script
# ----------------------------
if __name__ == "__main__":
    print("🔮 Solubility Prediction Tool (BigSolDB-trained)")
    solute_smiles = input("Enter solute SMILES: ")
    solvent_name = input("Enter solvent name (e.g., water, DMSO, ethanol, toluene): ")
    temperature = float(input("Enter temperature (K): "))

    try:
        pred = predict_solubility(model, encoder, solute_smiles, solvent_name, temperature)
        print(f"\nPredicted solubility in {solvent_name} at {temperature} K: {pred:.4f} log mol/L")

        # Add solubility threshold warning
        solubility_threshold = -2.0  # you can tune this
        if pred < solubility_threshold:
            print(f"⚠️ Warning: Poor solubility predicted in {solvent_name}. Consider alternative solvents.")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
