# predict_solubility_advanced.py
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load trained model and encoder
model = joblib.load("models/solubility_model1.pkl")
encoder = joblib.load("models/solvent_encoder1.pkl")

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

def predict_solubility(model, encoder, solute_smiles, solvent_name, temperature):
    # Compute solute descriptors
    mw, logp, tpsa = compute_descriptors(solute_smiles)

    # Encode solvent
    if solvent_name not in encoder.categories_[0]:
        print(f"⚠️ Warning: {solvent_name} was not seen during training. Prediction may be unreliable.")
        solvent_encoded = [0]*len(encoder.categories_[0])
    else:
        solvent_encoded = encoder.transform(pd.DataFrame({'Solvent':[solvent_name]})).tolist()[0]

    # Build input features
    feature_names = ["Temperature_K","MW","LogP","TPSA"] + encoder.get_feature_names_out(['Solvent']).tolist()
    X_input = pd.DataFrame([[temperature, mw, logp, tpsa] + solvent_encoded], columns=feature_names)

    # Predict solubility
    return model.predict(X_input)[0]

if __name__ == "__main__":
    print("Enter solute and solvent details:")
    solute_smiles = input("Solute SMILES: ")
    solvent_name = input("Solvent name: ")
    temperature = float(input("Temperature (K): "))

    solubility_threshold = -2  # log mol/L

    try:
        pred = predict_solubility(model, encoder, solute_smiles, solvent_name, temperature)
        print(f"Predicted solubility (log mol/L): {pred:.4f}")

        if pred < solubility_threshold:
            print(f"⚠️ This solute is poorly soluble in {solvent_name}.")

            # Suggest alternative solvents
            known_solvents = encoder.categories_[0]
            alt_solubility = {}
            for s in known_solvents:
                alt_solubility[s] = predict_solubility(model, encoder, solute_smiles, s, temperature)

            # Sort by solubility descending
            suggestions = sorted(alt_solubility.items(), key=lambda x: x[1], reverse=True)
            print("💡 You can try these solvents instead (higher solubility):")
            for s, val in suggestions[:3]:
                if s != solvent_name:
                    print(f"  {s}: {val:.4f}")

    except ValueError as e:
        print(f"Error: {e}")
