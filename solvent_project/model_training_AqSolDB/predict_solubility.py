# # predict_solubility.py
# import joblib
# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem import Descriptors

# # Load trained model and encoder
# model = joblib.load("models/solubility_model.pkl")
# encoder = joblib.load("models/solvent_encoder.pkl")

# def compute_descriptors(smiles):
#     """Compute MW, LogP, TPSA for a given SMILES string"""
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         raise ValueError("Invalid SMILES string")
#     mw = Descriptors.MolWt(mol)
#     logp = Descriptors.MolLogP(mol)
#     tpsa = Descriptors.TPSA(mol)
#     return mw, logp, tpsa

# def predict_solubility(model, encoder, solute_smiles, solvent_name, temperature):
#     # Compute solute descriptors
#     mw, logp, tpsa = compute_descriptors(solute_smiles)

#     # Encode solvent
#     solvent_encoded = encoder.transform(pd.DataFrame({'Solvent': [solvent_name]}))

#     # Build input feature row
#     feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(["Solvent"]).tolist()
#     X_input = pd.DataFrame([[temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]], columns=feature_names)

#     # Predict solubility
#     return model.predict(X_input)[0]

# if __name__ == "__main__":
#     # User input
#     solute_smiles = input("Enter solute SMILES: ")
#     solvent_name = input("Enter solvent name: ")
#     temperature = float(input("Enter temperature (K): "))

#     # Prediction
#     pred = predict_solubility(model, encoder, solute_smiles, solvent_name, temperature)
#     print(f"Predicted solubility (log mol/L): {pred:.4f}")


#########################################################################################################################


# predict_solubility_complex.py
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load trained model and encoder
model = joblib.load("models/solubility_model.pkl")
encoder = joblib.load("models/solvent_encoder.pkl")

def compute_descriptors(smiles):
    """Compute MW, LogP, TPSA for a given SMILES string"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    return mw, logp, tpsa

def predict_solubility(model, encoder, solute_smiles, solvent_name, temperature):
    # Compute solute descriptors
    mw, logp, tpsa = compute_descriptors(solute_smiles)

    # Encode solvent
    # For difficult solvents, make sure they were included in encoder.fit() earlier
    solvent_encoded = encoder.transform(pd.DataFrame({'Solvent': [solvent_name]}))

    # Build input feature row
    feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(["Solvent"]).tolist()
    X_input = pd.DataFrame([[temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]], columns=feature_names)

    # Predict solubility
    return model.predict(X_input)[0]

if __name__ == "__main__":
    print("Enter solute and organic solvent details:")
    solute_smiles = input("Solute SMILES: ")
    solvent_name = input("Solvent name (e.g., dimethylformamide, DMSO, toluene): ")
    temperature = float(input("Temperature (K): "))

    try:
        pred = predict_solubility(model, encoder, solute_smiles, solvent_name, temperature)
        print(f"Predicted solubility (log mol/L): {pred:.4f}")
    except ValueError as e:
        print(f"Error: {e}. Make sure the solvent is in the trained encoder list.")


solubility_threshold = -2  # log mol/L; adjust based on your dataset
if pred < solubility_threshold:
    print(f"⚠️ This solute is poorly soluble in {solvent_name}.")





###make its user interface
## in which the it depends on user to add temp 

## if user is not adding temp that means that it have to predict solubility at specific temp but if user is not adding temp so it should print the result like this solute (solute name extracted from smiles) is slouble in solvet at room temp 298K 