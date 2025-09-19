# # app.py
# from flask import Flask, request, jsonify
# import joblib
# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem import Descriptors
# import requests

# import os

# import joblib

# app = Flask(__name__)

# # --- Visitor Counter ---
# VISIT_FILE = "visit_count.txt"

# def get_visit_count():
#     if not os.path.exists(VISIT_FILE):
#         with open(VISIT_FILE, "w") as f:
#             f.write("0")
#         return 0
#     with open(VISIT_FILE, "r") as f:
#         return int(f.read().strip())

# def increment_visit_count():
#     count = get_visit_count() + 1
#     with open(VISIT_FILE, "w") as f:
#         f.write(str(count))
#     return count

# @app.route("/visit-count", methods=["GET"])
# def visit_count():
#     count = increment_visit_count()
#     return jsonify({"visits": count})

# # --- your existing solubility prediction routes continue below ---

# # Load model + encoder
# model = joblib.load("d:/AI_projects/DissolveAI/solvent_project/models/solubility_model_bigsol.pkl")
# encoder = joblib.load("d:/AI_projects/DissolveAI/solvent_project/models/solvent_encoder_bigsol.pkl")

# app = Flask(__name__)

# def compute_descriptors(smiles):
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         return None
#     return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

# def get_compound_name(smiles):
#     """Query PubChem API to get IUPAC or preferred name from SMILES"""
#     try:
#         url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON"
#         r = requests.get(url, timeout=5)
#         if r.status_code == 200:
#             data = r.json()
#             return data["PropertyTable"]["Properties"][0]["IUPACName"]
#         else:
#             return smiles  # fallback: return SMILES if not found
#     except Exception:
#         return smiles  # fallback if API fails

# @app.route("/predict", methods=["POST"])
# def predict():
#     data = request.json
#     solute_smiles = data.get("solute_smiles")
#     solvent = data.get("solvent")
#     temperature = data.get("temperature")  # may be None

#     # compute descriptors
#     desc = compute_descriptors(solute_smiles)
#     if desc is None:
#         return jsonify({"error": "Invalid SMILES"}), 400
#     mw, logp, tpsa = desc

#     # encode solvent
#     solvent_encoded = encoder.transform(pd.DataFrame({'Solvent': [solvent]}))

#     # set default temperature if missing
#     if temperature is None:
#         temperature = 298  # room temp

#     # build features
#     feature_names = ["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(['Solvent']).tolist()
#     X_input = pd.DataFrame([[temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]], columns=feature_names)

#     # predict
#     solubility = model.predict(X_input)[0]

#     # get compound name from PubChem
#     solute_name = get_compound_name(solute_smiles)

#     # message
#     if data.get("temperature") is None:
#         message = f"{solute_name} is predicted to be soluble in {solvent} at room temp (298K): {solubility:.4f} log mol/L"
#     else:
#         message = f"Predicted solubility of {solute_name} in {solvent} at {temperature}K: {solubility:.4f} log mol/L"

#     return jsonify({"prediction": solubility, "message": message})

# if __name__ == "__main__":
#     app.run(debug=True)










# # Load model + encoder
# model = joblib.load("models/solubility_model1.pkl")
# encoder = joblib.load("models/solvent_encoder1.pkl")

# # --- Helper: extract solute name from SMILES ---
# def smiles_to_name(smiles):
#     mol = Chem.MolFromSmiles(smiles)
#     if mol:
#         return Chem.MolToSmiles(mol)  # canonical SMILES (acts like a name)
#     return "Unknown Solute"

# @app.route("/predict", methods=["POST"])
# def predict():
#     data = request.get_json()
#     solute_smiles = data.get("solute")
#     solvent = data.get("solvent")
#     temperature = data.get("temperature")

#     # ✅ Compute descriptors
#     mol = Chem.MolFromSmiles(solute_smiles)
#     if mol is None:
#         return jsonify({"error": "Invalid SMILES"}), 400

#     MW = Descriptors.MolWt(mol)
#     LogP = Descriptors.MolLogP(mol)
#     TPSA = Descriptors.TPSA(mol)

#     # ✅ Handle temperature
#     if not temperature or str(temperature).strip() == "":
#         temperature = 298  # Default room temperature
#         user_provided_temp = False
#     else:
#         temperature = float(temperature)
#         user_provided_temp = True

#     # ✅ Encode solvent
#     import pandas as pd
#     solvent_encoded = encoder.transform([[solvent]])
#     solvent_cols = encoder.get_feature_names_out(["Solvent"])
#     solvent_df = pd.DataFrame(solvent_encoded, columns=solvent_cols)

#     # ✅ Create feature row
#     X_new = pd.DataFrame(
#         [[temperature, MW, LogP, TPSA] + list(solvent_encoded[0])],
#         columns=["Temperature_K", "MW", "LogP", "TPSA"] + list(solvent_cols)
#     )

#     # ✅ Predict
#     predicted_solubility = model.predict(X_new)[0]

#     # ✅ Extract solute name
#     solute_name = smiles_to_name(solute_smiles)

#     return jsonify({
#         "solute_name": solute_name,
#         "solvent": solvent,
#         "temperature": temperature,
#         "predicted_solubility": float(predicted_solubility),
#         "user_provided_temp": user_provided_temp
#     })

# if __name__ == "__main__":
#     app.run(debug=True)


########################-------------------------------------------------------
##### The above code is with visit count ..............


###### This is final code
from flask import Flask, request, jsonify
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from flask_cors import CORS

app = Flask(__name__)
CORS(app)  # allow frontend to connect

# Load trained model and encoder
model = joblib.load("d:/AI_projects/DissolveAI/solvent_project/models/solubility_model_bigsol.pkl")
encoder = joblib.load("d:/AI_projects/DissolveAI/solvent_project/models/solvent_encoder_bigsol.pkl")

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    solute_smiles = data.get("solute_smiles")
    solvent_name = data.get("solvent_name")
    temperature = data.get("temperature", 298)  # default room temp

    mw, logp, tpsa = compute_descriptors(solute_smiles)
    solvent_encoded = encoder.transform(pd.DataFrame({'Solvent': [solvent_name]}))

    features = pd.DataFrame(
        [[temperature, mw, logp, tpsa] + solvent_encoded.tolist()[0]],
        columns=["Temperature_K", "MW", "LogP", "TPSA"] + encoder.get_feature_names_out(["Solvent"]).tolist()
    )

    pred = model.predict(features)[0]
    return jsonify({"prediction": f"{pred:.4f} log mol/L"})

if __name__ == "__main__":
    app.run(debug=True)

