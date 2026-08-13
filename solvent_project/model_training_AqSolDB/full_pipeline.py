# full_pipeline.py
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
import joblib  # to save trained model

# ----------------------------
# 1. Load cleaned dataset
# ----------------------------
data_path = "data_processed/AqSolDB_clean.csv"  # change if needed
df = pd.read_csv(data_path)

# ----------------------------
# 2. Generate solute descriptors
# ----------------------------
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    return mw, logp, tpsa

df[['MW', 'LogP', 'TPSA']] = df['Solute'].apply(
    lambda s: pd.Series(compute_descriptors(s))
)

# Drop any rows where descriptors could not be computed
df = df.dropna(subset=['MW', 'LogP', 'TPSA'])

# ----------------------------
# 3. One-hot encode solvents
# ----------------------------
encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
solvent_encoded = encoder.fit_transform(df[['Solvent']])
solvent_cols = encoder.get_feature_names_out(['Solvent'])
solvent_df = pd.DataFrame(solvent_encoded, columns=solvent_cols, index=df.index)

# Combine features
df_features = pd.concat([df[['Temperature_K', 'MW', 'LogP', 'TPSA']], solvent_df], axis=1)
y = df['Solubility']

# ----------------------------
# 4. Train-test split
# ----------------------------
X_train, X_test, y_train, y_test = train_test_split(df_features, y, test_size=0.2, random_state=42)

# ----------------------------
# 5. Train Gradient Boosting model
# ----------------------------
model = GradientBoostingRegressor(
    n_estimators=200,
    learning_rate=0.1,
    max_depth=5,
    random_state=42
)
model.fit(X_train, y_train)

# ----------------------------
# 6. Predict & evaluate
# ----------------------------
y_pred = model.predict(X_test)
rmse = mean_squared_error(y_test, y_pred) ** 0.5
print(f"Root Mean Squared Error: {rmse:.4f}")

# ----------------------------
# 7. Save model & encoder
# ----------------------------
os.makedirs("models", exist_ok=True)
joblib.dump(model, "models/solubility_model.pkl")
joblib.dump(encoder, "models/solvent_encoder.pkl")

# ----------------------------
# 8. Optional: Save feature dataset
# ----------------------------
os.makedirs("features", exist_ok=True)
df_features.to_csv("features/AqSolDB_features_numeric.csv", index=False)
print("Pipeline completed: model and features saved.")



# import joblib
# from sklearn.preprocessing import OneHotEncoder

# # Fit encoder on the 'Solvent' column
# encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
# encoder.fit(X_train[['Solvent']])

# # Save both encoder and model
# joblib.dump(encoder, "features/solvent_encoder.pkl")
# joblib.dump(model, "features/solubility_model.pkl")

# print("✅ Model and encoder saved!")
