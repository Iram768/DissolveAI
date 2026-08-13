# train_full_pipeline.py
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
import joblib

# ----------------------------
# 1. Load cleaned dataset
# ----------------------------
data_path = "data_processed/AqSolDB_clean.csv"
df = pd.read_csv(data_path)

# ----------------------------
# 2. Compute solute descriptors
# ----------------------------
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.TPSA(mol)

df[['MW','LogP','TPSA']] = df['Solute'].apply(lambda s: pd.Series(compute_descriptors(s)))

# Drop rows with invalid descriptors
df = df.dropna(subset=['MW','LogP','TPSA'])

# ----------------------------
# 3. One-hot encode solvents (all solvents in dataset)
# ----------------------------
encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
encoder.fit(df[['Solvent']])  # important: fit on all known solvents
solvent_encoded = encoder.transform(df[['Solvent']])
solvent_cols = encoder.get_feature_names_out(['Solvent'])
solvent_df = pd.DataFrame(solvent_encoded, columns=solvent_cols, index=df.index)

# ----------------------------
# 4. Combine features and target
# ----------------------------
X = pd.concat([df[['Temperature_K','MW','LogP','TPSA']], solvent_df], axis=1)
y = df['Solubility']

# ----------------------------
# 5. Train-test split
# ----------------------------
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# ----------------------------
# 6. Train Gradient Boosting model
# ----------------------------
model = GradientBoostingRegressor(
    n_estimators=200,
    learning_rate=0.1,
    max_depth=5,
    random_state=42
)
model.fit(X_train, y_train)

# ----------------------------
# 7. Evaluate
# ----------------------------
y_pred = model.predict(X_test)
rmse = mean_squared_error(y_test, y_pred) ** 0.5
print(f"Root Mean Squared Error: {rmse:.4f}")

# ----------------------------
# 8. Save model and encoder
# ----------------------------
os.makedirs("models", exist_ok=True)
joblib.dump(model, "models/solubility_model1.pkl")
joblib.dump(encoder, "models/solvent_encoder1.pkl")
print("✅ Model and encoder saved successfully!")



print("Unique solvents in dataset:", df['Solvent'].unique())
print("Encoder categories:", encoder.categories_)
