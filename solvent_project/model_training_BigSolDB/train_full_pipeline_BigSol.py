# train_full_pipeline_bigsol.py
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
# 1. Load BigSolDB
# ----------------------------
data_path = "data_processed/BigSolDB_clean.csv"
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
df = df.dropna(subset=['MW','LogP','TPSA'])

# ----------------------------
# 3. One-hot encode solvents
# ----------------------------
encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore')
encoder.fit(df[['Solvent']])   # now includes many solvents
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
    n_estimators=300,
    learning_rate=0.05,
    max_depth=6,
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
joblib.dump(model, "models/solubility_model_bigsol.pkl")
joblib.dump(encoder, "models/solvent_encoder_bigsol.pkl")

print("✅ Model and encoder saved successfully!")
print("Unique solvents in dataset:", df['Solvent'].unique()[:10], "...")
print("Total solvents:", len(encoder.categories_[0]))
