import os
import json
import joblib
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FEATURES_PATH = os.path.join(BASE_DIR, "features", "BigSolDB_features_numeric.csv")
CLEAN_DATA_PATH = os.path.join(BASE_DIR, "data_processed", "BigSolDB_clean.csv")
MODEL_OUT_PATH = os.path.join(BASE_DIR, "models", "solubility_model_bigsol.pkl")
ENCODER_OUT_PATH = os.path.join(BASE_DIR, "models", "solvent_encoder_bigsol.pkl")
METADATA_OUT_PATH = os.path.join(BASE_DIR, "models", "metadata.json")

def main():
    print("Loading precomputed features...")
    df_features = pd.read_csv(FEATURES_PATH)
    df_clean = pd.read_csv(CLEAN_DATA_PATH)

    print("Fitting OneHotEncoder on all solvents...")
    encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    encoder.fit(df_clean[["Solvent"]])

    # Verify matching columns
    solvent_cols = encoder.get_feature_names_out(["Solvent"]).tolist()
    feature_cols = ["Temperature_K", "MW", "LogP", "TPSA"] + solvent_cols

    X = df_features[feature_cols]
    y = df_features["Solubility"]

    print(f"X shape: {X.shape}, y shape: {y.shape}")

    # Train/Test Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    print("Training Random Forest Regressor (500 trees)...")
    model = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=42, max_depth=20) # max_depth=20 to keep size reasonable and fast training
    model.fit(X_train, y_train)

    # Evaluate
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    print(f"Test RMSE: {rmse:.4f}")

    # Save
    print("Saving model, encoder, and metadata...")
    os.makedirs(os.path.dirname(MODEL_OUT_PATH), exist_ok=True)
    joblib.dump(model, MODEL_OUT_PATH)
    joblib.dump(encoder, ENCODER_OUT_PATH)

    # Build metadata
    temp_min = float(df_clean["Temperature_K"].min())
    temp_max = float(df_clean["Temperature_K"].max())
    sol_min = float(df_features["Solubility"].min())
    sol_max = float(df_features["Solubility"].max())

    metadata = {
        "temperature_min": temp_min,
        "temperature_max": temp_max,
        "solubility_lower_bound": sol_min,
        "solubility_upper_bound": sol_max,
        "solubility_threshold": -2.0,
        "total_solvents": len(solvent_cols),
        "unique_solvents_example": df_clean["Solvent"].unique().tolist()[:10]
    }

    with open(METADATA_OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print("Training and saving completed successfully!")

if __name__ == "__main__":
    main()
