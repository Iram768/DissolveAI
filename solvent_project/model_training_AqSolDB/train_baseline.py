########################################################## slightly small train model 

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# 1. Load the feature dataset
data = pd.read_csv("features/AqSolDB_features_numeric.csv")  # path to your generated features

# 2. Define inputs (X) and target (y)
X = data.drop(columns=["Solubility"])  # all columns except target
y = data["Solubility"]                 # target column

# 3. Split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 4. Initialize and train the model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# 5. Predict on test set
y_pred = model.predict(X_test)

# 6. Calculate RMSE (future-proof)
rmse = mean_squared_error(y_test, y_pred) ** 0.5

print(f"Root Mean Squared Error: {rmse:.4f}")



##################################### To check columns in features files



# import pandas as pd

# data = pd.read_csv("features/AqSolDB_features_numeric.csv")
# print(data.columns.tolist())





