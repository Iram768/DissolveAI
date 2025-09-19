import pandas as pd

# Load the main solubility dataset
bigsol = pd.read_csv("BigSolDBv2.0.csv")
# print("BigSolDB shape:", bigsol.shape)
# print(bigsol.head())



ethanol_data = bigsol[bigsol["Solvent"] == "ethanol"]
#print("Ethanol dataset shape:", ethanol_data.shape)
#print(ethanol_data.head())


print(bigsol.columns)   # See all column names
print(bigsol.info())    # See data types + missing values
print(bigsol.describe())  # Quick stats (for numeric columns)





