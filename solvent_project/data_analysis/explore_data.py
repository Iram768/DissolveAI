import pandas as pd
import matplotlib.pyplot as plt

# Load datasets
aqsol = pd.read_csv(r"../data_processed/AqSolDB_clean.csv")
bigsol = pd.read_csv(r"../data_processed/BigSolDB_clean.csv")

print("AqSolDB shape:", aqsol.shape)
print("BigSolDB shape:", bigsol.shape)
print("AqSolDB columns:", aqsol.columns.tolist())
print("BigSolDB columns:", bigsol.columns.tolist())

# Histogram of solubility
plt.hist(aqsol["Solubility"], bins=50, alpha=0.5, label="AqSolDB")
plt.hist(bigsol["Solubility"], bins=50, alpha=0.5, label="BigSolDB")
plt.xlabel("Solubility (mol/L)")
plt.ylabel("Frequency")
plt.legend()
plt.title("Distribution of Solubility Values")
plt.show()

# Count unique solvents in BigSolDB
print("Unique solvents in BigSolDB:", bigsol["Solvent"].nunique())
print("First 10 solvents:", bigsol["Solvent"].unique()[:10])

# Scatter plot: Solubility vs Temperature for BigSolDB
sample = bigsol.sample(1000)  # sample to speed up plotting
plt.scatter(sample["Temperature_K"], sample["Solubility"], alpha=0.5)
plt.xlabel("Temperature (K)")
plt.ylabel("Solubility (mol/L)")
plt.title("Solubility vs Temperature (sampled BigSolDB)")
plt.show()










