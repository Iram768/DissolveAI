import pandas as pd

df = pd.read_csv("solvents_data.csv")

#print(df.info())
#print(df.describe())
print(df["Polarity"].value_counts())
