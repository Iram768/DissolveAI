import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("titanic.csv")
head = df.info()


print(head)