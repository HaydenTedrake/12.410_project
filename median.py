import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "15.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()
# print(e)
print("first")
print(np.median(e))
print(np.mean(e))

CSV_PATH = "13.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()

print("second")
print(np.median(e))
print(np.mean(e))

CSV_PATH = "14.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()
# print(e)
print("last")
print(np.median(e))
print(np.mean(e))