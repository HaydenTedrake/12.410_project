import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "g17.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()
# print(e)
print("first")
print(np.median(e))
print(np.mean(e))

CSV_PATH = "newg18.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()

print("second")
print(np.median(e))
print(np.mean(e))

CSV_PATH = "g19.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()
# print(e)
print("last")
print(np.median(e))
print(np.mean(e))

# notes
# for 1010, use r19 and g18, for 1024 use 18 and newg18