import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "r'.csv"

df = pd.read_csv(CSV_PATH)
e = df["MagErr"].astype(float).to_numpy()
# print(e)
print("first")
print(np.min(e))
print(np.mean(e))
print(np.median(e))

# notes
# for 1010, use r19 and g18, for 1024 use 18 and newg18
