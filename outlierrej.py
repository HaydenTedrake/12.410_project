import numpy as np
import pandas as pd
from scipy.stats import norm


r = "r'.csv"
g = "g'.csv"

dfr = pd.read_csv(r)
dfg = pd.read_csv(g)

tr = dfr["JD"].astype(float).to_numpy()
mr = dfr["Mag"].astype(float).to_numpy()
er = dfr["MagErr"].astype(float).to_numpy()

tg = dfg["JD"].astype(float).to_numpy()
mg = dfg["Mag"].astype(float).to_numpy()
eg = dfg["MagErr"].astype(float).to_numpy()

meanr = np.mean(mr)
print(meanr)
meang = np.mean(mg)
print(meang)

stdr = np.std(mr)
print(stdr)
stdg = np.std(mg)
print(stdg)

print(len(mr))
print(len(mg))
print(mr)
print(mg)

print("hi")

for m in mr:
    t = abs(m - meanr) / stdr
    p = 2 - 2 * norm.cdf(t)
    l = p * 91
    if l < 0.5:
        print(m)

print("next")
for m in mg:
    t = abs(m - meang) / stdg
    p = 2 - 2 * norm.cdf(t)
    l = p * 56
    if l < 0.5:
        print(m)