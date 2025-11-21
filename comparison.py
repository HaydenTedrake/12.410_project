import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip
import pickle

PATH1 = "1010rcomp.csv"
PATH2 = "1010gcomp.csv"
PATH3 = "1024rcomp.csv"
PATH4 = "1024gcomp.csv"

df1 = pd.read_csv(PATH1)
df2 = pd.read_csv(PATH2)
df3 = pd.read_csv(PATH3)
df4 = pd.read_csv(PATH4)

t1 = df1["JD"].astype(float).to_numpy()
t2 = df2["JD"].astype(float).to_numpy()
t3 = df3["JD"].astype(float).to_numpy()
t4 = df4["JD"].astype(float).to_numpy()

m1 = df1["Mag"].astype(float).to_numpy()
m2 = df2["Mag"].astype(float).to_numpy()
m3 = df3["Mag"].astype(float).to_numpy()
m4 = df4["Mag"].astype(float).to_numpy()

e1 = df1["MagErr"].astype(float).to_numpy()
e2 = df2["MagErr"].astype(float).to_numpy()
e3 = df3["MagErr"].astype(float).to_numpy()
e4 = df4["MagErr"].astype(float).to_numpy()

t = []
m = []
e = []
for i in range(len(t2)):
    t.append((t1[i]+t2[i])/2)
    m.append(m1[i]-m2[i])
    e.append(np.sqrt(e1[i]**2+e2[i]**2))
for i in range(len(t3)):
    t.append((t3[i]+t4[i])/2)
    m.append(m3[i]-m4[i])
    e.append(np.sqrt(e3[i]**2+e4[i]**2))
t = np.array(t)
m = np.array(m)
e = np.array(e)
print(t)

JD0 = 2460958.5558045
P_HOURS = 3.066
P_days = P_HOURS / 24.0
gap_days = 0.50
phi = ((t - JD0) / P_days) % 1.0
title = "g'-r' color index of comparison stars"

dt = np.diff(t, prepend=t[0])
night_id = np.zeros_like(t, dtype=int)
for i in range(1, len(t)):
    night_id[i] = night_id[i-1] + (dt[i] > gap_days)

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = [
    "#009E73",  # bluish green
    "#E69F00",  # orange
]
markers = ["o", "s", "^", "D"]
dates = ["20251010 UT", "20251024 UT"]

unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    w_night = w[sel]
    m_night = m[sel]
    
    mean_night = np.average(m_night, weights=w_night)
    var_night = np.average((m_night - mean_night)**2, weights=w_night)
    weight_sum_night = np.sum(w_night)
    mean_err_night = np.sqrt(var_night / weight_sum_night)
    
    plt.axhline(mean_night, color='k', linestyle='-', linewidth=1.1, alpha=0.8)
    if nid == 0:
        plt.text(0.5, mean_night + 0.03, f'Mean: {mean_night:.2f} $\\pm$ {mean_err_night:.2f}', 
             ha='center', va='bottom', alpha=0.9, fontsize=12)
    else:
        plt.text(0.5, mean_night, f'Mean: {mean_night:.2f} $\\pm$ {mean_err_night:.2f}', 
             ha='center', va='bottom', alpha=0.9, fontsize=12)
    
    plt.errorbar(
        phi[sel], m[sel], yerr=e[sel],
        fmt=markers[i % len(markers)],
        markersize=4.5,
        elinewidth=0.7,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
        label=f"Comp star {nid+1} ({dates[nid]})",
    )

period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'
plt.text(0.01, 0.02, period_text, 
         transform=plt.gca().transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase of 4217 Engelhardt", fontsize=12)
plt.ylabel("Apparent magnitude (g'-r')", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12) 
plt.title(title, fontsize=12)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)
plt.tight_layout()
plt.savefig("Figures/colorcurve.png")
plt.show()