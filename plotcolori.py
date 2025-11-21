import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip
import pickle

with open("rtenten.pkl", "rb") as f:
    rtenten = pickle.load(f)
with open("rtentwentyfour.pkl", "rb") as f:
    rtentwentyfour = pickle.load(f)
with open("gtenten.pkl", "rb") as f:
    gtenten = pickle.load(f)
with open("gtentwentyfour.pkl", "rb") as f:
    gtentwentyfour = pickle.load(f)
    gtentwentyfour = gtentwentyfour[:13]

t = []
m = []
e = []
for i in range(len(rtenten)):
    t.append((rtenten[i][0]+gtenten[i][0])/2)
    m.append(gtenten[i][1]-rtenten[i][1])
    e.append(np.sqrt(rtenten[i][2]**2+gtenten[i][2]**2))
for i in range(len(rtentwentyfour)):
    t.append((rtentwentyfour[i][0]+gtentwentyfour[i][0])/2)
    m.append(gtentwentyfour[i][1]-rtentwentyfour[i][1])
    e.append(np.sqrt(rtentwentyfour[i][2]**2+gtentwentyfour[i][2]**2))
t = np.array(t)
m = np.array(m)
e = np.array(e)

P_HOURS = 3.066
P_days = P_HOURS / 24.0
gap_days = 0.50
JD0 = 2460958.55063
phi = ((t - JD0) / P_days) % 1.0
title = "g'-r' color index of 4217 Engelhardt (WAO 14-in data)"

dt = np.diff(t, prepend=t[0])
night_id = np.zeros_like(t, dtype=int)
for i in range(1, len(t)):
    night_id[i] = night_id[i-1] + (dt[i] > gap_days)

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = [
    "#009E73",  # bluish green
    "#CC79A7",  # purple/magenta
]
markers = ["o", "s", "^", "D"]
dates = ["20251010 UT", "20251024 UT"]

unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    plt.errorbar(
        phi[sel], m[sel], yerr=e[sel],
        fmt=markers[i % len(markers)],
        markersize=4.5,
        elinewidth=0.7,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
        label=f"{dates[nid]}",
    )

mean_val = np.average(m, weights=w)
variance = np.average((m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

plt.axhline(mean_val, color='k', linestyle='-', linewidth=1.1, alpha=0.8)
plt.text(0.5, mean_val + 0.015, f'Mean: {mean_val:.2f} $\\pm$ {mean_err:.2f}', 
         ha='center', va='bottom', alpha=0.9, fontsize=12)

period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'
plt.text(0.01, 0.02, period_text, 
         transform=plt.gca().transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase", fontsize=12)
plt.ylabel("Apparent magnitude (g'-r')", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12) 
plt.title(title, fontsize=12)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper left", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)
plt.tight_layout()
plt.savefig("Figures/colorcurve.png")
plt.show()