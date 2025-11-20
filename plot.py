import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "r'.csv"
JD0 = 2460937.567619
P_HOURS = 3.066
gap_days = 0.50
title = "r' light curve of 4217 Engelhardt (WAO 14-in data)"

df = pd.read_csv(CSV_PATH)
t = df["JD"].astype(float).to_numpy()
m = df["Mag"].astype(float).to_numpy()
e = df["MagErr"].astype(float).to_numpy()

P_days = P_HOURS / 24.0
phi = ((t - JD0) / P_days) % 1.0

dt = np.diff(t, prepend=t[0])
night_id = np.zeros_like(t, dtype=int)
for i in range(1, len(t)):
    night_id[i] = night_id[i-1] + (dt[i] > gap_days)

def three_harmonic(phi, a0, a1, b1, a2, b2, a3, b3):
    return (
        a0
        + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
        + b1 * np.sin(2 * np.pi * phi)
        + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
        + b2 * np.sin(4 * np.pi * phi)
        + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
        + b3 * np.sin(6 * np.pi * phi)
    )

popt, pcov = curve_fit(
    three_harmonic,
    phi,
    m,
    sigma=e,
    absolute_sigma=True,
    p0=[np.mean(m), 0, 0, 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = m - three_harmonic(phi, *popt)
chi2 = np.sum((residuals / e) ** 2)
dof = len(m) - len(popt)
print(f"Reduced χ² = {chi2/dof:.2f}")

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#CC79A7",  # purple/magenta
]
markers = ["o", "s", "^", "D"]
dates = ["20250919 UT", "20251003 UT", "20251010 UT", "20251024 UT"]

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

tenten = []
tentwentyfour = []
for nid in unique_nights:
    sel = (night_id == nid)
    if nid == 2:
        tenten.extend(zip(t[sel], m[sel], e[sel]))
    elif nid == 3:
        tentwentyfour.extend(zip(t[sel], m[sel], e[sel]))
with open("rtenten.pkl", "wb") as f:
    pickle.dump(tenten, f)
with open("rtentwentyfour.pkl", "wb") as f:
    pickle.dump(tentwentyfour, f)

phi_fit = np.linspace(0, 1, 600)
m_fit = three_harmonic(phi_fit, *popt)
np.save('r_fitted_curve.npy', m_fit)
plt.plot(phi_fit, m_fit, 'k-', lw=2, label="3rd Order Fit")

mean_val = np.average(m, weights=w)
variance = np.average((m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

plt.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
plt.text(0.515, mean_val, f'Mean: {mean_val:.2f} $\\pm$ {mean_err:.2f}', 
         ha='center', va='bottom', alpha=0.9, fontsize=12)

period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'
plt.text(0.01, 0.02, period_text, 
         transform=plt.gca().transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase", fontsize=12)
plt.ylabel("Apparent magnitude (r')", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12) 
plt.title(title, fontsize=12)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper left", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)
plt.tight_layout()
plt.savefig("Figures/rcurve.png")
plt.show()