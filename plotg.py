import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "g'2.csv"
JD0 = 2460958.552769
P_HOURS = 3.066
gap_days = 0.50
title = "g' light curve of 4217 Engelhardt (WAO 14-in data)"

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

def two_harmonic(phi, a0, a1, b1, a2, b2):
    return (
        a0
        + a1 * np.cos(2 * np.pi * phi)
        + b1 * np.sin(2 * np.pi * phi)
        + a2 * np.cos(4 * np.pi * phi)
        + b2 * np.sin(4 * np.pi * phi)
    )

popt, pcov = curve_fit(
    two_harmonic,
    phi,
    m,
    sigma=e,
    absolute_sigma=True,
    p0=[np.mean(m), 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = m - two_harmonic(phi, *popt)
chi2 = np.sum((residuals / e) ** 2)
dof = len(m) - len(popt)
print(f"Reduced χ² = {chi2/dof:.2f}")

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]
dates = ["20251009 UT", "20251022 UT", "20251023 UT"]

unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    plt.errorbar(phi[sel], m[sel], yerr=e[sel],
                 fmt='.', ms=3, elinewidth=0.7, capsize=0,
                 alpha=0.95, color=colors[i % len(colors)],
                 label=f"{dates[nid]}")

phi_fit = np.linspace(0, 1, 600)
m_fit = two_harmonic(phi_fit, *popt)
np.save('g_fitted_curve.npy', m_fit)
plt.plot(phi_fit, m_fit, 'k-', lw=2, label="2nd Order Fit")

mean_val = np.average(m, weights=w)
variance = np.average((m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

plt.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
plt.text(0.5, mean_val, f'Mean: {mean_val:.3f} $\\pm$ {mean_err:.3f}', 
         ha='center', va='bottom', alpha=0.9)

period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'
plt.text(0.01, 0.02, period_text, 
         transform=plt.gca().transAxes,
         ha='left', va='bottom', alpha=0.9)

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase")
plt.ylabel("Apparent magnitude (g')")
plt.title(title)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5)
plt.tight_layout()
plt.savefig("Figures/gcurve.png")
plt.show()