import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "g'2.csv"
JD0 = 2460958.552769
P_HOURS = 3.066
gap_days = 0.50
title = "Light curve of 4217 Engelhardt (WAO 14-in data)"

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

residuals = m - two_harmonic(phi, *popt)
mask = ~sigma_clip(residuals, sigma=2.5, maxiters=5).mask
popt, pcov = curve_fit(
    two_harmonic,
    phi[mask],
    m[mask],
    sigma=e[mask],
    absolute_sigma=True,
    p0=[np.mean(m[mask]), 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = m - two_harmonic(phi, *popt)
chi2 = np.sum((residuals / e) ** 2)
dof = len(m) - len(popt)
print(f"Reduced χ² = {chi2/dof:.2f}")

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    plt.errorbar(phi[sel], m[sel], yerr=e[sel],
                 fmt='.', ms=3, elinewidth=0.7, capsize=0,
                 alpha=0.95, color=colors[i % len(colors)],
                 label=f"Night {nid+1}")

phi_fit = np.linspace(0, 1, 600)
m_fit = two_harmonic(phi_fit, *popt)
plt.plot(phi_fit, m_fit, 'k-', lw=2, label="2nd Order Fit")

plt.axhline(np.average(m, weights=w), color='k', linestyle=':', linewidth=1.1, alpha=0.8, label="Mean")

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase")
plt.ylabel("Apparent magnitude (r′)")
plt.title(title)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()

periods = np.linspace(3.05, 3.08, 200)
chi2_vals = []

for P in periods:
    P_days = P / 24.0
    phi_test = ((t - JD0) / P_days) % 1.0
    popt_test, _ = curve_fit(two_harmonic, phi_test, m, sigma=e, absolute_sigma=True, p0=popt)
    residuals = m - two_harmonic(phi_test, *popt_test)
    chi2_vals.append(np.sum((residuals / e)**2))

best_P = periods[np.argmin(chi2_vals)]
print("Refined period:", best_P)

CSV_PATH = "r'.csv"
JD0 = 2460937.569719
P_HOURS = best_P
gap_days = 0.50
title = "Light curve of 4217 Engelhardt (WAO 14-in data)"

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

residuals = m - two_harmonic(phi, *popt)
mask = ~sigma_clip(residuals, sigma=2.5, maxiters=5).mask
popt, pcov = curve_fit(
    two_harmonic,
    phi[mask],
    m[mask],
    sigma=e[mask],
    absolute_sigma=True,
    p0=[np.mean(m[mask]), 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = m - two_harmonic(phi, *popt)
chi2 = np.sum((residuals / e) ** 2)
dof = len(m) - len(popt)
print(f"Reduced χ² = {chi2/dof:.2f}")

plt.figure(figsize=(9, 5.2))
w = 1.0 / np.maximum(e, 1e-6)**2

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    plt.errorbar(phi[sel], m[sel], yerr=e[sel],
                 fmt='.', ms=3, elinewidth=0.7, capsize=0,
                 alpha=0.95, color=colors[i % len(colors)],
                 label=f"Night {nid+1}")

phi_fit = np.linspace(0, 1, 600)
m_fit = two_harmonic(phi_fit, *popt)
plt.plot(phi_fit, m_fit, 'k-', lw=2, label="2nd Order Fit")

plt.axhline(np.average(m, weights=w), color='k', linestyle=':', linewidth=1.1, alpha=0.8, label="Mean")

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase")
plt.ylabel("Apparent magnitude (r′)")
plt.title(title)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()