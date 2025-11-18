import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "r'.csv"
JD0 = 2460937.569719
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

phi_fit = np.linspace(0, 1, 600)
m_fitr = two_harmonic(phi_fit, *popt)

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

m_fitg = two_harmonic(phi_fit, *popt)

plt.figure(figsize=(9, 5.2))
plt.plot(phi_fit, m_fitg - m_fitr, 'k-', lw=2, label="g'-r'")
plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase")
plt.ylabel("Apparent magnitude (r′)")
plt.title(title)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()