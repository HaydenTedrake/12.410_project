import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip
from scipy.interpolate import interp1d

PATHr = "r'.csv"
JD0r = 2460937.569719

dfr = pd.read_csv(PATHr)
tr = dfr["JD"].astype(float).to_numpy()
mr = dfr["Mag"].astype(float).to_numpy()
er = dfr["MagErr"].astype(float).to_numpy()

PATHg = "g'2.csv"
JD0g = 2460958.552769

dfg = pd.read_csv(PATHg)
tg = dfg["JD"].astype(float).to_numpy()
mg = dfg["Mag"].astype(float).to_numpy()
eg = dfg["MagErr"].astype(float).to_numpy()

gap_days = 0.50
P_HOURS = 3.066
title = "Light curve of 4217 Engelhardt (WAO 14-in data)"

P_days = P_HOURS / 24.0
phir = ((tr - JD0r) / P_days) % 1.0
phig = ((tg - JD0g) / P_days) % 1.0

print(phir)
print(phig)

# phi = np.linspace(0, 1, 200)
# interp1 = interp1d(phir, mr, kind='cubic', bounds_error=False, fill_value='extrapolate')
# einterp1 = interp1d(phir, er, kind='cubic', bounds_error=False, fill_value='extrapolate')
# interp2 = interp1d(phig, mg, kind='cubic', bounds_error=False, fill_value='extrapolate')
# einterp2 = interp1d(phig, eg, kind='cubic', bounds_error=False, fill_value='extrapolate')

# mr_new = interp1(phi)
# mg_new = interp2(phi)
# er_new = einterp1(phi)
# eg_new = einterp2(phi)

# e = np.sqrt(er_new**2 + eg_new**2)
# m = mg_new - mr_new

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
dates = ["20250918 UT", "20251002 UT", "20251009 UT", "20251023 UT"]

phi_fit = np.linspace(0, 1, 600)
m_fit = two_harmonic(phi_fit, *popt)
plt.plot(phi_fit, m_fit, 'k-', lw=2, label="2nd Order Fit")

mean_val = np.average(m, weights=w)
variance = np.average((m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

plt.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
plt.text(0.5, mean_val, f'Mean: {mean_val:.3f} $\\pm$ {mean_err:.3f}', 
         ha='center', va='bottom', alpha=0.9)

period_text = f'Period: {P_HOURS:.3f} hrs'
plt.text(0.02, 0.98, period_text, 
         transform=plt.gca().transAxes,
         ha='left', va='top', alpha=0.9)

plt.gca().invert_yaxis()
plt.xlim(0, 1)
plt.xlabel("Rotational Phase")
plt.ylabel("Apparent magnitude (r′)")
plt.title(title)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()