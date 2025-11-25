import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "g'.csv"
JD0 = 2460958.552769
P_HOURS = 3.0661
gap_days = 0.50
title = "g' light curve of 4217 Engelhardt (WAO 14-in data)"

df = pd.read_csv(CSV_PATH)
t = df["JD"].astype(float).to_numpy()
m = df["Mag"].astype(float).to_numpy()
e = df["MagErr"].astype(float).to_numpy()

P_days = P_HOURS / 24.0
phi = ((t - JD0) / P_days) % 1.0

new_m = []
dt = np.diff(t, prepend=t[0])
night_id = np.zeros_like(t, dtype=int)
for i in range(0, len(t)):
    night_id[i] = night_id[i-1] + (dt[i] > gap_days)
    if night_id[i] == 1:
        new_m.append(m[i]+0.05179233333)
    else:
        new_m.append(m[i])
new_m = np.array(new_m)

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
    new_m,
    sigma=e,
    absolute_sigma=True,
    p0=[np.mean(new_m), 0, 0, 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = new_m - three_harmonic(phi, *popt)
chi2 = np.sum((residuals / e) ** 2)
dof = len(new_m) - len(popt)
print(f"Reduced χ² = {chi2/dof:.2f}")

# Create single plot
fig, ax1 = plt.subplots(1, 1, figsize=(9, 6))

w = 1.0 / np.maximum(e, 1e-6)**2

colors = [
    "#009E73",  # bluish green
    "#CC79A7",  # purple/magenta
    "#E69F00",  # orange
]
markers = ["o", "s", "^", "D"]
dates = ["20251010 UT", "20251023 UT", "20251024 UT"]

# Plot light curve on main plot
unique_nights = np.unique(night_id)
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    ax1.errorbar(
        phi[sel], new_m[sel], yerr=e[sel],
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
    if nid == 0:
        tenten.extend(zip(t[sel], new_m[sel], e[sel]))
    elif nid == 2:
        tentwentyfour.extend(zip(t[sel], new_m[sel], e[sel]))
with open("gtenten.pkl", "wb") as f:
    pickle.dump(tenten, f)
with open("gtentwentyfour.pkl", "wb") as f:
    pickle.dump(tentwentyfour, f)

phi_fit = np.linspace(0, 1, 600)
new_m_fit = three_harmonic(phi_fit, *popt)
ax1.plot(phi_fit, new_m_fit, 'k-', lw=2, label="3rd Order Fit")

mean_val = np.average(new_m, weights=w)
variance = np.average((new_m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

ax1.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
ax1.text(0.775, mean_val, f'Mean: {mean_val:.2f} $\\pm$ {mean_err:.2f}', 
         ha='center', va='bottom', alpha=0.9, fontsize=12)

period_text = f'Period: {P_HOURS} $\\pm$ 0.0002 hrs'
ax1.text(0.01, 0.02, period_text, 
         transform=ax1.transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

ax1.invert_yaxis()
ax1.set_xlim(0, 1)
ax1.tick_params(axis='both', which='major', labelsize=12) 
ax1.set_ylabel("Apparent magnitude (g')", fontsize=12)
ax1.set_xlabel("Rotational Phase", fontsize=12)
ax1.set_title(title, fontsize=12)
ax1.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
ax1.legend(loc="upper left", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)

# Create inset axes for residuals with exact position control
# [left, bottom, width, height] in figure coordinates (0-1)
left = 0.47    # Distance from left edge
bottom = 0.155  # Distance from bottom edge  
width = 0.19   # Width of inset
height = 0.19  # Height of inset

ax2 = fig.add_axes([left, bottom, width, height])

# Plot residuals on inset
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    ax2.errorbar(
        phi[sel], residuals[sel], yerr=e[sel],
        fmt=markers[i % len(markers)],
        markersize=3.0,
        elinewidth=0.5,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
    )

# Calculate mean and error of residuals
residual_mean = np.average(residuals, weights=w)
residual_variance = np.average((residuals - residual_mean)**2, weights=w)
residual_mean_err = np.sqrt(residual_variance / weight_sum)

ax2.axhline(0, color='k', linestyle='-', linewidth=0.8, alpha=0.8)
ax2.axhline(residual_mean, color='k', linestyle=':', linewidth=0.8, alpha=0.8)

ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_title("Residuals", fontsize=12)  # Same font size as main plot

# Set consistent font sizes for all ticks
ax2.tick_params(axis='both', which='major', labelsize=12)  # Same as main plot

ax2.grid(True, linestyle=':', linewidth=0.5, alpha=0.7)

plt.tight_layout()
plt.savefig("Figures/gcurve_with_residuals.png", dpi=300, bbox_inches='tight')
plt.show()