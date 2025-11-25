# import numpy as np
# import pickle
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
# from astropy.stats import sigma_clip

# CSV_PATH = "r'.csv"
# JD0 = 2460937.567619
# P_HOURS = 3.0661
# gap_days = 0.50
# title = "r' light curve of 4217 Engelhardt (WAO 14-in data)"

# df = pd.read_csv(CSV_PATH)
# t = df["JD"].astype(float).to_numpy()
# m = df["Mag"].astype(float).to_numpy()
# e = df["MagErr"].astype(float).to_numpy()

# P_days = P_HOURS / 24.0
# phi = ((t - JD0) / P_days) % 1.0

# new_m = []
# dt = np.diff(t, prepend=t[0])
# night_id = np.zeros_like(t, dtype=int)
# for i in range(0, len(t)):
#     night_id[i] = night_id[i-1] + (dt[i] > gap_days)
#     if night_id[i] == 1:
#         new_m.append(m[i]-0.1295071438)
#     else:
#         new_m.append(m[i])
# new_m = np.array(new_m)

# def two_harmonic(phi, a0, a1, b1, a2, b2):
#     return (
#         a0
#         + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
#         + b1 * np.sin(2 * np.pi * phi)
#         + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
#         + b2 * np.sin(4 * np.pi * phi)
#     )

# def three_harmonic(phi, a0, a1, b1, a2, b2, a3, b3):
#     return (
#         a0
#         + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
#         + b1 * np.sin(2 * np.pi * phi)
#         + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
#         + b2 * np.sin(4 * np.pi * phi)
#         + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
#         + b3 * np.sin(6 * np.pi * phi)
#     )

# def four_harmonic(phi, a0, a1, b1, a2, b2, a3, b3, a4, b4):
#     return (
#         a0
#         + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
#         + b1 * np.sin(2 * np.pi * phi)
#         + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
#         + b2 * np.sin(4 * np.pi * phi)
#         + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
#         + b3 * np.sin(6 * np.pi * phi)
#         + a4 * np.cos(8 * np.pi * phi)   # 4th harmonic
#         + b4 * np.sin(8 * np.pi * phi)
#     )

# def five_harmonic(phi, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5):
#     return (
#         a0
#         + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
#         + b1 * np.sin(2 * np.pi * phi)
#         + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
#         + b2 * np.sin(4 * np.pi * phi)
#         + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
#         + b3 * np.sin(6 * np.pi * phi)
#         + a4 * np.cos(8 * np.pi * phi)   # 4th harmonic
#         + b4 * np.sin(8 * np.pi * phi)
#         + a5 * np.cos(10 * np.pi * phi)  # 5th harmonic
#         + b5 * np.sin(10 * np.pi * phi)
#     )

# # popt, pcov = curve_fit(
# #     two_harmonic,
# #     phi,
# #     new_m,
# #     sigma=e,
# #     absolute_sigma=True,
# #     p0=[np.mean(new_m), 0, 0, 0, 0]
# # )
# # popt, pcov = curve_fit(
# #     three_harmonic,
# #     phi,
# #     new_m,
# #     sigma=e,
# #     absolute_sigma=True,
# #     p0=[np.mean(new_m), 0, 0, 0, 0, 0, 0]
# # )
# popt, pcov = curve_fit(
#     four_harmonic,
#     phi,
#     new_m,
#     sigma=e,
#     absolute_sigma=True,
#     p0=[np.mean(new_m), 0, 0, 0, 0, 0, 0, 0, 0]
# )
# # popt, pcov = curve_fit(
# #     five_harmonic,
# #     phi,
# #     new_m,
# #     sigma=e,
# #     absolute_sigma=True,
# #     p0=[np.mean(new_m), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# # )

# print("Best-fit parameters:", popt)

# # residuals = new_m - two_harmonic(phi, *popt)
# # residuals = new_m - three_harmonic(phi, *popt)
# residuals = new_m - four_harmonic(phi, *popt)
# # residuals = new_m - five_harmonic(phi, *popt)
# chi2 = np.sum((residuals / e) ** 2)
# dof = len(new_m) - len(popt)
# print(f"Reduced χ² = {chi2/dof:.2f}")

# plt.figure(figsize=(9, 5.2))
# w = 1.0 / np.maximum(e, 1e-6)**2

# colors = [
#     "#0072B2",  # blue
#     "#CC79A7",  # purple/magenta
#     "#009E73",  # bluish green
#     "#E69F00",  # orange
# ]
# markers = ["o", "s", "^", "D"]
# dates = ["20250919 UT", "20251003 UT", "20251010 UT", "20251024 UT"]

# unique_nights = np.unique(night_id)
# for i, nid in enumerate(unique_nights):
#     sel = (night_id == nid)
#     plt.errorbar(
#         phi[sel], new_m[sel], yerr=e[sel],
#         fmt=markers[i % len(markers)],
#         markersize=4.5,
#         elinewidth=0.7,
#         capsize=0,
#         alpha=0.7,
#         color=colors[i % len(colors)],
#         linestyle='none',
#         label=f"{dates[nid]}",
#     )

# tenten = []
# tentwentyfour = []
# for nid in unique_nights:
#     sel = (night_id == nid)
#     if nid == 2:
#         tenten.extend(zip(t[sel], new_m[sel], e[sel]))
#     elif nid == 3:
#         tentwentyfour.extend(zip(t[sel], new_m[sel], e[sel]))
# with open("rtenten.pkl", "wb") as f:
#     pickle.dump(tenten, f)
# with open("rtentwentyfour.pkl", "wb") as f:
#     pickle.dump(tentwentyfour, f)

# phi_fit = np.linspace(0, 1, 600)
# # new_m_fit = two_harmonic(phi_fit, *popt)
# # new_m_fit = three_harmonic(phi_fit, *popt)
# new_m_fit = four_harmonic(phi_fit, *popt)
# # new_m_fit = five_harmonic(phi_fit, *popt)
# plt.plot(phi_fit, new_m_fit, 'k-', lw=2, label="4th Order Fit")

# mean_val = np.average(new_m, weights=w)
# variance = np.average((new_m - mean_val)**2, weights=w)
# weight_sum = np.sum(w)
# mean_err = np.sqrt(variance / weight_sum)

# plt.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
# plt.text(0.515, mean_val, f'Mean: {mean_val:.2f} $\\pm$ {mean_err:.2f}', 
#          ha='center', va='bottom', alpha=0.9, fontsize=12)

# period_text = f'Period: {P_HOURS} $\\pm$ 0.0002 hrs'
# plt.text(0.01, 0.02, period_text, 
#          transform=plt.gca().transAxes,
#          ha='left', va='bottom', alpha=0.9, fontsize=12)

# plt.gca().invert_yaxis()
# plt.xlim(0, 1)
# plt.xlabel("Rotational Phase", fontsize=12)
# plt.ylabel("Apparent magnitude (r')", fontsize=12)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12) 
# plt.title(title, fontsize=12)
# plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
# plt.legend(loc="upper left", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)
# plt.tight_layout()
# plt.savefig("Figures/rcurve.png")
# plt.show()

# plt.plot(phi, residuals, 'ko', markersize=4.5, alpha=0.7)
# plt.show()

import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip

CSV_PATH = "r'.csv"
JD0 = 2460937.567619
P_HOURS = 3.0661
gap_days = 0.50
title = "r' light curve of 4217 Engelhardt (WAO 14-in data)"

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
        new_m.append(m[i]-0.1295071438)
    else:
        new_m.append(m[i])
new_m = np.array(new_m)

def two_harmonic(phi, a0, a1, b1, a2, b2):
    return (
        a0
        + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
        + b1 * np.sin(2 * np.pi * phi)
        + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
        + b2 * np.sin(4 * np.pi * phi)
    )

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

def four_harmonic(phi, a0, a1, b1, a2, b2, a3, b3, a4, b4):
    return (
        a0
        + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
        + b1 * np.sin(2 * np.pi * phi)
        + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
        + b2 * np.sin(4 * np.pi * phi)
        + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
        + b3 * np.sin(6 * np.pi * phi)
        + a4 * np.cos(8 * np.pi * phi)   # 4th harmonic
        + b4 * np.sin(8 * np.pi * phi)
    )

def five_harmonic(phi, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5):
    return (
        a0
        + a1 * np.cos(2 * np.pi * phi)   # 1st harmonic
        + b1 * np.sin(2 * np.pi * phi)
        + a2 * np.cos(4 * np.pi * phi)   # 2nd harmonic
        + b2 * np.sin(4 * np.pi * phi)
        + a3 * np.cos(6 * np.pi * phi)   # 3rd harmonic
        + b3 * np.sin(6 * np.pi * phi)
        + a4 * np.cos(8 * np.pi * phi)   # 4th harmonic
        + b4 * np.sin(8 * np.pi * phi)
        + a5 * np.cos(10 * np.pi * phi)  # 5th harmonic
        + b5 * np.sin(10 * np.pi * phi)
    )

popt, pcov = curve_fit(
    four_harmonic,
    phi,
    new_m,
    sigma=e,
    absolute_sigma=True,
    p0=[np.mean(new_m), 0, 0, 0, 0, 0, 0, 0, 0]
)

print("Best-fit parameters:", popt)

residuals = new_m - four_harmonic(phi, *popt)

# Create subplots
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(9, 8),
    gridspec_kw={"height_ratios": [2, 1]}  # top : bottom
)

w = 1.0 / np.maximum(e, 1e-6)**2

colors = [
    "#0072B2",  # blue
    "#CC79A7",  # purple/magenta
    "#009E73",  # bluish green
    "#E69F00",  # orange
]
markers = ["o", "s", "^", "D"]
dates = ["20250919 UT", "20251003 UT", "20251010 UT", "20251024 UT"]

# Plot light curve on top subplot
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
    if nid == 2:
        tenten.extend(zip(t[sel], new_m[sel], e[sel]))
    elif nid == 3:
        tentwentyfour.extend(zip(t[sel], new_m[sel], e[sel]))
with open("rtenten.pkl", "wb") as f:
    pickle.dump(tenten, f)
with open("rtentwentyfour.pkl", "wb") as f:
    pickle.dump(tentwentyfour, f)

phi_fit = np.linspace(0, 1, 600)
new_m_fit = four_harmonic(phi_fit, *popt)
ax1.plot(phi_fit, new_m_fit, 'k-', lw=2, label="4th Order Fit")

mean_val = np.average(new_m, weights=w)
variance = np.average((new_m - mean_val)**2, weights=w)
weight_sum = np.sum(w)
mean_err = np.sqrt(variance / weight_sum)

ax1.axhline(mean_val, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
ax1.text(0.515, mean_val, f'Mean: {mean_val:.2f} $\\pm$ {mean_err:.2f}', 
         ha='center', va='bottom', alpha=0.9, fontsize=12)

period_text = f'Period: {P_HOURS} $\\pm$ 0.0002 hrs'
ax1.text(0.01, 0.02, period_text, 
         transform=ax1.transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

ax1.invert_yaxis()
ax1.set_xlim(0, 1)
fig.supylabel("Apparent magnitude (r')", fontsize=12)
ax1.set_title(title, fontsize=12)
ax1.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
ax1.legend(loc="upper left", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)

# Plot residuals on bottom subplot
for i, nid in enumerate(unique_nights):
    sel = (night_id == nid)
    ax2.errorbar(
        phi[sel], residuals[sel], yerr=e[sel],
        fmt=markers[i % len(markers)],
        markersize=4.5,
        elinewidth=0.7,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
    )

# Calculate mean and error of residuals
residual_mean = np.average(residuals, weights=w)
residual_variance = np.average((residuals - residual_mean)**2, weights=w)
residual_mean_err = np.sqrt(residual_variance / weight_sum)

ax2.axhline(0, color='k', linestyle='-', linewidth=1, alpha=0.8)
ax2.axhline(residual_mean, color='k', linestyle=':', linewidth=1.1, alpha=0.8)
xtext = 0.5
ytext = residual_mean
ax2.text(xtext, ytext,
         f'Mean: {-residual_mean:.2f} $\\pm$ {residual_mean_err:.2f}',
         ha='center', va='bottom', alpha=0.9, fontsize=12)
ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_title("Residuals", fontsize=12)
ax2.set_xlabel("Rotational Phase", fontsize=12)
ax2.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)

plt.tight_layout()
plt.savefig("Figures/rcurve_with_residuals.png")
plt.show()
