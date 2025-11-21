# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
# from astropy.stats import sigma_clip
# import pickle

# PATH1 = "1010rcomp.csv"
# PATH2 = "1010gcomp.csv"
# PATH3 = "1024rcomp.csv"
# PATH4 = "1024gcomp.csv"

# df1 = pd.read_csv(PATH1)
# df2 = pd.read_csv(PATH2)
# df3 = pd.read_csv(PATH3)
# df4 = pd.read_csv(PATH4)

# t1 = df1["JD"].astype(float).to_numpy()
# t2 = df2["JD"].astype(float).to_numpy()
# t3 = df3["JD"].astype(float).to_numpy()
# t4 = df4["JD"].astype(float).to_numpy()

# m1 = df1["Mag"].astype(float).to_numpy()
# m2 = df2["Mag"].astype(float).to_numpy()
# m3 = df3["Mag"].astype(float).to_numpy()
# m4 = df4["Mag"].astype(float).to_numpy()

# e1 = df1["MagErr"].astype(float).to_numpy()
# e2 = df2["MagErr"].astype(float).to_numpy()
# e3 = df3["MagErr"].astype(float).to_numpy()
# e4 = df4["MagErr"].astype(float).to_numpy()

# t = []
# m = []
# e = []
# for i in range(len(t2)):
#     t.append((t1[i]+t2[i])/2)
#     m.append(m1[i]-m2[i])
#     e.append(np.sqrt(e1[i]**2+e2[i]**2))
# for i in range(len(t3)):
#     t.append((t3[i]+t4[i])/2)
#     m.append(m3[i]-m4[i])
#     e.append(np.sqrt(e3[i]**2+e4[i]**2))
# t = np.array(t)
# m = np.array(m)
# e = np.array(e)
# print(t)

# JD0 = 2460958.5558045
# P_HOURS = 3.066
# P_days = P_HOURS / 24.0
# gap_days = 0.50
# phi = ((t - JD0) / P_days) % 1.0
# title = "g'-r' color index of comparison stars"

# dt = np.diff(t, prepend=t[0])
# night_id = np.zeros_like(t, dtype=int)
# for i in range(1, len(t)):
#     night_id[i] = night_id[i-1] + (dt[i] > gap_days)

# plt.figure(figsize=(9, 5.2))
# w = 1.0 / np.maximum(e, 1e-6)**2

# colors = [
#     "#009E73",  # bluish green
#     "#E69F00",  # orange
# ]
# markers = ["o", "s", "^", "D"]
# dates = ["20251010 UT", "20251024 UT"]

# unique_nights = np.unique(night_id)
# for i, nid in enumerate(unique_nights):
#     sel = (night_id == nid)
#     w_night = w[sel]
#     m_night = m[sel]
    
#     mean_night = np.average(m_night, weights=w_night)
#     var_night = np.average((m_night - mean_night)**2, weights=w_night)
#     weight_sum_night = np.sum(w_night)
#     mean_err_night = np.sqrt(var_night / weight_sum_night)
    
#     plt.axhline(mean_night, color='k', linestyle='-', linewidth=1.1, alpha=0.8)
#     if nid == 0:
#         plt.text(0.5, mean_night + 0.03, f'Mean: {mean_night:.2f} $\\pm$ {mean_err_night:.2f}', 
#              ha='center', va='bottom', alpha=0.9, fontsize=12)
#     else:
#         plt.text(0.5, mean_night, f'Mean: {mean_night:.2f} $\\pm$ {mean_err_night:.2f}', 
#              ha='center', va='bottom', alpha=0.9, fontsize=12)
    
#     plt.errorbar(
#         phi[sel], m[sel], yerr=e[sel],
#         fmt=markers[i % len(markers)],
#         markersize=4.5,
#         elinewidth=0.7,
#         capsize=0,
#         alpha=0.7,
#         color=colors[i % len(colors)],
#         linestyle='none',
#         label=f"Comp star {nid+1} ({dates[nid]})",
#     )

# period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'
# plt.text(0.01, 0.02, period_text, 
#          transform=plt.gca().transAxes,
#          ha='left', va='bottom', alpha=0.9, fontsize=12)

# plt.gca().invert_yaxis()
# plt.xlim(0, 1)
# plt.xlabel("Rotational Phase of 4217 Engelhardt", fontsize=12)
# plt.ylabel("Apparent magnitude (g'-r')", fontsize=12)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12) 
# plt.title(title, fontsize=12)
# plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
# plt.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)
# plt.tight_layout()
# plt.savefig("Figures/colorcurve.png")
# plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip
import pickle

# =========================
# 1. ASTEROID COLOR CURVE
# =========================

with open("rtenten.pkl", "rb") as f:
    rtenten = pickle.load(f)
with open("rtentwentyfour.pkl", "rb") as f:
    rtentwentyfour = pickle.load(f)
with open("gtenten.pkl", "rb") as f:
    gtenten = pickle.load(f)
with open("gtentwentyfour.pkl", "rb") as f:
    gtentwentyfour = pickle.load(f)
    gtentwentyfour = gtentwentyfour[:13]

t_ast = []
m_ast = []
e_ast = []

for i in range(len(rtenten)):
    t_ast.append((rtenten[i][0] + gtenten[i][0]) / 2)
    m_ast.append(gtenten[i][1] - rtenten[i][1])
    e_ast.append(np.sqrt(rtenten[i][2]**2 + gtenten[i][2]**2))
for i in range(len(rtentwentyfour)):
    t_ast.append((rtentwentyfour[i][0] + gtentwentyfour[i][0]) / 2)
    m_ast.append(gtentwentyfour[i][1] - rtentwentyfour[i][1])
    e_ast.append(np.sqrt(rtentwentyfour[i][2]**2 + gtentwentyfour[i][2]**2))

t_ast = np.array(t_ast)
m_ast = np.array(m_ast)
e_ast = np.array(e_ast)

P_HOURS = 3.066
P_days = P_HOURS / 24.0
gap_days = 0.50
JD0_ast = 2460958.55063  # your asteroid JD0
phi_ast = ((t_ast - JD0_ast) / P_days) % 1.0

dt_ast = np.diff(t_ast, prepend=t_ast[0])
night_id_ast = np.zeros_like(t_ast, dtype=int)
for i in range(1, len(t_ast)):
    night_id_ast[i] = night_id_ast[i-1] + (dt_ast[i] > gap_days)

w_ast = 1.0 / np.maximum(e_ast, 1e-6)**2

colors = [
    "#009E73",  # bluish green
    "#E69F00",  # orange
]
markers = ["o", "s", "^", "D"]
dates = ["20251010 UT", "20251024 UT"]

# Global mean for asteroid
mean_ast = np.average(m_ast, weights=w_ast)
var_ast = np.average((m_ast - mean_ast)**2, weights=w_ast)
weight_sum_ast = np.sum(w_ast)
mean_err_ast = np.sqrt(var_ast / weight_sum_ast)

period_text = f'Period: {P_HOURS:.3f} $\\pm$ 0.001 hrs'

# =========================
# 2. COMPARISON STAR COLOR CURVE
# =========================

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

t_comp = []
m_comp = []
e_comp = []

for i in range(len(t2)):
    t_comp.append((t1[i] + t2[i]) / 2)
    m_comp.append(m1[i] - m2[i])
    e_comp.append(np.sqrt(e1[i]**2 + e2[i]**2))
for i in range(len(t3)):
    t_comp.append((t3[i] + t4[i]) / 2)
    m_comp.append(m3[i] - m4[i])
    e_comp.append(np.sqrt(e3[i]**2 + e4[i]**2))

t_comp = np.array(t_comp)
m_comp = np.array(m_comp)
e_comp = np.array(e_comp)

JD0_comp = 2460958.5558045  # your JD0 for comp stars (use same as asteroid if you prefer)
phi_comp = ((t_comp - JD0_comp) / P_days) % 1.0

dt_comp = np.diff(t_comp, prepend=t_comp[0])
night_id_comp = np.zeros_like(t_comp, dtype=int)
for i in range(1, len(t_comp)):
    night_id_comp[i] = night_id_comp[i-1] + (dt_comp[i] > gap_days)

w_comp = 1.0 / np.maximum(e_comp, 1e-6)**2

# =========================
# 3. MAKE THE SUBPLOTS
# =========================

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(9, 8),
    gridspec_kw={"height_ratios": [2, 1]}  # top : bottom
)

# --- Top: 4217 Engelhardt color curve ---
unique_nights_ast = np.unique(night_id_ast)
for i, nid in enumerate(unique_nights_ast):
    sel = (night_id_ast == nid)
    ax1.errorbar(
        phi_ast[sel], m_ast[sel], yerr=e_ast[sel],
        fmt=markers[i % len(markers)],
        markersize=4.5,
        elinewidth=0.7,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
        label=f"{dates[nid]}",
    )

ax1.axhline(mean_ast, color='k', linestyle='-', linewidth=1.1, alpha=0.8)
ax1.text(0.5, mean_ast + 0.017, f'Mean: {mean_ast:.2f} $\\pm$ {mean_err_ast:.2f}',
         ha='center', va='bottom', alpha=0.9, fontsize=12)

ax1.text(0.01, 0.02, period_text,
         transform=ax1.transAxes,
         ha='left', va='bottom', alpha=0.9, fontsize=12)

ax1.invert_yaxis()
ax1.set_xlim(0, 1)
ax1.set_title("4217 Engelhardt", fontsize=12)
ax1.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
ax1.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)

# --- Bottom: comparison star color curve ---
unique_nights_comp = np.unique(night_id_comp)
for i, nid in enumerate(unique_nights_comp):
    sel = (night_id_comp == nid)
    w_night = w_comp[sel]
    m_night = m_comp[sel]

    mean_night = np.average(m_night, weights=w_night)
    var_night = np.average((m_night - mean_night)**2, weights=w_night)
    weight_sum_night = np.sum(w_night)
    mean_err_night = np.sqrt(var_night / weight_sum_night)

    # mean line for that night's comp-star color
    ax2.axhline(mean_night, color='k', linestyle='-', linewidth=1.1, alpha=0.8)

    # offset labels slightly so they don't overlap
    xtext = 0.5
    ytext = mean_night + (0.05 if nid == 0 else 0.07)
    ax2.text(xtext, ytext,
             f'Mean: {mean_night:.2f} $\\pm$ {mean_err_night:.2f}',
             ha='center', va='bottom', alpha=0.9, fontsize=12)

    ax2.errorbar(
        phi_comp[sel], m_comp[sel], yerr=e_comp[sel],
        fmt=markers[i % len(markers)],
        markersize=4.5,
        elinewidth=0.7,
        capsize=0,
        alpha=0.7,
        color=colors[i % len(colors)],
        linestyle='none',
        label=f"Comp star {nid+1}",
    )

ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_xlabel("Rotational Phase", fontsize=12)
ax2.set_title("Comparison Stars (phased using 4217 Engelhardt period)", fontsize=12)
ax2.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
ax2.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5, fontsize=12)

fig.suptitle("g'-r' color index (WAO 14-in data)",
             fontsize=12, x=0.535)
fig.supylabel("Apparent magnitude (g'-r')", fontsize=12)
plt.tight_layout()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("Figures/colorcurve_with_comp.png", dpi=300)
plt.show()