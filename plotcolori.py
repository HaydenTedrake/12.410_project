import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phi_fit = np.linspace(0, 1, 600)
m_fitr = np.load('r_fitted_curve.npy')
m_fitg = np.load('g_fitted_curve.npy')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), 
                               gridspec_kw={'height_ratios': [2, 1]})

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

ax1.plot(phi_fit, m_fitr, '-', linewidth=1.5, alpha=0.8, label="r' fit", color=colors[0])
ax1.plot(phi_fit, m_fitg, '-', linewidth=1.5, alpha=0.8, label="g' fit", color=colors[1])

ax1.invert_yaxis()
ax1.set_xlim(0, 1)
ax1.set_ylabel("Apparent Magnitude")
ax1.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5)
plt.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)
ax1.set_title("Rotational color variation of 4217 Engelhardt (WAO 14-in data)")

color_curve = m_fitg - m_fitr
ax2.plot(phi_fit, color_curve, 'k-', linewidth=2, label="g'-r' color")
ax2.invert_yaxis()
ax2.set_xlim(0, 1)
ax2.set_xlabel("Rotational Phase")
ax2.set_ylabel("g'-r' Color")
ax2.legend(loc="upper right", frameon=True, handlelength=1.5, handletextpad=0.5)
ax2.grid(True, linestyle=':', linewidth=0.7, alpha=0.7)

period_text = f'Period: 3.066 $\\pm$ 0.001 hrs'
ax1.text(0.01, 0.02, period_text, 
         transform=ax1.transAxes,
         ha='left', va='bottom', alpha=0.9)

plt.tight_layout()
plt.savefig("Figures/colorcurve.png")
plt.show()