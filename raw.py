import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.gridspec import GridSpec
from astropy.io import fits
import matplotlib as mpl
from matplotlib.patches import Circle
import math

# ----------------------------------------------------------------------
# Global plotting style
# ----------------------------------------------------------------------
mpl.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
})

# ----------------------------------------------------------------------
# Input file and pointing info
# ----------------------------------------------------------------------
fits_file = "Raw/4217/4217_r'_180secs_1.fit"

# Asteroid is centered in the frame at:
#   RA  = 00 49 25.70
#   Dec = +15 33 05.5
ra_h = 0 + 49/60 + 25.70/3600.0
ra0_deg = ra_h * 15.0                         # central RA in degrees
dec0_deg = 15 + 33/60 + 5.5/3600.0            # central Dec in degrees

# Plate scale: 0.2"/px native, but 2x2 binning on this night -> 0.4"/px
platescale_arcsec = 0.4
platescale_deg = platescale_arcsec / 3600.0   # degrees per pixel

# ----------------------------------------------------------------------
# Load image
# ----------------------------------------------------------------------
data = fits.getdata(fits_file)
ny, nx = data.shape

# Assume asteroid is at the geometric center of the frame
x0 = (nx - 1) / 2.0
y0 = (ny - 1) / 2.0

# ----------------------------------------------------------------------
# Build approximate RA/Dec coordinates across the field (in degrees)
# ----------------------------------------------------------------------
dec_rad = math.radians(dec0_deg)
cos_dec = math.cos(dec_rad)

x = np.arange(nx)
y = np.arange(ny)

dx_deg = (x - x0) * platescale_deg / cos_dec   # RA offset in degrees
dy_deg = (y - y0) * platescale_deg            # Dec offset in degrees

ra_deg = ra0_deg + dx_deg
dec_deg = dec0_deg + dy_deg

ra_min_deg, ra_max_deg = ra_deg.min(), ra_deg.max()
dec_min_deg, dec_max_deg = dec_deg.min(), dec_deg.max()

# ----------------------------------------------------------------------
# Tick label formatters
# ----------------------------------------------------------------------
def ra_formatter(x_deg, pos):
    """Format RA tick (input in degrees) as hh:mm:ss."""
    ra_hours = (x_deg / 15.0) % 24.0

    h = int(ra_hours)
    m_frac = (ra_hours - h) * 60.0
    m = int(m_frac)
    s = int((m_frac - m) * 60.0 + 0.5)

    # carry rounding
    if s == 60:
        s = 0
        m += 1
    if m == 60:
        m = 0
        h = (h + 1) % 24

    return f"{h:02d}h{m:02d}m{s:02d}s"

def dec_formatter(y_deg, pos):
    """Format Dec tick (input in degrees) as dd:mm."""
    sign = "+" if y_deg >= 0 else "-"
    y_abs = abs(y_deg)
    d = int(y_abs)
    m = int((y_abs - d) * 60.0 + 0.5)
    if m == 60:
        d += 1
        m = 0
    return f"{sign}{d:02d}°{m:02d}'"

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(7.0, 4.5))  # wider figure
gs = GridSpec(1, 2, width_ratios=[20, 1], wspace=0.05)

ax = fig.add_subplot(gs[0, 0])   # main image axis
cax = fig.add_subplot(gs[0, 1])  # skinny colorbar axis

data_pos = np.where(data > 0, data, np.nan)

im = ax.imshow(
    data_pos,
    origin="lower",
    cmap="gray",
    norm=LogNorm(
        vmin=np.nanpercentile(data_pos, 5),
        vmax=np.nanpercentile(data_pos, 99.5)
    ),
    extent=[ra_max_deg, ra_min_deg, dec_min_deg, dec_max_deg]  # RA reversed: east left
)

# --------- RED CIRCLE AROUND THE ASTEROID (ADJUSTABLE) ---------

# Change these three numbers to tweak the circle:
radius_arcsec = 16.0          # circle radius in arcsec (try 8, 10, 12, ...)
d_ra_arcsec   = 95           # RA offset in arcsec (+ = to the left)
d_dec_arcsec  = -30           # Dec offset in arcsec (+ = up)

# Convert offsets to degrees
radius_deg  = radius_arcsec / 3600.0
d_ra_deg    = d_ra_arcsec  / 3600.0 / math.cos(math.radians(dec0_deg))
d_dec_deg   = d_dec_arcsec / 3600.0

# Final center of the circle in RA/Dec (degrees)
ra_center_deg  = ra0_deg  + d_ra_deg
dec_center_deg = dec0_deg + d_dec_deg

circ = Circle(
    (ra_center_deg, dec_center_deg),
    radius_deg,
    edgecolor="red",
    facecolor="none",
    linewidth=2,
)
ax.add_patch(circ)

# ---------------------------------------------------------------

cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Counts (log scale)")

ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
ax.yaxis.set_major_formatter(FuncFormatter(dec_formatter))
ax.xaxis.set_major_locator(MaxNLocator(4))
ax.yaxis.set_major_locator(MaxNLocator(4))

ax.set_xlabel("Right Ascension (hours)")
ax.set_ylabel("Declination (deg)")
ax.set_title("Raw $r'$ frame of (4217) Engelhardt, 180 s, 2025-09-19, 02:46 UT")

ax.set_aspect("equal")

fig.subplots_adjust(left=0.10, right=0.92, top=0.90, bottom=0.15)

fig.savefig("raw_20250919.png", dpi=300, bbox_inches="tight")
plt.show()