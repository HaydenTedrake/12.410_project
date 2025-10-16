from astropy.io import fits
import numpy as np
from photutils.aperture import CircularAperture
import matplotlib.pyplot as plt

### QUESTION 5 ###
# open file
file = fits.open("421720251002/Light/4217_r'_1x1_180.000secs_1.fit")

# separate into data and header
data = file[0].data
header = file[0].header
print(header)

# make the plot with labels and a color map
plt.figure(figsize=(8,6))
im = plt.imshow(np.log10(data), origin='lower', cmap='viridis')
plt.xlabel("X [pixels]")
plt.ylabel("Y [pixels]")
plt.colorbar(im, label="log10(Counts)")
title = f"4217 - {header['FILTER']}-band, {header['EXPTIME']} s exposure, 2025-10-02"
plt.title(title)
plt.savefig('workshop_log.png')
plt.show()

file.close()