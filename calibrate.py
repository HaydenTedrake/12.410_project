from pathlib import Path
import numpy as np
from astropy.io import fits

date = "1009"

script_dir = Path(__file__).parent
dark_dir = script_dir / "Raw" / f"42172025{date}" / "Dark"

files = sorted(p for p in dark_dir.glob("*.fit*") if p.is_file())
assert files, f"No dark FITS files found in {dark_dir.resolve()}"

stack = []
for f in files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(np.float32)
        stack.append(data)
cube = np.stack(stack, axis=0)  # shape: (N, H, W)

master_dark = np.median(cube, axis=0)

hdr = fits.Header()
hdr['NCOMBINE'] = (len(files), 'number of frames combined')
hdr['COMBINE'] = ('median', 'combination method')
hdr['IMAGETYP'] = ('DARK', 'master dark')
hdr['DARKLVL'] = (float(np.nanmedian(master_dark)), 'median level of master dark (ADU)')
fits.writeto(f"master-dark-{date}.fits", master_dark.astype(np.float32), hdr, overwrite=True)
print(f"Wrote master-dark-{date}.fits")