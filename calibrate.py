from pathlib import Path
import numpy as np
from astropy.io import fits

date = "1023"

script_dir = Path(__file__).parent
bias_dir = script_dir / "Raw" / f"2025{date}" / "Bias" / "1x1"

files = sorted(p for p in bias_dir.glob("*.fit*") if p.is_file())
assert files, f"No bias FITS files found in {bias_dir.resolve()}"

stack = []
for f in files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(np.float32)
        stack.append(data)
cube = np.stack(stack, axis=0)  # shape: (N, H, W)

master_bias = np.median(cube, axis=0)

hdr = fits.Header()
hdr['NCOMBINE'] = (len(files), 'number of frames combined')
hdr['COMBINE'] = ('median', 'combination method')
hdr['IMAGETYP'] = ('BIAS', 'master bias')
hdr['BIASLVL'] = (float(np.nanmedian(master_bias)), 'median level of master bias (ADU)')
fits.writeto(f"master-bias-{date}.fits", master_bias.astype(np.float32), hdr, overwrite=True)
print(f"Wrote master-bias-{date}.fits")