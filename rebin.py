#!/usr/bin/env python3
from pathlib import Path
import numpy as np
from astropy.io import fits

# ---- edit this if you like ----
in_path = Path("/Users/haydentedrake/Desktop/Y2/12.410/Data/Given/master-r-flat-20250919.fits")
out_path = in_path.with_name(in_path.stem + "-2x2.fits")
normalize_to_median_one = True   # keep True for master flats
# -------------------------------

def bin2x2_sum(arr: np.ndarray) -> np.ndarray:
    """Bin 2D image by 2x2, summing electrons (emulates on-chip binning).
       Crops odd edge if needed."""
    if arr.ndim != 2:
        raise ValueError("Expected a 2D image in the primary HDU.")
    H, W = arr.shape
    H2, W2 = (H // 2) * 2, (W // 2) * 2  # ensure even dims
    cropped = arr[:H2, :W2]
    # reshape to (H/2, 2, W/2, 2) and sum over the 2x2 axes
    binned = cropped.reshape(H2 // 2, 2, W2 // 2, 2).sum(axis=(1, 3))
    return binned

def main():
    assert in_path.exists(), f"Input FITS not found: {in_path}"

    with fits.open(in_path) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header.copy()

    if data is None:
        raise RuntimeError("Primary HDU has no data.")

    data = np.asarray(data, dtype=np.float32)
    binned = bin2x2_sum(data)

    if normalize_to_median_one:
        med = float(np.nanmedian(binned))
        if med != 0 and np.isfinite(med):
            binned = binned / med
        else:
            print("Warning: median is zero/NaN/inf; skipping normalization.")
            med = np.nan
    else:
        med = np.nan

    # Update/insert helpful header cards
    hdr['HISTORY'] = "Rebinned 2x2 from 1x1 (sum of 2x2 blocks)."
    if normalize_to_median_one:
        hdr['HISTORY'] = "Renormalized to median = 1.0 after binning."
        hdr['NORMMED'] = (1.0, 'Post-binning normalization target median')
        hdr['PREVMED'] = (med, 'Median before normalization (ADU-sum units)')
    # Common binning keywords across cameras vary; set a few standard ones
    hdr['XBINNING'] = (2, 'Binning factor in X')
    hdr['YBINNING'] = (2, 'Binning factor in Y')
    hdr['BINNING']  = ('2 2', 'Binning factor (X Y)')
    # Clean shape-specific keywords if present (optional, safe to leave as-is)
    for k in ('NAXIS1', 'NAXIS2'):
        if k in hdr:
            del hdr[k]  # astropy will set these from the new data shape

    fits.writeto(out_path, binned.astype(np.float32), hdr, overwrite=True)
    print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
