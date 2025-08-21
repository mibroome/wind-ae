#!/usr/bin/env python3

# Source paper: Linsky et al. 2014 (2014ApJ...780...61L)

import numpy as np

_band_edges = np.array([10, 20, 30, 40, 50, 60, 70, 80, 91.2, 117])
_band_widths = np.diff(_band_edges)
_band_centers = (_band_edges[1:]+_band_edges[:-1])/2

def f_uv_bins(f_Lya, Mstar=False):
    """
    Description:
        Uses Linsky et al. 2014 to calculate the uv spectrum binned into
        9 bands at 1 au given a stellar Lyman alpha flux (also at 1 au).
        Bands are [10-20], [20-30], [30-40], [40-50], [50-60], [60-70],
        [70-80], [80-91.2], and [91.2-117] nm.

    Arguments:
        f_Lya: The Lyman-alpha flux at 1 a.u. (in ergs/cm^2/s)
        
    Keyword arguments:
        Mstar: Boolean if the star is an M star
        
    Returns:
        uv spectrum in 9 bins at 1 au (in erg/cm^2/s)
    """
    f_Lya = np.asarray(f_Lya)
    m = np.array([0.344, 0.309, 0.000, 0.258,
                  0.572, 0.240, 0.518, 0.764, 0.065])
    b = np.array([1.357, 1.300, 0.882, 2.294,
                  2.098, 1.920, 1.894, 1.811, 1.004])
    if Mstar:
        m[:3] = [0.000, 0.000, 0.000]
        b[:3] = [0.491, 0.548, 0.602]
    if f_Lya.size > 1:
        return np.array([(10**-b)*f**(1+m) for f in f_Lya])
    else:
        return f_Lya*10**(-b+m*np.log10(f_Lya))


def f_uv(f_Lya, Mstar=False, wl_min=40, wl_max=91.2):
    """
    Description:
        Uses Linsky et al. 2014 to calculate a uv flux at 1 au given a
        Lyman alpha flux (also at 1 au).

    Arguments:
        f_Lya: The Lyman-alpha flux at 1 a.u. (in ergs/cm^2/s)
        
    Keyword arguments:
        Mstar: Boolean if the star is an M star
        wl_min: Lower bound of integrated spectrum (in nm)
        wl_max: Upper bound of integrated spectrum (in nm)
        
    Returns:
        Integrated uv flux at 1 au (in erg/cm^2/s)
    """
    if wl_min < _band_edges[0] or wl_max > _band_edges[-1]:
        print("ERROR: wl_min outside Linsky's domain [10-117] nm.")
        return
    bmin = next(b for b, edge in enumerate(_band_edges) if edge > wl_min)-1
    bmax = next(b for b, edge in enumerate(_band_edges) if edge >= wl_max)
    f_Lya = np.asarray(f_Lya)
    f_bins = f_uv_bins(f_Lya, Mstar)
    if f_Lya.size > 1:
        f_uv = [None,]*len(f_Lya)
        for i in range(len(f_Lya)):
            f_uv[i] = f_bins[i][bmin:bmax].sum()
            f_uv[i] -= (f_bins[i][bmin]
                        *(wl_min-_band_edges[bmin])/_band_widths[bmin])
            f_uv[i] -= (f_bins[i][bmax]
                        *(_band_edges[bmax]-wl_max)/_band_widths[bmax])
    else:
        f_uv = f_bins[bmin:bmax].sum()
        f_uv -= f_bins[bmin]*(wl_min-_band_edges[bmin])/_band_widths[bmin]
        f_uv -= f_bins[bmax]*(_band_edges[bmax]-wl_max)/_band_widths[bmax]
    return np.asarray(f_uv)