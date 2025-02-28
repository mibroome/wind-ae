#!/usr/bin/env python3

import numpy as np

from .high_energy_models import (Chadney2015, Jackson2012, Linsky2014,
                                 SanzForcada2011, Tu2015, Wright2018)
from .black_body import BlackBody

__all__ = []

def Chadney_euv_surface_flux(F_x, updated=False):
    return Chadney2015.euv_surface_flux(F_x, updated=updated)
Chadney_euv_surface_flux.__doc__ = (Chadney2015.euv_surface_flux.__doc__)
__all__ += ['Chadney_euv_surface_flux']


def Jackson_xray_fraction(BV0, stellar_age):
    return Jackson2012.xray_fraction(BV0, stellar_age)
Jackson_xray_fraction.__doc__ = (Jackson2012.xray_fraction.__doc__)
__all__ += ['Jackson_xray_fraction']


def SanzForcada_euv_luminosity(L_x):
    return SanzForcada2011.euv_luminosity(L_x)
SanzForcada_euv_luminosity.__doc__ = (
    SanzForcada2011.euv_luminosity.__doc__)
__all__ += ['SanzForcada_euv_luminosity']


def SanzForcada_xray_luminosity(L_bol, stellar_age):
    return SanzForcada2011.xray_luminosity(L_bol, stellar_age)
SanzForcada_xray_luminosity.__doc__ = (
    SanzForcada2011.xray_luminosity.__doc__)
__all__ += ['SanzForcada_xray_luminosity']


def Tu_xray_luminosity(rot_rate, stellar_age):
    return Tu2015.xray_luminosity(rot_rate, stellar_age)
Tu_xray_luminosity.__doc__ = (Tu2015.xray_luminosity.__doc__)
__all__ += ['Tu_xray_luminosity']


def Wright_xray_fraction(stellar_mass, rot_rate, use_results='2018'):
    return Wright2018.xray_fraction(stellar_mass, rot_rate,
                                    use_results=use_results)
Wright_xray_fraction.__doc__ = (Wright2018.xray_fraction.__doc__)
__all__ += ['Wright_xray_fraction']


def Linsky_euv(f_Lya):
    return Linsky2014.f_uv(f_Lya)
Linsky_euv.__doc__ = (Linsky2014.f_uv.__doc__)
__all__ += ['Linsky_euv']


def bb_euv_surface_flux(T_eff, lob=10e-7, upb=91.2e-7):
    """
    Description:
        Calculates the uv emission from a the blackbody spectrum. Models
        the blackbody temperature that emits the uv spectrum by fitting
        PHOENIX models of the given T_eff in the uv region.
        
    Arguments:
        T_eff: The effective temperature of the star (in K)
        
    Keyword arguments:
        lob: Lower bound of spectral integral (in cm)
        upb: Upper bound of spectral integral (in cm)
        
    Returns:
        The integrated uv flux at the stellar surface (in erg/cm^2/s)
    """
    T_eff = np.asarray(T_eff)
    # Fits for T_euv gotten from fitting PHOENIX models
    m0, b0, T0 = 0.655635, 1106.224278, 2300.000000
    m1, b1, T1 = 0.733707, 3204.257429, 5500.000000
    m2, b2, T2 = 0.915299, 5625.490934, 8800.000000
    m3, b3, T3 = 0.744415, 9012.095703, 12500.000000
    T_euv = np.where(T_eff < T0, b0*T_eff/T0,
                     np.where(T_eff < T1, m0*(T_eff-T0)+b0,
                              np.where(T_eff < T2, m1*(T_eff-T1)+b1,
                                       np.where(T_eff < T3, m2*(T_eff-T2)+b2,
                                                m3*(T_eff-T3)+b3))))
    Fuv = []
    print(T_eff, T_euv)
    if T_euv.size == 1:
        Fuv.append(BlackBody(T_euv).spec_integral(lob=lob, upb=upb))
    else:
        for Temp in T_euv:
            Fuv.append(BlackBody(Temp).spec_integral(lob=lob, upb=upb))
    return np.asarray(Fuv)