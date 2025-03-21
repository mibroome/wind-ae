#!/usr/bin/env python3

# Source paper: Eker et al. 2018 (2018MNRAS.479.5491E)

import numpy as np

from wind_ae.McAstro.utils import constants as const

# Stephan-Boltzmann constant scaled for luminosity equation in solar units,
# i.e., SB_solar * (R/Rsun)**2 * T**4 = L/Lsun, e.g., SB_solar * 5772**4 = 1
_SB_solar = 4.*np.pi*const.sig_SB*(const.Rsun**2/const.Lsun)

_M_cut = [0.179, 0.450, 0.720, 1.050, 2.400, 7.000, 31.000]
_L_shift = [0,
           2.028*np.log10(_M_cut[1])-0.976-(4.572*np.log10(_M_cut[1])-0.102),
           4.572*np.log10(_M_cut[2])-0.102-(5.743*np.log10(_M_cut[2])-0.007),
           5.743*np.log10(_M_cut[3])-0.007-(4.329*np.log10(_M_cut[3])+0.010),
           4.329*np.log10(_M_cut[4])+0.010-(3.967*np.log10(_M_cut[4])+0.093),
           3.967*np.log10(_M_cut[5])+0.093-(2.865*np.log10(_M_cut[5])+1.105)
           ]

def _MLR(mass, shift=False):
    """
    Description:
        For a given mass use a mass-luminosity relationship (MLR) to
        prescribe a luminosity. Uses Table 4 from Eker et al. 2018.
    
    Arguments:
        mass: mass of star (in solar masses)
        
    Keyword arguments:
        shift: shift fits to be continuous (boolean)
        
    Returns:
        Luminosity of star (in solar luminosities)
        
    Source paper:
        Eker et al. 2018 (2018MNRAS.479.5491E)
    """
    log_mass = np.log10(mass)
    if mass < _M_cut[0]:
        print("Error: Mass smaller than Eker's emperical limits–no data.")
        return -1
    elif mass <= _M_cut[1]:
        log_Luminosity = 2.028*log_mass - 0.976
        if shift:
            log_Luminosity += _L_shift[0]*(_M_cut[1]-mass)/(_M_cut[1]-_M_cut[0])
    elif mass <= _M_cut[2]:
        log_Luminosity = 4.572*log_mass - 0.102
        if shift:
            log_Luminosity += _L_shift[1]*(_M_cut[2]-mass)/(_M_cut[2]-_M_cut[1])
    elif mass <= _M_cut[3]:
        log_Luminosity = 5.743*log_mass - 0.007
        if shift:
            log_Luminosity += _L_shift[2]*(_M_cut[3]-mass)/(_M_cut[3]-_M_cut[2])
    elif mass <= _M_cut[4]:
        log_Luminosity = 4.329*log_mass + 0.010
        if shift:
            log_Luminosity += _L_shift[3]*(_M_cut[4]-mass)/(_M_cut[4]-_M_cut[3])
    elif mass <= _M_cut[5]:
        log_Luminosity = 3.967*log_mass + 0.093
        if shift:
            log_Luminosity += _L_shift[4]*(_M_cut[5]-mass)/(_M_cut[5]-_M_cut[4])
    elif mass <= _M_cut[6]:
        log_Luminosity = 2.865*log_mass + 1.105
        if shift:
            log_Luminosity += _L_shift[5]*(_M_cut[6]-mass)/(_M_cut[6]-_M_cut[5])
    else:
        print("Error: Mass greater than Eker's emperical limits–no data.")
        return -1
    return 10**log_Luminosity
MLR = np.vectorize(_MLR)

_L_cut = [np.log10(_MLR(_M_cut[0], shift=True)),
          np.log10(_MLR(_M_cut[1], shift=True)),
          np.log10(_MLR(_M_cut[2], shift=True)),
          np.log10(_MLR(_M_cut[3], shift=True)),
          np.log10(_MLR(_M_cut[4], shift=True)),
          np.log10(_MLR(_M_cut[5], shift=True)),
          np.log10(_MLR(_M_cut[6], shift=True))]

def _LMR(lum):
    log_Luminosity = np.log10(lum)
    if log_Luminosity < _L_cut[0]:
        print("Error: Luminosity smaller than Eker's emperical limits–no data.")
        return -1
    elif log_Luminosity <= _L_cut[1]:
        log_mass = (log_Luminosity+0.976)/2.028
    elif log_Luminosity <= _L_cut[2]:
        log_mass = (log_Luminosity+0.102)/4.572
    elif log_Luminosity <= _L_cut[3]:
        log_mass = (log_Luminosity+0.007)/5.743
    elif log_Luminosity <= _L_cut[4]:
        log_mass = (log_Luminosity-0.010)/4.329
    elif log_Luminosity <= _L_cut[5]:
        log_mass = (log_Luminosity-0.093)/3.967
    elif log_Luminosity <= _L_cut[6]:
        log_mass = (log_Luminosity-1.105)/2.865
    else:
        print(_L_cut)
        print("Error: Luminosity greater than Eker's emperical limits–no data.")
        return -1
    return 10**log_mass
LMR = np.vectorize(_LMR)
    

def _MRR(mass, shift=False):
    """
    Description:
        For a given mass use a mass-radius relationship (MRR) to prescribe a
        radius. Uses Table 5 from Eker et al. 2018.
    
    Arguments:
        mass: mass of star (in solar masses)
        
    Keyword arguments:
        shift: shift fits to be continuous (boolean)
    
    Returns:
        Radius of star (in solar radii)
        
    Source paper:
        Eker et al. 2018 (2018MNRAS.479.5491E)
    """
    if mass < 0.179:
        print("Error: Mass smaller than Eker's emperical limits–no data.")
        return -1
    elif mass <= 1.5:
        Radius = 0.438*mass**2 + 0.479*mass + 0.075
    elif mass <= 31:
        Radius = np.sqrt(MLR(mass, shift=shift)/
                         (_SB_solar*MTR(mass, shift=shift)**4))
    else:
        print("Error: Mass greater than Eker's emperical limits–no data.")
        return -1
    return Radius
MRR = np.vectorize(_MRR)

_T_shift = [0,
            (_MLR(1.5, shift=True)/(
                _SB_solar*(0.438*(1.5)**2+0.479*(1.5)+0.075)**2))**(0.25)-
            10**(-0.170*np.log10(1.5)**2+0.888*np.log10(1.5)+3.671)
           ]
def _MTR(mass, shift=False):
    """
    Description:
        For a given mass use a mass-temperature relationship (MTR) to
        prescribe an effective temperature. Uses Table 5 from Eker et al.
        2018.
    
    Arguments:
        mass: mass of star (in solar masses)
        
    Keyword arguments:
        shift: shift fits to be continuous (boolean)
        
    Returns:
        Temperature of star (in Kelvin)
    
    Source paper:
        Eker et al. 2018 (2018MNRAS.479.5491E)
    """
    if mass < 0.179:
        print("Error: Mass smaller than Eker's emperical limits–no data.")
        return -1
    elif mass <= 1.5:
        Temperature = (MLR(mass, shift=shift)/
                       (_SB_solar*MRR(mass, shift=shift)**2))**(0.25)
    elif mass <= 31:
        Temperature = -0.170*np.log10(mass)**2 + 0.888*np.log10(mass) + 3.671
        Temperature = 10**Temperature
        if shift:
            Temperature += _T_shift[1]*(31-mass)/(31-1.5)
    else:
        print("Error: Mass greater than Eker's emperical limits–no data.")
        return -1
    return Temperature
MTR = np.vectorize(_MTR)
