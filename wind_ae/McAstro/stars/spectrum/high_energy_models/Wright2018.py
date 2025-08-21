#!/usr/bin/env python3

# Source paper: Wright et al. 2018 (2018MNRAS.479.2351W)

import numpy as np

# Tuples of beta, log10(L_XoB_sat), Ro_sat (See end of section 3.1)
_new_results = (-2.3, -3.05, 0.14) # Wright et al. 2018
_old_results = (-2.7, -3.13, 0.16) # Wright et. al 2011
_canonical = (-2.0, -3.13, 0.13)

def xray_fraction(Mass, P_rot, use_results='2018'):
    """
    Description:
        Calculates the x-ray luminosity fraction of the total bolometric
        luminosity as a function of stellar mass and rotation rate.
        Assuming that x-ray luminosity is ROAST [0.1-2.4] keV band.
    
    Arguments:
        Mass: stellar mass (in solar masses)
        P_rot: rotation rate of star (in days)
        
    Keyword arguments:
        use_results: results to use: '2018', '2011', or 'canonical'

    Returns:
        X-ray luminosity in terms of the bolometric luminosity
    
    Source paper:
        Wright et al. 2018 (2018MNRAS.479.2351W)
    """
    if use_results == '2018':
        beta, L_XoB_sat, Ro_sat = _new_results
    elif use_results == '2011':
        beta, L_XoB_sat, Ro_sat = _old_results
    elif use_results == 'canonical':
        beta, L_XoB_sat, Ro_sat = _canonical
    else:
        print(f"ERROR: Unrecognized use_result: {use_results}.\n"
              "       Use: '2018', '2011', or 'canonical'")
        return
    Mass = np.asarray(Mass)
    P_rot = np.asarray(P_rot)
    tau = 10**(2.33-1.5*Mass+0.31*Mass**2)
    P_sat = Ro_sat*tau
    with np.errstate(all='ignore'): #ignore np.where internally complaints
        return (10**(L_XoB_sat)
                *np.where(P_rot <= P_sat, 1., (P_rot/P_sat)**beta))