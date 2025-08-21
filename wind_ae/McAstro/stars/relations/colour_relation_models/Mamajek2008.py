#!/usr/bin/env python3

# Source paper: Mamajek & Hillenbrand 2008 (2008ApJ...687.1264M)

import numpy as np

from wind_ae.McAstro.utils import constants as const

def rotation_rate(BV0, age, verbose=False):
    """
    Description:
        For a given age, turns a B-V colour into a stellar rotation rate.
        
    Notes:
        ^Results were fitted for 0.5 < B-V < 0.9
        ^The fit has a 'color singularity' at 0.495, below which results are
         imaginary. We've patched to give 0, but probably should keep to
         the range of the data fitted
         
    Arguments:
        BV0: intrinsic B-V colour of the star (in magnitudes)
        age: age of star (in years)
                
    Returns:
        Stellar rotation rate (in seconds)
        
    Source paper:
        Mamajek & Hillenbrand 2008 (2008ApJ...687.1264M)
    """
    BV0 = np.asarray(BV0)
    age_Myr = np.asarray(age)/1e6
    Mamajek_warning = False
    if np.any(BV0 <= 0.495) or np.any(BV0 >= 0.9):
        Mamajek_warning = True
        if verbose:
            print('WARNING: When calculating stellar rotation rates some B-V\n'
                  '         colors where outside valid range [0.5, 0.9].\n'
                  '         Be wary of results that use those rotation rate.')
    with np.errstate(all='ignore'): #ignore np.where internally complaints
        return (Mamajek_warning,
                np.where(BV0 <= 0.495, 0,
                         0.407*const.day*(BV0-0.495)**(0.325)*age_Myr**(0.566)))


def stellar_age(BV0, P_rot, verbose=False):
    """
    Description:
        Solve's Mamajek's relationship for the age given rotation rate.

    Arguments:
        BV0: intrinsic B-V colour of the star (in magnitudes)
        P_rot: stellar rotation rate (in seconds)

    Returns:
        Stellar age (in years)
    """
    BV0 = np.asarray(BV0)
    P_rot = np.asarray(P_rot)
    Mamajek_warning = False
    if np.any(BV0 <= 0.495) or np.any(BV0 >= 0.9):
        Mamajek_warning = True
        if verbose:
            print('WARNING: When calculating stellar ages from rotation rates\n'
                  '         some B-V colors outside valid range [0.5, 0.9].\n'
                  '         Be wary of results that use those stellar ages.')
    with np.errstate(all='ignore'): #ignore np.where internally complaints
        return (Mamajek_warning,
                np.where(BV0 <= 0.495, 0,
                         (P_rot/(0.407*const.days)
                          *(BV0-0.495)**(-0.325))**(1./0.566)*1e6))