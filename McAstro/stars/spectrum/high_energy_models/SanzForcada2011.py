#!/usr/bin/env python3

# Source paper: Sanz-Forcada et al. 2011 (2011A&A...532A...6S)

import numpy as np

def euv_luminosity(L_x):
    """
    Description:
        Uses Equation 3 of Sanz-Forcada et al. 2011 to calculate an euv
        luminosity given the x-ray luminosity.
        
    Notes:
        ^euv defined as wavelengths from 100–920 Angstrom
    
    Arguments:
        L_x: x-ray luminosity of star (in erg/s)

    Returns:
        euv luminosity (in erg/s)
    
    Source paper:
        Sanz-Forcada et al. 2011 (2011A&A...532A...6S)
    """
    L_euv = 10**(4.8 + 0.86*np.log10(L_x))
    return L_euv


def xray_luminosity(L_bol, stellar_age):
    """
    Description:
        Uses Equation 5 of Sanz-Forcada et al. 2011 to calculate a x-ray
        luminosity given the bolometric luminosity and age of the star.
        
    Notes:
        ^The relation was calibrated with late F to early M dwarfs.
        ^X-ray defined as wavelengths from 5–100 Angstrom
    
    Arguments:
        L_bol: bolometric luminosity of star (in erg/s)
        stellar_age: age of star (in years)

    Returns:
        x-ray luminosity (in erg/s)
    
    Source paper:
        Sanz-Forcada et al. 2011 (2011A&A...532A...6S)
            Credited to Garcés et al. in prep, but appears to never have
            been published
    """
    L_bol = np.asarray(L_bol)
    stellar_age = np.asarray(stellar_age)
    tau = stellar_age/1e9
    tau_i = 2.03e20/(L_bol**(0.65))
    return np.where(tau < tau_i, 6.3e-4*L_bol, 1.89e28/(tau**(-1.55)))