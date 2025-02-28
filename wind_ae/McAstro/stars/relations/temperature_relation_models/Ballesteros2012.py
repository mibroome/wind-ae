#!/usr/bin/env python3

# Source paper: Ballesteros 2012 (2012EL.....9734008B)

import numpy as np

def Temp2BV(T_eff):
    """
    Description:
        Given the effective temperature of the star, return a B-V colour.
        Uses quasi-analytic black-body calculations.
         
    Arguments:
        T_eff: effective temperature of star (in Kelvin)
        
    Returns:
        B-V colour of the star (in magnitudes)
        
    Source paper:
        Ballesteros 2012 (2012EL.....9734008B)
    """
    a = 0.92
    b = 1.70
    c = 0.62
    d = np.asarray(T_eff)/4600
    # Inversion of Equation 14 of Ballesteros 2012
    return (2.-d*(b+c)+np.sqrt(4.+(b-c)**2*d**2))/(2*a*d)