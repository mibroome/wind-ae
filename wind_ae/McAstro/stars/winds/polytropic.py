#!/bin/env python3

import numpy as np
from scipy.optimize import fsolve

class polytropic():
    def __init__(self, Gamma, rad_crit, vel_crit,
                 x_npts=101, x_lo=-2, x_hi=2):
        self.Gamma = Gamma
        self.rad_crit = rad_crit
        self.vel_crit = vel_crit
        
        self.x = np.logspace(x_lo, x_hi, x_npts)
        self.rad = self.x*self.rad_crit
        
        self.polytropic_profile()


    def polytropic_profile(self):
        self.w = polytropic_vel_norm(self.x, self.Gamma)
        self.vel = self.w*self.vel_crit
        self.T = self.w**(1.-self.Gamma) * self.x**(2.-2.*self.Gamma)
        self.Mach = self.w/np.sqrt(self.T)
        
    
def _polytropic_momentum_eq(w, x, Gamma):
    """
    Equation 4.52 of Introduction to Stellar Winds (Lamers & Cassinelli)
    """
    # Avoid imaginary solutions (vel > 0 everywhere)
    if w <= 0:
        w = 0
    zero = (w**(Gamma+1.)
            -w**(Gamma-1.)*(4./x+(5.-3.*Gamma)/(Gamma-1.))
            +x**(2.-2.*Gamma)*(2./(Gamma-1.)))
    return zero


def _polytropic_vel_norm(x, Gamma):
    w0 = 0
    if x >= 1.0:
        w0 = 1.0
    return fsolve(_polytropic_momentum_eq, w0, args=(x, Gamma))
polytropic_vel_norm = np.vectorize(_polytropic_vel_norm)