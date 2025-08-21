#!/bin/env python3

import numpy as np
from scipy.optimize import fsolve

from wind_ae.McAstro.utils import constants as const

class isentropic:
    def __init__(self, Mstar, R0, T0, v0, n0, mu0, gamma=5/3):
        self.Mstar = Mstar
        self.R0 = R0
        self.T0 = T0
        self.v0 = v0
        self.n0 = n0
        self.mu0 = mu0
        self.gamma = gamma
        
        self.b_min = self.bern_const(self.R0, self.v0, self.T0)
        r = np.linspace(1, 50, 1000)*R0
        self.vel = fsolve(self.usolve, v0, args=(r))
        return


    def bern_const(self, r, v, T):
        """
        Assumes constant ionization
        """
        bern = (0.5*v**2
                +self.gamma/(self.gamma-1)*const.kB*T/self.mu0
                -const.G*self.Mstar/r)
        return bern


    def usolve(self, u, r):
        return (0.5*u**2
                +self.gamma/(self.gamma-1)*const.kB*self.T0/self.mu0
                *((self.v0*self.R0**2)/(u*r))**(self.gamma-1.)
                -(const.G*self.Mstar/r +self.b_min))