#!/usr/bin/env python3

import numpy as np
from scipy import integrate, optimize

import wind_ae.McAstro.utils.constants as const

class planet_atmosphere:
    """
    Arguments:
        Mp: Mass of planet
        R_obs: Observed optical radius
        T_skin: Skin (equilibrium) temperature
        
    Keyword arguments:
        R_trunc: Truncation radius of atmosphere (in units of R_obs)
        model: isothermal or plane-paralell
    """
    def __init__(self, M_p, R_obs, T_skin, mu, R_trunc=10, model='iso',
                 gamma=5/3):
        self.M_p = M_p
        self.R_obs = R_obs
        self.T_skin = T_skin
        self.mu = mu
        self.R_trunc = R_trunc
        self.model = model
        self.gamma = gamma
        
        self.cs2_iso = const.kB*self.T_skin/self.mu
        # C is phi0/h0, gravity over enthalpy
        self.C = (const.G*self.M_p/self.R_obs)/self.cs2_iso
        
        if self.model == 'adia':
            self.C = (self.gamma-1)*self.C/self.gamma
            self.R_trunc = self.C/(self.C-1.)
        
        return
        
        
    def profile(self, r):
        if self.model == 'iso':
            power = -self.C*(1.-1./r)
        elif self.model == 'pp':
            power = -self.C*(r-1.)
        elif self.model == 'adia':
            return (1.+self.C*(1./r-1.))**(1./(self.gamma-1.))
        else:
            SystemExit('Profile kind not known.')
            
        return np.exp(power)


    def column(self, r):
        x = np.linspace(float(r), self.R_trunc, 4096, endpoint=False)
        N = self.R_obs*integrate.simps(self.profile(x), x)
        N *= self.rho_obs/self.mu
        
        return N
    
    
    def tau_1(self, x, kappa):
        tau = self.mu*kappa*self.column(x)
        
        return tau-1
    
    
    def dtau(self, x, kappa):
        rho = self.rho_obs*self.profile(x)
        dtau = -kappa*rho*self.R_obs
        
        return dtau


    def planet_radius(self, kappa_opt, kappa_uv):
        s_max = np.sqrt(self.R_trunc**2-1.0)
        s = np.linspace(0.0, s_max, 4096, endpoint=False)
        
        self.rho_obs = 1/(2.*kappa_opt*self.R_obs
                          *integrate.simps(self.profile(np.sqrt(s**2+1)), s))
        self.R_opt = optimize.fsolve(self.tau_1, 0.9, args=(kappa_opt),
                                     fprime=self.dtau)[0]
        self.R_uv = optimize.fsolve(self.tau_1, 0.9, args=(kappa_uv),
                                    fprime=self.dtau)[0]
        
        return