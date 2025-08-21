#!/usr/bin/env python3

import numpy as np
from scipy import integrate, optimize

from . import constants as const

class atmosphere:
    def __init__(self, system, T_obs, mu, kappa_opt=None, P_obs=None,
                 R_trunc=10, model='iso', gamma=5/3):
        """
        DEFUNCT CLASS - no longer used, but kept for future reference

        Arguments:
            system: planetary system object (see system.py)
            T_obs: Isothermal temperature of upper stratosphere (K)
            mu: Mean molecular mass of upper stratosphere (g)

        Keyword arguments:
            kappa_opt: Optical oppacity of upper stratosphere (cm^3/g)
            P_obs: Pressure at observed radius (erg/cm^3)
            R_trunc: Truncation radius of atmosphere (system.Rp)
            model: isothermal or plane-paralell
            gamma: Adiabatic index of gas
        """
        self.system = system
        self.T_obs = T_obs
        self.mu = mu
        self.R_trunc = R_trunc
        self.model = model
        self.gamma = gamma
        # isothermal sound speed
        self.cs2_iso = const.kB*self.T_obs/self.mu
        # C is phi0/h0, gravity over enthalpy
        self.C = (const.G*self.system.Mp/self.system.Rp)/self.cs2_iso
        if self.model == 'adia':
            self.C = (self.gamma-1)*self.C/self.gamma
            self.R_trunc = self.C/(self.C-1.)
        if kappa_opt is not None:
            self.kappa_opt = kappa_opt
            # calculates density at observed radius
            self.rho_obs = self.observed_density(self.kappa_opt)
            self.P_obs = self.rho_obs*self.cs2_iso
        elif P_obs is not None:
            self.P_obs = P_obs
            self.rho_obs = self.P_obs/self.cs2_iso
            # Reverse engineer kappa_opt
            self.kappa_opt = self.observed_density(1)/self.rho_obs
        else:
            print("ERROR: Unable to determine density at observed radius.\n"
                  "       Supply either a Pressure at the observed radius "
                  "(kwarg: P_obs) or\n"
                  "       supply the optical opacity of the upper stratosphere "
                  "(kwarg: kappa_opt)")
            return
        self.last_rho = None

        return


    def atmosphere_tuple(self):
        return (self.system, self.T_obs, self.mu, self.kappa_opt)


    def print_atmosphere(self):
        print('Atmosphere parameters (cgs)           # Normalized units\n'
              '  T_obs:     {:e} K           #\n'
              '  mu:        {:e} g           # %8.2f amu\n'
              '  rho_obs:   {:e} g/cm^3      #\n'
              '  P_obs:     {:e} erg/cm^3    # %8.2e bar\n'
              '  kappa_opt: {:e} cm^3/g      #\n'
              .format(self.T_obs, self.mu, self.rho_obs, self.P_obs,
                      self.kappa_opt)
              % (self.mu/const.amu, self.P_obs/const.bar))
        return


    def value(self, var=None):
        if var is None:
            return self.atmosphere_tuple()
        elif var == "T_obs":
            return self.T_obs
        elif var == "mu":
            return self.mu
        elif var == "rho_obs":
            return self.rho_obs
        elif var == "P_obs":
            return self.P_obs
        elif var == "kappa_opt":
            return self.kappa_opt
        elif var == "Rtrunc":
            return self.Rtrunc
        elif var == "gamma":
            return self.gamma
        elif var == "model":
            return self.model
        elif var == "C":
            return self.C
        elif var == "cs2_iso":
            return self.cs2_iso
        return


    def assign(self, var, value):
        if var == "T_obs":
            self.T_obs = value
            self.cs2_iso = const.kB*self.T_obs/self.mu
            self.C = (const.G*self.system.Mp/self.system.Rp)/self.cs2_iso
        elif var == "mu":
            self.mu = value
            self.cs2_iso = const.kB*self.T_obs/self.mu
            self.C = (const.G*self.system.Mp/self.system.Rp)/self.cs2_iso
        elif var == "P_obs":
            self.P_obs = value
            self.rho_obs = self.P_obs/self.cs2_iso
            self.kappa_opt = self.observed_density(1)/self.rho_obs
        elif var == "kappa_opt":
            self.kappa_opt = value
            self.rho_obs = self.observed_density(value)
            self.P_obs = self.rho_obs*self.cs2_iso
        else:
            print(f"ERROR: {var:s} is not a recognized variable for altering")
            return
        if self.last_rho is not None:
            self.windbase_radius(self.last_rho)

        return


    def deepcopy(self):
        new = self.__class__(np.copy.deepcopy(self.system),
                             np.copy.deepcopy(self.T_obs),
                             np.copy.deepcopy(self.mu),
                             np.copy.deepcopy(self.kappa_opt))
        return new


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


    def observed_density(self, kappa_opt):
        s_max = np.sqrt(self.R_trunc**2-1.0)
        s = np.linspace(0.0, s_max, 4096, endpoint=False)
        return (1/(2.*kappa_opt*self.system.Rp
                   *integrate.simpson(self.profile(np.sqrt(s**2+1)), s)))


    def windbase_radius(self, rho):
        """
        Description:
            Calculates the radius of the wind base by finding the radius at
            which rho_rmin matches onto our isothermal toy model atmosphere
        """
        self.last_rho = rho
        # Finds Rmin in units of system.Rp such that rho is in agreement
        self.Rmin = optimize.fsolve(self.rad_rho_diff, 0.9, args=(rho))[0]
        return


    def rad_rho_diff(self, r, rho):
        return self.rho_obs*self.profile(r)-rho