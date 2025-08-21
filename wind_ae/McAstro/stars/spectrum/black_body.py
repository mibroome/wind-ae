#!/usr/bin/env python3

import numpy as np
import scipy.integrate as integrate

from wind_ae.McAstro.utils import constants as const

class BlackBody():
    def __init__(self, Temp):
        self.Temp = Temp


    def spec_plot(self, ax, **kwargs):
        x = np.logspace(-10, 2, 5000)
        ax.plot(x, np.pi*_Planck_wl(x, self.Temp),
                label=(r'B$_\lambda$(%d K)' % (self.Temp)), **kwargs)


    def spec_integral(self, lob=1e-8, upb=912e-8):
        """
        Description:
            Integrates a black bodies spectrum to return the irradiance.
        
        Notes:
            ^Mimumum lower bound is 1e-8 cm.
            ^Assume uniform radiance from surface of black body, such that
             irradiance is pi*radiance.
        
        Keyword arguments:
            lob: lower bound (in centimeters)
            upb: upper bound (in centimeters)
            
        Returns:
            irradiance (in erg per second per square centimeter)
        """
        result = 0
        x_lob = 1e-8 if lob < 1e-8 else lob
        x_upb = 1e-6 if upb > 1e-6 else upb
        result += integrate.quad(lambda x: _Planck_fq(x, self.Temp),
                                 const.c/x_upb, const.c/x_lob)[0]
        while x_upb < upb:
            x_lob = x_upb
            x_upb *= 1e2
            x_upb = x_upb if upb > x_upb else upb
            result += integrate.quad(lambda x: _Planck_fq(x, self.Temp),
                                     const.c/x_upb, const.c/x_lob)[0]
        return np.pi*result


def _Planck_wl(wl, T):
    """
    Description:
        For a given temperature and wavelength returns the spectral radiance
        
    Arguments:
        wl: wavelength (in centimeters)
        T: temperature (in Kelvins)
        
    Returns:
        spectral radiance (in erg per second per cubic centimeter per
                           steradian)
    """
    B = 2*const.hc*const.c/wl**5
    return B*n_Bose_Einstein(const.hc/wl, 0, T)


def Planck_irradiance_wl(wl, T):
    """
    Description:
        For a given temperature and wavelength returns the spectral
        irradiance.
        
    Arguments:
        wl: wavelength (in centimeters)
        T: temperature (in Kelvins)
        
    Returns:
        spectral irradiance (in erg per second per cubic centimeter)
    """
    return np.pi*_Planck_wl(wl, T)


def _Planck_fq(nu, T):
    """
    Description:
        For a given temperature and wavelength returns the spectral radiance
        
    Arguments:
        nu: frequency (in Hertzs)
        T: temperature (in Kelvins)
        
    Returns:
        spectral radiance (in erg per second per square centimeter per Hertz
                           per steradian)
    """
    B = 2*const.h*nu**3/const.c**2
    return B*n_Bose_Einstein(const.h*nu, 0, T)


def n_Bose_Einstein(e, mu, T, g=1):
    """
    Description:
        For a system of a given temperature, calcualtes the expected number
        of particles in a state of a given energy, chemical potential, and
        degencery.
        
    Arguments:
        e: energy (in ergs)
        mu: chemical potential (in ergs)
        T: temperature (in Kelvins)
    
    Keyword arguments:
        g: degeneracy
        
    Returns:
        expected number of particles in a given state
    """
    x = (e-mu)/(const.kB*T)
    with np.errstate(all='ignore'): #ignore np.where internally complaints
        return np.where(x > 100, g*np.exp(-x),
                        np.where(x < 1e-10, g/x, g/(np.exp(x)-1.)))