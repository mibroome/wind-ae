import numpy as np
from scipy.special import wofz

from wind_ae.wrapper.wrapper_utils import constants as const


"""
CII 2s2 2p 2P1/2 - 2s 2p2 2D3/2 from the Chianti database
"""
_wl = 1334.532e-8 # cm
_nu = const.c/_wl
_f12 = 0.3107  # not sure if there's a factor of g missing
_gamma = 2.423e8 # 1/s


def CII_1_nu(nu, n_CII, T, mu):
    """
    Description:
        Calculates the Lyman-alpha attenuation coefficient of the gas due to
        neutral hydrogen. Takes into account the frequnecy dependence of the
        cross section, and thermal spontaneous emission.

    Arguments:
        nu: freuqnecy (in inverse seconds)
        n_CII: number density of CII (in inverse cubic centimeters)
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)
        
    Returns:
        Lyman-alpha attenuation coefficient
    """
    # Calculate intrinsic cross-section
    # N.B. this is not sigma at nu_0 since Voigt(nu_0) \neq 1.
    sigma_0 = _f12*const.pi*const.e**2./(const.me*const.c)
    # Calculate frequency dependent cross section
    sigma = sigma_0*CII_1_Voigt(n_CII, nu, T, mu)
    # Calculate attenuation coefficient
    alpha_nu = n_CII*sigma
    # Take spontaneous emission into account
    alpha_nu *= (1.-np.exp(-const.h*nu/(const.kB*T)))
    
    return alpha_nu

def CII_1_Voigt(n_CII, nu, T, mu):
    """
    Description:
        Calculates the Voigt profile, taking into consideration the Doppler,
        and nautral broadening of the line profile due to the properties of
        the gas. Uses a numerical evaluation of Faddeeva's function (wofz)
        to calculate the Voigt profile.

    Arguments:
        n_CII: number density of CII (in inverse cubic centimeters)
        nu: freuqnecy (in inverse seconds)
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)

    Returns:
        Voigt line profile given local properties of the gas
        
    Reference:
        1. https://en.wikipedia.org/wiki/Voigt_profile
        2. Draine's Physics of the ISM/IGM, pg. 58
    """
    # Doppler frequency shift
    x = float(nu-_nu)
    # Lorentzian broadening terms (natural)
    y = float(_gamma/(4.*const.pi))
    # Doppler broadening parameter
    sigma_rt2 = float(np.sqrt(2.*const.kB*T/mu)*_nu/const.c)
    z = (x+1j*y)/sigma_rt2
    return np.real(wofz(z))/(sigma_rt2*np.sqrt(const.pi))
