import numpy as np
from scipy.special import wofz

from wind_ae.wrapper.wrapper_utils import constants as const


"""
Mean values for Lyman-alpha doublet, i.e., marginalizing over both 
2P_{3/2} and 2P_{1/2} to 1S_{1/2} transitions
"""
_Lya_wl = 1215.67e-8 # cm
_Lya_nu = const.c/_Lya_wl
_Lya_f12 = .4164
_Lya_gamma = 6.265e8 # 1/s


def Lya_alpha_nu(nu, n_HI, T, mu, _nu, _f12, _A21, species=''):
    """
    Description:
        Calculates the Lyman-alpha attenuation coefficient of the gas due to
        neutral hydrogen. Takes into account the frequnecy dependence of the
        cross section, and thermal spontaneous emission.

    Arguments:
        nu: freuqnecy (in inverse seconds)
        n_HI: neutral hydrogen number density (in inverse cubic centimeters)
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)
        
    Returns:
        Lyman-alpha attenuation coefficient
    """
    # Calculate intrinsic cross-section
    # N.B. this is not sigma at nu_0 since Voigt(nu_0) \neq 1.
    sigma_0 = _Lya_f12*const.pi*const.e**2./(const.me*const.c)
    # Calculate frequency dependent cross section
    sigma = sigma_0*Lya_Voigt(n_HI, nu, T, mu)
    # Calculate attenuation coefficient
    alpha_nu = n_HI*sigma
    # Take spontaneous emission into account
    alpha_nu *= (1.-np.exp(-const.h*nu/(const.kB*T)))
    
    return alpha_nu

def Lya_Voigt(n_HI, nu, T, mu):
    """
    Description:
        Calculates the Voigt profile, taking into consideration the Doppler,
        and nautral broadening of the line profile due to the properties of
        the gas. Uses a numerical evaluation of Faddeeva's function (wofz)
        to calculate the Voigt profile.

    Arguments:
        n_HI: neutral hydrogen number density (in inverse cubic centimeters)
        nu: freuqnecy (in inverse seconds)
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)

    Returns:
        Lyman-alpha Voigt line profile given local properties of the gas
        
    Reference:
        1. https://en.wikipedia.org/wiki/Voigt_profile
        2. Draine's Physics of the ISM/IGM, pg. 58
    """
    # Doppler frequency shift
    x = float(nu-_Lya_nu)
    # Lorentzian broadening terms (natural)
    y = float(_Lya_gamma/(4.*const.pi))
    # Doppler broadening parameter
    sigma_rt2 = float(np.sqrt(2.*const.kB*T/mu)*_Lya_nu/const.c)
    z = (x+1j*y)/sigma_rt2
    return np.real(wofz(z))/(sigma_rt2*np.sqrt(const.pi))
