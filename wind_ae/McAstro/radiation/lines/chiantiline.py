import numpy as np
from scipy.special import wofz

from wind_ae.wrapper.wrapper_utils import constants as const


#_nu is the line center frequency s^-1


def line_alpha_nu(nu, n_abs, T, mu, _nu, _f12, _A21):
    """
    Description:
        Calculates the line attenuation coefficient of the gas due to
        an absorbing species. Takes into account the frequnecy dependence of the
        cross section, and thermal spontaneous emission.

    Arguments:
        nu: freuqnecy (in inverse seconds)
        n_sptot: total number density of species that produces the line (in inverse cubic centimeters)
                the code assigns level populations according to ionization equlibrium
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)
        
    Returns:
        line attenuation coefficient
    """
    # Calculate intrinsic cross-section
    # N.B. this is not sigma at nu_0 since Voigt(nu_0) \neq 1.
    sigma_0 = _f12*const.pi*const.e**2./(const.me*const.c)
    # Calculate frequency dependent cross section
    sigma = sigma_0*get_Voigt(nu, T, mu, _nu, _f12, _A21)
    # Calculate attenuation coefficient
    alpha_nu = n_abs*sigma
    # Take spontaneous emission into account
    alpha_nu *= (1.-np.exp(-const.h*nu/(const.kB*T)))
    
    return alpha_nu

def get_Voigt(nu, T, mu, _nu, _f12, _A21):
    """
    Description:
        Calculates the Voigt profile, taking into consideration the Doppler,
        and nautral broadening of the line profile due to the properties of
        the gas. Uses a numerical evaluation of Faddeeva's function (wofz)
        to calculate the Voigt profile.

    Arguments:
        nu: frequency (in inverse seconds)
        T: temperature (in Kelvin)
        mu: mean molecular weight (in grams)

    Returns:
        Voigt line profile given local properties of the gas
        
    Reference:
        1. https://en.wikipedia.org/wiki/Voigt_profile
        2. Draine's Physics of the ISM/IGM, pg. 58
    """

    _gamma = _A21 # just natural broadening for now
    
    # Doppler frequency shift
    x = float(nu-_nu)
    # Lorentzian broadening terms (natural)
    y = float(_gamma/(4.*const.pi))
    # Doppler broadening parameter
    sigma_rt2 = float(np.sqrt(2.*const.kB*T/mu)*_nu/const.c)
    #print(sigma_rt2, _nu)
    if(sigma_rt2 == 0): print(sigma_rt2, T, _nu)
    z = (x+1j*y)/sigma_rt2
    return np.real(wofz(z))/(sigma_rt2*np.sqrt(const.pi))
