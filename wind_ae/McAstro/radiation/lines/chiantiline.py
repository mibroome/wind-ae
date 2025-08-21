import numpy as np
from scipy.special import wofz

from wind_ae.wrapper.wrapper_utils import constants as const
import wind_ae.McAstro.atoms.atomic_species as McAtom


#_nu is the line center frequency s^-1


def continuum_alpha_nu(sim,rs):
    """
    Description:
        Calculates the continuum attenuation coefficient of the gas due to
        an absorbing species. Assumes that photoionization is the dominant 
        absorber of photons, thus the absorption cross section is the
        photoionization sigma.

    Arguments:
        nu: frequency (in inverse seconds)
        n_abs: total number density of species that produces the line (in inverse cubic centimeters)
                the code assigns level populations according to ionization equlibrium
        T: temperature (in Kelvin)
        
    Returns:
        continuum attenuation coefficient
    """

    # Calculate frequency dependent ionization cross section
    spaced_list = McAtom.formatting_species_list(sim.windsoln.species_list)
    alpha_nu = np.zeros_like(nu)
    for sp,species_name in enumerate(spaced_list):

        n0_fit = CubicSpline(sim.windsoln.soln_norm['r'],sim.windsoln.soln['n_'+species_name.replace(' ','')])
        for j in range(len(alpha_nu)):
            # attenuated alpha_nu (no spontaneous emission term)
            alpha_nu += sim.windsoln.sigma_wl[j,sp]*n0_fit(rs)
    # # Take spontaneous emission into account -- NEED?
    # alpha_nu *= (1.-np.exp(-const.h*nu/(const.kB*T)))
    
    return alpha_nu #as a function of frequency and radius


def line_alpha_nu(nu, n_abs, T, _mu, _nu, _f12, _A21, _species):
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
