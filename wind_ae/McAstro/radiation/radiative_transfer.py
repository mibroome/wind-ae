import numpy as np

from wind_ae.wrapper.wrapper_utils import constants as const
import wind_ae.McAstro.atoms.atomic_species as McAtom


class rt_ray:
    """
    Ray used for radiative transfer
    
    Attributes:
        I: spectral irradiance (in in ergs per second per square centimeters
                                per Hertz)
        v_LOS: line of sight velocity (in centimeters per second)
        n_abs: number density of absorbers (in inverse cubic centimeters)
        ds: path length (in centimeters)
        T: temperature (in Kelvins)
        mu: mean molecular weight (in grams per cubic centimeter)
    """
    def __init__(self, n_cells):
        """
        Arguments:
            n_cells: number of cells in ray
        """
        self.I = np.zeros(n_cells+1)
        self.v_LOS = np.zeros(n_cells)
        self.n_abs = np.zeros(n_cells)
        self.ds = np.ones(n_cells)
        self.T = np.zeros(n_cells)
        self.mu = np.zeros(n_cells)
        self.nu0 = 0.
        self.f12 = 0.
        self.A21 = 0.

def radiative_transfer(nu_ref, ray, alpha_func, species='', break_tau=230):
    """
    Description:
        Calculates the total optical depth along the chosen ray. User will
        supply a function for calculating the attenuation coefficient, and
        the ray object contains 1-D arrays of all required variables to
        accurately perform radiative transfer.

    Arguments:
        nu_ref: referene frequency (in inverse seconds)
        ray: object ray along which we perform radiative transfer
        alpha_nu: function which determines the attenuation coefficient
        
    Keyword arguments:
        break_tau: optical depth after which we can safely stop integrating
        
    Note:
        ^Also solves for the intensity along the ray
        
    Returns:
        Total optical depth along ray
    """
    s_tau = 0.0
    print('raylenght',len(ray.ds))

    for t in range(len(ray.ds)):
        cell_nu = shifted_nu(nu_ref, ray.v_LOS[t])
        cell_alpha_nu = alpha_func(cell_nu, ray.n_abs[t], ray.T[t], ray.mu[t], ray.nu0, ray.f12, ray.A21, species)
        print(np.shape(cell_alpha_nu))
        ray.I[t+1], dtau = rad_step(cell_nu, ray.I[t], ray.ds[t], cell_alpha_nu,
                                    ray.T[t])
        s_tau += dtau
        if s_tau > break_tau:
            if ray.I[t+1] < 0:
                print("ERROR: Receiving a negative intensity, I[{:d}]={:.4e}"
                      .format(t+1, ray.I[t+1]))
            break
    return s_tau


def rad_step(nu, I0, ds, alpha_nu, T):
    """
    Description:
        Calculates the spectral irradiance of thermal emission in accordance
        with Planck's law.
        
    Notes:
        ^Assumes that variables are constant over distance ds

    Arguments:
        nu: frequency (in inverse seconds)
        I0: incoming spectral irradiance (in ergs per second
                                          per square centimeters per Hertz)
        ds: path length (in centimeters)
        alpha_nu: attenuation coefficient (in inverse centimeters)
        T: temperature (in Kelvin)
        
    Returns:
        Outgoing spectral irradiance & optical depth over distance ds
        
    Source:
        1. https://en.wikipedia.org/wiki/Radiative_transfer
    """
    # If variables are constant than optical depth reduces to alpha*ds
    #make earlier alpa_nu = n*sigma_j(nu) + ...
    #replace with dtau = (number dens * cross section) *ds - sum to get continuum
    dtau = alpha_nu*ds
    exp_dtau = np.exp(-dtau)
    # Attenuate intensity in accordance with attenuation coefficient,
    # and add black-body emission integral which again simplfies with
    # constants variables.
    # N.B. multipled by pi to turn radiance into irradiance
    I = I0*exp_dtau + const.pi*B_nu(nu, T)*(1.-exp_dtau)
    return I, dtau


def B_nu(nu, T):
    """
    Description:
        Calculates the spectral radiance of thermal emission in accordance
        with Planck's law.

    Arguments:
        nu: freuqnecy (in inverse seconds)
        T: temperature (in Kelvin)
        
    Returns:
        Spectral radiance due to thermal emission
        
    Source:
        1. https://en.wikipedia.org/wiki/Planck's_law
    """
    B = 2.*const.h*nu**3./const.c**2
    exponent = const.h*nu/(const.kB*T)
    if exponent > 500:
        return B/exponent
    else:
        return B/(np.exp(const.h*nu/(const.kB*T))-1.)


def shifted_nu(nu_0, v_LOS, obs_infront=True):
    """
    Description:
        Handles the sign convention of the velocity used for Doppler
        shifting the observed frequency. If the observer (the gas) is in
        front of the source (the star), the the line of sight velocity from
        the lab is negative of what the Doppler shifted velocity is defined
        as. If the gas is behind the source, then the line of sight velocity
        and the Doppler shifted velocity are the same thing.
        
    Arguments:
        nu_0: frequnecy in source frame (in inverse seconds)
        v_LOS: line of sight velocity (in centimeters per second)
        
    Keyword arguments:
        obs_infront: If the observers is in front of the emitter (boolean)
        
    Notes:
        ^v_LOS is positive if the gas is moving away from the lab
        ^v_obs is positive if the gas and star are receding from eachother
    
    Returns:
        Frequnecy in observer's (the gas') frame
    """
    if obs_infront:
        v_obs = -v_LOS
    else:
        v_obs = v_LOS
    return Doppler_shifted_nu(nu_0, v_obs)


def Doppler_shifted_nu(nu_0, v_obs):
    """
    Description:
        Returns the relativistic Doppler shift of the frequency between two
        frames that are moving at a velocity relative to one another.
        
    Arguments:
        nu_0: frequnecy in source frame (in inverse seconds)
        v_obs: velocity of observer in the emitter's frame (in centimeters
                                                            per second)
        
    Notes:
        ^v_obs is positive if the gas and star are receding from eachother

    Returns:
        Frequnecy in observer's frame
        
    Sources:
        1. https://en.wikipedia.org/wiki/Relativistic_Doppler_effect
    """
    beta_obs = v_obs/const.c
    return nu_0*np.sqrt((1.-beta_obs)/(1.+beta_obs))
