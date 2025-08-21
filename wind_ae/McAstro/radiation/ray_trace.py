import numpy as np

from .radiative_transfer import radiative_transfer

def ray_trace(IP, ds, pop_ray, alpha_func, nu_ref, r_star=6.96e10, I0_star=1e0,
              savefile=None):
    """
    Description:
        Given an image plane and data structure, this function will
        calculate the optical depth and spectral irradiance after performing
        radiative transfer along the image plane's rays thru the data
        structure. Users must supply a function that generates rt_rays from
        the data structure, a function that calculates the attenuation
        coefficient, and the reference frequency observations are being made
        at.
        
    Arguments:
        IP: image plane object
        ds: data structure that contains simulation result
        pop_ray: a function that populations a rt_ray from the ds
        alpha_func: function which determines the attenuation coefficient
        nu_ref: referene frequency (in inverse seconds)
        
    Keyword arguments:
        r_star: radius of the star, outside of which I[0]=0 (in centimeters)
        I0_star: surface spectral irradiance of star (in ergs per second
                                                      per square centimeters
                                                      per Hertz)
        savefile: compressed numpy file used for saving ray tracing results
    
    Returns:
        void
    """
    if IP.ray_enter is None or IP.ray_exit is None:
        print("ERROR: Must calculate rays of image plane first.\n"
              "       Run image_plane.calc_rays().")
        return
    IP.tau, IP.I = (np.empty(IP.n_x1*IP.n_x2) for i in range(2))
    for j in range(IP.x2_min, IP.x2_max):
        for i in range(IP.x1_min, IP.x1_max):
            # Populate a rt_ray
            ray = pop_ray(i-IP.x1_min, j-IP.x2_min, ds, IP)
            # Set ray's intial flux
            if np.sqrt((j*IP.dx2)**2+(i*IP.dx1)**2) <= r_star:
                ray.I[0] = I0_star
            else:
                ray.I[0] = 0.0
            # Perform the radiative transfer
            index = (i-IP.x1_min)*IP.n_x1+(j-IP.x2_min)
            IP.tau[index] = radiative_transfer(nu_ref, ray, alpha_func)
            IP.I[index] = ray.I[-1]
    if savefile is not None:
        np.savez_compressed(savefile, x1=IP.x1, x2=IP.x2, tau=IP.tau, I=IP.I)
    return