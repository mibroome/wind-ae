#!/usr/bin/env python3

import numpy as np

from wind_ae.McAstro.utils import constants as const
from wind_ae.McAstro.stars.spectrum import high_energy
from wind_ae.McAstro.stars.relations import mass_relation
from wind_ae.McAstro.stars.relations import temperature_relation
from wind_ae.McAstro.stars.relations import colour_relation

def integrated_uv(mass, semimajor, stellar_age=1e9, P_rot=None,
                  xray='Jackson', euv='Chadney', updated=False,
                  verbose=False, extrapolate=False, shift=False):
    """
    Description:
        For a given stellar mass and age calculates the euv flux at the
        given semimajor axis. Several methods avaiable for calculating
        the euv luminoisty of a star, but currently all first calculate
        the x-ray luminosity.
    
    Arguments:
        mass: mass of star (in solar masses)
        semimajor: distance where flux is calculated (in centimeters)
        
    Keyword arguments:
        stellar_age: age of star (in years)
        xray: method for determining x-ray luminosity
            ('Jackson', 'Owen', 'SanzForcada', 'Tu', 'Wright')
        euv: method for determining euv luminosity
            ('blackbody', 'Chadney', 'SanzForcada')
        verbose: print debugging info (boolean)
        extrapolate: if methods extrapolate outside their fitted data
        
    Returns:
        Integrated uv flux ([13.6-100] eV) at orbital distance of planet
        (in erg/s/cm**2 = mW/m**2)
    """
    if P_rot is not None and xray != 'Wright':
        print("WARNING: only xray='Wright' uses P_rot.\n"
              "         To use P_rot retry w/keyword argument: xray='Wright'.")
    mass = np.asarray(mass)
    semimajor = np.asarray(semimajor)
    stellar_age = np.asarray(stellar_age)
    P_rot = np.asarray(P_rot)
    # Generate an Eker star given a mass
    radius = mass_relation.Eker_mass_radius(mass, shift=shift)*const.Rsun
    bolometric_lum = mass_relation.Eker_mass_luminosity(mass, shift=shift)
    bolometric_lum *= const.Lsun
    T_eff = mass_relation.Eker_mass_temperature(mass, shift=shift)
    # Calculate x-ray luminosity
    if xray == 'Jackson':
        # Caclulate x-ray luminosity using Jackson's result
        BV0 = temperature_relation.Ballesteros_teff_to_BV(T_eff)
        xray_lum = high_energy.Jackson_xray_fraction(BV0, stellar_age,
                                                     extrapolate=extrapolate)
        xray_lum *= bolometric_lum
    elif xray == 'Owen':
        # Uses Jackson's model with rough round numbers that work across colour
        # bins and solar fitted time dependent power-law index
        xray_lum = 10**(-3.5)*bolometric_lum
        if stellar_age > 1e8:
            xray_lum *= (stellar_age/1e8)**(-1.5)
    elif xray == 'SanzForcada':
        xray_lum = (
            high_energy.SanzForcada_xray_luminosity(bolometric_lum, stellar_age)
        )
    elif xray == 'Tu':
        # Get rotation rate at 1 Megayear
        BV0 = temperature_relation.Ballesteros_teff_to_BV(T_eff)
        rot_period_Myr = (
            colour_relation.Mamajek_rotation_rate(BV0, 1e6,
                                                  extrapolate=extrapolate)
        )
        rot_rate_Myr = (2*const.pi)/(rot_period_Myr)
        xray_lum = high_energy.Tu_xray_luminosity(rot_rate_Myr, stellar_age)
    elif xray == 'Wright':
        xray_lum = high_energy.Wright_xray_fraction(mass, P_rot)
        xray_lum *= bolometric_lum
    else:
        print(f"ERROR: Xray method '{xray}' not recognized.\n"
              "        Use: 'Jackson', 'Owen', 'SanzForcada', 'Tu', 'Wright'.")
        return -1
    # Calculate euv luminosity
    if euv == 'blackbody':
        lob = const.hc/(100*const.eV)
        upb = const.hc/const.RyH
        euv_surface_flux = high_energy.bb_euv_surface_flux(T_eff, lob=lob,
                                                           upb=upb)
        euv_lum = (4*const.pi*radius**2)*euv_surface_flux
    elif euv == 'Chadney':
        xray_surface_flux = xray_lum/(4*const.pi*radius**2)
        euv_surface_flux = (
            high_energy.Chadney_euv_surface_flux(xray_surface_flux,
                                                 updated=updated)
        )
        euv_lum = (4*const.pi*radius**2)*euv_surface_flux
    elif euv == 'SanzForcada':
        euv_lum = high_energy.SanzForcada_euv_luminosity(xray_lum)
    else:
        print(f"ERROR: uv method '{euv}' not recognized.\n"
              "        Use: 'blackbody', 'Chadney', 'SanzForcada'")
        return -1
    # Check if user wants debugging print out
    if verbose:
        print('Rad: {:e} Rsun\nBol: {:e} Lsun\nTemp: {:e} K\nB–V: {:e}\n'
              'L_x: {:e} Lsun\nL_euv: {:e}'
              .format(radius/const.Rsun, bolometric_lum/const.Lsun, T_eff, BV0,
                      xray_lum/const.Lsun, euv_lum/const.Lsun))
    # return euv flux at the planet's location
    return euv_lum/(4*const.pi*semimajor**2)