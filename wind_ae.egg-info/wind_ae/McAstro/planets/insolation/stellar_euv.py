#!/usr/bin/env python3

import numpy as np

from McAstro.utils import constants as const
from McAstro.stars.spectrum import high_energy
from McAstro.stars.relations import mass_relation
from McAstro.stars.relations import temperature_relation
from McAstro.stars.relations import colour_relation

class stellar_euv:
    def __init__(self, mass, semimajor, stellar_age=None, P_rot=None,
                 T_eff=None, radius=None, L_bol=None, Eker=True,
                 shift=True, verbose=True):
        """
        Description:
            Object for calculating the integrate uv flux of a star. For a
            given stellar mass and age (or rotation period) capable of
            estimating an euv flux at the given semimajor axis.

            As the euv of stars is hard to observe, models rely on using the
            x-ray luminosity. Each model has their limitions, and notably
            around the Kraft break the x-ray relations no longer hold.

            Units are currently not uniform, so care should be employed when
            using the methods.

        Arguments:
            mass: mass of star (in solar masses)
            semimajor: distance where flux is calculated (in centimeters)

        Keyword arguments:
            stellar_age: stellar age (in years)
            P_rot: rotation period (in seconds)
            T_eff: effective surface temperature (in Kelvin)
            radius: radius of star (in solar radii)
            L_bol: Bolometric luminosity (in solar luminosities)
            Eker: generate Eker star: radius, Lbol, Teff from mass (boolean)
            shift: shift Eker relations to be continuous (boolean)
            verbose: complains about extrapolations (boolean)

        Methods:
            calc_xray_luminosity()
            calc_euv_luminosity()
            integrate_uv()
        """
        self.mass = np.asarray(mass)
        self.semimajor = np.asarray(semimajor)
        # Warning flags
        self.verbose = verbose
        self.age_warning = False
        self.rot_warning = False
        self.stellar_age_estimted = False
        self.P_rot_estimated = False
        # Generate an Eker star given a mass
        if not Eker and (radius is None or T_eff is None or L_bol is None):
            print('ERROR: Must either generate Eker star or,'
                  '       provide stellar radius, effective temperature '
                  'and bolometric luminosity.')
            return
        if Eker:
            self.radius = mass_relation.Eker_mass_radius(self.mass, shift=shift)
            self.radius *= const.Rsun
            self.bolometric_lum = const.Lsun*(
                mass_relation.Eker_mass_luminosity(self.mass, shift=shift)
            )
            self.T_eff = mass_relation.Eker_mass_temperature(self.mass,
                                                             shift=shift)
        else:
            self.radius = radius*const.Rsun
            self.bolometric_lum = L_bol*const.Lsun
            self.T_eff = T_eff
        # Calculate BV colour based off Sung's results
        self.BV0 = temperature_relation.Sung_teff_to_BV(self.T_eff)
        # Attempt to set or calculate stellar ages and rotation rates
        if stellar_age is not None:
            self.stellar_age = np.asarray(stellar_age)
        else:
            self.stellar_age = None
        if P_rot is not None:
            self.P_rot = np.asarray(P_rot)
            if stellar_age is None:
                self.stellar_age_estimted = True
                self.age_warning, self.stellar_age = (
                    colour_relation.Mamajek_stellar_age(self.BV0, self.P_rot)
                )
        elif stellar_age is not None:
            self.P_rot_estimated = True
            self.rot_warning, self.P_rot = (
                colour_relation.Mamajek_rotation_rate(self.BV0,
                                                      self.stellar_age)
            )
        else:
            self.P_rot = None
        return


    def calc_xray_luminosity(self, xray_method):
        """
        Description:
            Calculates the x-ray luminosity. Some papers either were not
            forthright or I didn't spend enough time searching to find
            everyone's definition of x-ray band. If found universally used
            ROSAT's [0.1-2.4] keV band.

        Argument:
            xray_method: which model used to estimate an x-ray luminosity

        Models:
            Jackson: Uses stellar B-V color and age to estimate x-ray
            Owen: Uses Jackson's model with round numbers and solar fitted 
                  time dependent power-law index.
            SanzForcada:
            Tu:
            Wright: Uses stellar rotation rate to estimate the x-ray. As
                    stellar rotation rate is easier to observe than age this
                    a useful relation.
        """
        if xray_method == 'Jackson':
            if self.stellar_age is None:
                print('ERROR: Must specify a stellar age to use Jackson.')
            elif self.age_warning and self.verbose:
                print('WARNING: Age estimates are extrapolations.\n'
                      '         Results may not be valid, proceed w/caution.')
            # Caclulate x-ray luminosity using Jackson's result
            self.Jackson_warning, self.xray_lum = (
                high_energy.Jackson_xray_fraction(self.BV0, self.stellar_age)
            )
            self.xray_lum *= self.bolometric_lum
            if self.Jackson_warning and self.verbose:
                print('WARNING: Some B-V colours are outside of result range.\n'
                      '         Results interpolated, be wary of results.')
        elif xray_method == 'Owen':
            if self.stellar_age is None:
                print('ERROR: Must specify a stellar age to use Owen.')
            elif self.age_warning and self.verbose:
                print('WARNING: Age estimates are extrapolations.\n'
                      '         Results may not be valid, proceed w/caution.')
            self.xray_lum = 10**(-3.5)*self.bolometric_lum
            with np.errstate(all='ignore'): #ignore np.where internally complaints
                self.xray_lum *= np.where(self.stellar_age > 1e8,
                                          (self.stellar_age/1e8)**(-1.5), 1)
        elif xray_method == 'SanzForcada':
            if self.stellar_age is None:
                print('ERROR: Must specify a stellar age to use SanzForcada.')
            elif self.age_warning and self.verbose:
                print('WARNING: Age estimates are extrapolations.\n'
                      '         Results may not be valid, proceed w/caution.')
            self.xray_lum = (
                high_energy.SanzForcada_xray_luminosity(self.bolometric_lum,
                                                        self.stellar_age)
            )
        elif xray_method == 'Tu':
            if self.stellar_age is None:
                print('ERROR: Must specify a stellar age to use Tu.')
            elif self.age_warning and self.verbose:
                print('WARNING: Age estimates are extrapolations.\n'
                      '         Results may not be valid, proceed w/caution.')
            # Get rotation rate at 1 Megayear
            rot_rate_Myr = (2*const.pi)/(
                colour_relation.Mamajek_rotation_rate(self.BV0, 1e6)[1]
            )
            self.xray_lum = high_energy.Tu_xray_luminosity(rot_rate_Myr,
                                                           self.stellar_age)
        elif xray_method == 'Wright':
            if self.P_rot is None:
                print('ERROR: Must specify a stellar P_rot to use Wright.')
            elif self.rot_warning and self.verbose:
                print('WARNING: Rotation estimates are extrapolations.\n'
                      '         Results may not be valid, proceed w/caution.')
            self.xray_lum = high_energy.Wright_xray_fraction(self.mass,
                                                             self.P_rot
                                                             /const.day)
            self.xray_lum *= self.bolometric_lum
        else:
            print(f"ERROR: Xray method '{xray_method}' not recognized.\n"
                  "        Use: 'Jackson', 'Owen', 'SanzForcada', 'Tu', "
                  "'Wright'.")
            return -1
        self.xray_method = xray_method
        return


    def calc_euv_luminosity(self, euv_method):
        """
        Description:
            Calculates the euv luminosity. This is defined as the photons
            with energies between 100eV and RyH (hydrogen Rydberg constant).

        Argument:
            euv_method: which model used to estimate an euv luminosity

        Models:
            blackbody: Uses interpolated PHOENIX models
            Chadney: Uses solar observations to fit surface flux relation
            FISM2: Reproduction of Chadney (SEE observations) with FISM2
            SanzForcada:
        """
        # Calculate euv luminosity
        if euv_method == 'blackbody':
            lob = const.hc/(100*const.eV)
            upb = const.hc/const.RyH
            euv_surface_flux = high_energy.bb_euv_surface_flux(self.T_eff,
                                                               lob=lob, upb=upb)
            self.euv_lum = (4*const.pi*self.radius**2)*euv_surface_flux
        elif euv_method == 'Chadney' or euv_method == 'FISM2':
            if self.xray_lum is None:
                print('ERROR: Must first determine x-ray luminosity to use '
                      'Chadney.'
                      '       Use stellar_euv.calc_xray_luminosity(method).')
            xray_surface_flux = self.xray_lum/(4*const.pi*self.radius**2)
            if euv_method == 'Chadney':
                euv_surface_flux = (
                    high_energy.Chadney_euv_surface_flux(xray_surface_flux,
                                                         updated=False)
                )
            else:
                euv_surface_flux = (
                    high_energy.Chadney_euv_surface_flux(xray_surface_flux,
                                                         updated=True)
                )
            self.euv_lum = (4*const.pi*self.radius**2)*euv_surface_flux
        elif euv_method == 'SanzForcada':
            if self.xray_lum is None:
                print('ERROR: Must first determine x-ray luminosity to use '
                      'SanzForcada.'
                      '       Use stellar_euv.calc_xray_luminosity(method).')
            self.euv_lum = high_energy.SanzForcada_euv_luminosity(self.xray_lum)
        else:
            print(f"ERROR: uv method '{euv_method}' not recognized.\n"
                  "        Use: 'blackbody', 'Chadney', 'SanzForcada'")
            return -1
        self.euv_method = euv_method
        return


    def integrate_uv(self):
        """
        Description:
            Returns the total intgrated euv flux [100-13.6] eV at the
            location of the planet (in erg/s/cm**2 = mW/m**2).
        """
        return self.euv_lum/(4.*const.pi*self.semimajor**2)