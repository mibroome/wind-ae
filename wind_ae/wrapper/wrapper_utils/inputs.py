#!/bin/env python3

from . import constants as const
import pandas as pd
import numpy as np
import sys
import wind_ae.McAstro.atoms.atomic_species as McAtom
import importlib.resources as pkg_resources

class input_handler:
    """
    Class for handling input files (in self.path+'inputs/') for relaxation code
    """
    path = str(pkg_resources.files('wind_ae'))+'/'
    def __init__(self):
        self.path = str(pkg_resources.files('wind_ae'))+'/'
        return


    def write_planet_params(self, Mp, Ruv, Mstar, semimajor, Ftot, Lstar):
        '''
        Writes planet parameters to wind_ae/inputs/planet_params.inp input file.

        Args:
            Mp (float): Mass of the planet in g
            Ruv (float): Radius of the planet in cm
            Mstar (float): Mass of the star in g
            semimajor (float): Semimajor axis in cm
            Ftot (float): Total incident flux in erg/cm^2/s
            Lstar (float): Luminosity of the star in erg/s
        '''
        with open(self.path+'inputs/planet_params.inp', 'w') as f:
            f.write('#<parameter>: <value>           #<units>;     <comment>\n')
            f.write(f'Mp:        {Mp:.12e}   #g;          '
                    f'{Mp/const.Mjupiter:8.2f} MJ\n')
            f.write('Ruv:       %.12e   #cm;         %8.2f RJ\n'
                    % (Ruv, Ruv/const.Rjupiter))
            f.write('Mstar:     %.12e   #g;          %8.2f Msun\n'
                    % (Mstar, Mstar/const.Msun))
            f.write('semimajor: %.12e   #cm;         %8.2f au\n'
                    % (semimajor, semimajor/const.au))
            f.write('Ftot:      %.12e   #erg/cm^2/s; %8.2e FuvEarth\n'
                    % (Ftot, Ftot/const.FuvEarth))
            f.write('Lstar:     %.12e   #erg/s;      %8.2f Lstar\n'
                    % (Lstar, Lstar/const.Lsun))
            f.close()
            
            
    def write_additional_params(self, num_additional_params, additional_param_names=[0], additional_param_vals=[0]):
        '''
        Writes additional parameters to wind_ae/inputs/add_params.inp input file.
        This file is intended to allow users to add custom input parameters to the simulation without 
        having to modify the C source code as heavily.

        Args:
            num_additional_params (int): Number of additional parameters
            additional_param_names (list): List of additional parameter names
            additional_param_vals (list): List of additional parameter values
        '''
        with open(self.path+'inputs/add_params.inp', 'w') as f:
            f.write('#<parameter>: <value>           #<units>;     <comment>\n')
            f.write(f'Num_additional_params:        {num_additional_params:d}  ; \n')

            for i in range(num_additional_params):
                f.write('%s:       %.12e   #cgs units;\n'
                        % (additional_param_names[i], additional_param_vals[i])) #name of additional var + var itself
            f.close()



    def write_physics_params(self, mass_fraction, species_list, molec_adjust,
                             atomic_masses=np.array([0]),phys_file=path+'inputs/phys_params.inp'):
        '''
        Writes physics parameters to wind_ae/inputs/phys_params.inp input file.

        Args:
            mass_fraction (np.array): MASS fractions of each species
            species_list (list of str): List of species names (e.g., 'Fe III', 'fe4', etc.)
            molec_adjust (float): Molecular adjustment factor (default = 2.3 for molec. hydrogen)
            atomic_masses (np.array): Atomic masses of each species in g
            phys_file (str): Path to the physics parameters input file
        '''
        with open(phys_file, 'w') as f:
            f.write('#<parameter>: <value>         #<units>;         '
                    '<comment>\n')
            filepath = pkg_resources.files('wind_ae.McAstro.atoms').joinpath('atomic_table.txt')
            atomictable = pd.read_csv(filepath,comment='#')
            if len(mass_fraction) != len(species_list):
                sys.exit("ERROR: Mass fraction and species list must be the same length.")
            if np.round(np.sum(mass_fraction),5) != 1:
                print('WARNING: Total Mass Fraction must sum to 1. sum(ZX) = %.3f' %np.sum(mass_fraction))
            try:
                species_list = McAtom.formatting_species_list(species_list) #adds spaces and converts, e.g., fe6 to Fe IX
            except IndexError:
                print("ERROR: Species must include an ionization number either in Roman or Arabic numerals (e.g., H1, heIV, Fe 6)", 
                      file=sys.stderr) # Python 3.x
                sys.exit(1)
            Nspecies = len(mass_fraction)
            self.nspecies = Nspecies
            species_list_tight = [(species_list[i]).replace(' ','') for i in range(Nspecies)]
            
            HXstring = '%.12e' %mass_fraction[0]
            for i in range(1,Nspecies):
                HXstring+=', %.12e' %mass_fraction[i]
            Zstring=''
            Nestring=''
            massstring=''
            namestring=''
            for j in range(Nspecies):
                num1 = atomictable['Z'][atomictable['Name']==(species_list_tight[j])].iloc[0]
                num2 = atomictable['N_electrons'][atomictable['Name']==(species_list_tight[j])].iloc[0]
                if atomic_masses[0] == 0:
                    if species_list_tight[j] == 'HI':
                        num3 = 1.67330000000000013e-24
                    else:
                        num3 = atomictable['Atomic Mass'][atomictable['Name']==(species_list_tight[j])].iloc[0] #Verner masses
                else:
                    num3 = atomic_masses[j]
                name = species_list_tight[0]
                if j==0:
                    Zstring+='%d' %num1
                    Nestring+='%d' %num2
                    massstring+='%.12e' %1.67330000000000013e-24
                    namestring = '%s' %name
                else:
                    Zstring+=', %d' %num1
                    Nestring+=', %d' %num2
                    massstring+=', %.12e' %num3
                    namestring+=',%s' %(species_list_tight[j])
            f.write("HX:            "+HXstring+"       #;  Mass fractions\n")
            f.write("species_name:"+namestring+"#; Species name in Roman Numeral format\n") #must have no spaces for C readability
            f.write("atomic_mass:      "+massstring+"    #; Atomic mass in gs\n")
            f.write("Z:             "+Zstring+"                                                   #; Atomic number\n")
            f.write("Ne:            "+Nestring+"                                                   #; Number of electrons in species\n")
            f.write("Molec_adjust:            %.3f   #; Weighting factor to account for molecules below atomic wind" %molec_adjust)
            f.close()


    def write_spectrum(self, npts, nspecies, spec_date, spec_file, spec_kind, window,
                       resolved, normalized, ion_pot, E_wl, wPhi_wl, sigma_wl,
                       species):
        """
        Writes spectrum data to wind_ae/inputs/spectrum.inp input file. 
        NOTE: Windows are in NANOMETERS

        Args:
            npts (int): Number of points in the smoothed spectrum
            nspecies (int): Number of species
            spec_date (str): 'YYYY-MM-DD'; Date of the spectrum if using a Lisird spectrum ('0000-00-00' otherwise)
            spec_file (str): File name of spectrum if using custom (default='scaled-solar')
            spec_kind (str): Kind of the spectrum (options: 'mono', 'multi')
            window (list/array): Window for the spectrum IN NANOMETERS 
            resolved (list/array): Resolved values for the spectrum IN NANOMETERS (somewhat DEPRECATED, most functionalities can be covered by "window")
            normalized (list/array): Normalized values for the spectrum IN NANOMETERS(somewhat DEPRECATED, most functionalities can be covered by "window")
            ion_pot (list of floats): Ionization potentials for each species in ergs
            E_wl (array): Spectrum's photon energies (wavelengths) in ergs, array of length npts
            wPhi_wl (array): Number of photons per wavelength bin over Ftot, array of length npts
            sigma_wl (array): Ionization cross-sections for each species at E_wl, array [npts, nspecies]
            species (list of str): List of species names
        """
        try:
            if (len(E_wl) != npts or len(wPhi_wl) != npts or len(sigma_wl) != npts
                or len(ion_pot) != nspecies or len(sigma_wl[0]) != nspecies):
                print("ERROR: Spectrum write error, incorrect input sizes.")
                return
        except TypeError:
            E_wl = np.array([E_wl])
            wPhi_wl = np.array([wPhi_wl])
            sigma_wl = np.array([sigma_wl])
            
            if npts != 1:
                print("ERROR: Spectrum write error, incorrect input sizes. Npts should be 1")
                return
        with open(self.path+'inputs/spectrum.inp', 'w') as f:
            f.write(f'# NPTS: {npts:d}\n'
                    f'# NSPECIES: {nspecies:d}\n'
                    f'# DATE: {spec_date:s}\n'
                    f'# FILE: {spec_file:s}\n'
                    f'# KIND: {spec_kind:s}\n'
                    f'# WINDOW: {window[0]:.17e},{window[1]:.17e}\n'
                    f'# RESOLVED: {resolved[0]:.17e},{resolved[1]:.17e}\n'
                    f'# NORMALIZED: {normalized[0]:.17e},{normalized[1]:.17e}\n')
            f.write(f'# IONPOTS: {ion_pot[0]:.17e}')
            for s in range(1, nspecies):
                f.write(f',{ion_pot[s]:.17e}')
            f.write('\n'
                    r'# $hc/\lambda_i$, $w_i\Phi_{\lambda_i}/F_{tot}$')
            for sp in species:
                sp = McAtom.formatting_species_list([sp])[0].replace(' ','')
                f.write(r', $\sigma_{\lambda_i,'+sp+r'}$')
            f.write('\n')
            for b in range(npts):
                f.write(f'{E_wl[b]:.17e},{wPhi_wl[b]:.17e}')
                for s in range(nspecies):
                    f.write(f',{sigma_wl[b][s]:.17e}')
                f.write('\n')


    def write_bcs(self, Rmin, Rmax, rho_rmin, T_rmin, Ys_rmin, Ncol_sp, erf_drop):
        """
        Write boundary condition parameters to the wind_ae/inputs/bcs.inp file.

        Args:
            Rmin (float): Minimum radius in units of Rp
            Rmax (float): Maximum radius in units of Rp
            rho_rmin (float): Density at Rmin in units of RHO0 (default = 1.0E-15 g cm$^{-3}$)
            T_rmin (float): Temperature at Rmin in units of T0 (default = 1e4 K)
            Ys_rmin (list/array): Neutral fraction of each species at Rmin (default [1.0,...], where 1.0 is completely neutral)
            Ncol_sp (list/array): Column densities in units of cm$^{-2}$ at R_sp
            erf_drop (list/array): Drop-off parameters for the erf function that transitions from molecular to atomic and bolometric-heating-dominated to photoionization-heating-dominated regions (to compute, run sim.erf_velocity()).
                Ignored if bolo_heat_cool flag = 0 in inputs/flags.inp.
        """
        with open(self.path+'inputs/bcs.inp', 'w') as f:
            f.write('#<parameter>: <value>          #<units>; <comment>\n')
            f.write('Rmin:     %.12e   #R0;\n' % Rmin)
            f.write('Rmax:     %.12e   #R0;\n' % Rmax)
            f.write('rho_rmin: %.12e   #RHO0;\n' % rho_rmin)
            f.write('T_rmin:   %.12e   #T0;\n' % T_rmin)
            f.write('Ys_rmin:  %.12e' % Ys_rmin[0])
            for s in range(1,len(Ncol_sp)):
                f.write(',   %.12e' %Ys_rmin[s])
            f.write('#;\n')
            f.write('Ncol_sp:  %.12e' % Ncol_sp[0])
            for s in range(1,len(Ncol_sp)):
                f.write(',   %.12e' %Ncol_sp[s])
            f.write('#;\n')
            f.write('erf_drop:   %.12e,   %.12e#V0;' %(erf_drop[0],erf_drop[1]))


    def write_tech(self, breezeparam, rapidity, erfn, mach_limit):
        """
        Write technical parameters to the wind_ae/inputs/tech_params.inp file. Rarely needed.

        Args:
            breezeparam (float): Fraction of sound speed at critical point. If 1.0, is Mach 1 at the sonic point and solution is a transonic wind. If <1.0, is a subsonic breeze.
            rapidity (float): Rapidity of linearization transition at critical point.
            erfn (float): Used with rapidity to calculate analytic dvdr weighting (DISTINCT FROM bolometric/molecular complementary erf)
            mach_limit (float): Limit past which dvdr is purely linearized
        """
        if breezeparam > 1 or breezeparam < 0:
            print("breezeparam must be between 0 and 1")
            return 1
        if rapidity <= 0.:
            print("rapidity must be strictly positive")
            return 1
        if erfn <= 0:
            print("erfn must be strictly positive")
            return 1
        if mach_limit > 1 or mach_limit < 0:
            print("mach_limit must be between 0 and 1")
            return 1
        with open(self.path+'inputs/tech_params.inp', 'w') as f:
            f.write('#<parameter>: <value>                #<comment>\n')
            f.write('breezeparam:    %.12e   #Fraction of sound speed '
                    'at critical point\n'
                    % breezeparam)
            f.write("rapidity:       %.12e   "
                    "#Rapidity of linearization transition\n"
                    % rapidity)
            f.write("erfn:           %.12e   "
                    "#Closest non-unity weight is erf(erfn)\n"
                    % erfn)
            f.write("mach_limit:     %.12f       "
                    "#Limit past which dvdr is purely linearized\n"
                    % mach_limit)


    def write_flags(self, lyacool, tidalforce, bolo_heat_cool,
                    integrate_outward,integrate_out=True):
        """Write flags to the wind_ae/inputs/term_ind.inp file. 

        Args:
            lyacool (int): Flag for Lya and atomic line cooling (0 or 1)
            tidalforce (float): Tidal force scaling factor (1.0 = on, 0 = off)
            bolo_heat_cool (float): Bolometric heating/cooling scaling factor (1.0 = on, 0 = off). 
                                    If 0, complementary error function is not used and there is no
                                    bolometric heating/cooling included in the simulation and there is no
                                    molecular layer below wind (mu is purely atomic).
            integrate_outward (int): If 1, integrates out from sonic point to coriolis radius (0 or 1)
            integrate_out (bool): Whether to integrate outward (default: True) 
                                (implemented to avoid issues with integrate_outward flag being overwritten elsewhere)
        """
        with open(self.path+'inputs/term_ind.inp', 'w') as f:
            f.write('#<flag>: <boolean or scaling>  #<comment>\n')
            f.write('lyacool:            {:d}      #\n'.format(lyacool))
            f.write('tidalforce:         {:.5f}      #\n'.format(tidalforce))
            f.write('bolo_heat_cool:     {:.5f} #Can also be used to ramp in bolometric heating/cooling \n'.format(bolo_heat_cool))
            if integrate_out == False:    
                f.write('integrate_outward:  {:d}      #\n'.format(0))
            else:
                f.write('integrate_outward:  {:d}      #\n'.format(1))