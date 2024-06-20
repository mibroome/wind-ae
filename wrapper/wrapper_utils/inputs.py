#!/bin/env python3

from . import constants as const
import pandas as pd
import numpy as np
import sys
import McAstro.atoms.atomic_species as McAtom

class input_handler:
    """
    Class for handling input files (in 'inputs/') for relaxation code
    """
    def __init__(self):
        return


    def write_planet_params(self, Mp, Ruv, Mstar, semimajor, Ftot, Lstar):
        with open('inputs/planet_params.inp', 'w') as f:
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
        with open('inputs/add_params.inp', 'w') as f:
            f.write('#<parameter>: <value>           #<units>;     <comment>\n')
            f.write(f'Num_additional_params:        {num_additional_params:d}  ; \n')

            for i in range(num_additional_params):
                f.write('%s:       %.12e   #cgs units;\n'
                        % (additional_param_names[i], additional_param_vals[i])) #name of additional var + var itself
            f.close()



    def write_physics_params(self, mass_fraction, species_list, molec_adjust,
                             atomic_masses=np.array([0]),phys_file='inputs/phys_params.inp'):
        with open(phys_file, 'w') as f:
            f.write('#<parameter>: <value>         #<units>;         '
                    '<comment>\n')
            atomictable = pd.read_csv('McAstro/atoms/atomic_table.txt',comment='#')
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
                num1 = atomictable['Z'][atomictable['Name']==(species_list_tight[j])]
                num2 = atomictable['N_electrons'][atomictable['Name']==(species_list_tight[j])]
                if atomic_masses[0] == 0:
                    if species_list_tight[j] == 'HI':
                        num3 = 1.67330000000000013e-24
                    else:
                        num3 = atomictable['Atomic Mass'][atomictable['Name']==(species_list_tight[j])] #Verner masses
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
        with open('inputs/spectrum.inp', 'w') as f:
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
                    f'# $hc/\lambda_i$, $w_i\phi_{{\lambda_i}}$')
            for s in species:
                f.write(f', $\sigma_{{\lambda_i,{s}}}$')
            f.write('\n')
            for b in range(npts):
                f.write(f'{E_wl[b]:.17e},{wPhi_wl[b]:.17e}')
                for s in range(nspecies):
                    f.write(f',{sigma_wl[b][s]:.17e}')
                f.write('\n')


    def write_bcs(self, Rmin, Rmax, rho_rmin, T_rmin, Ys_rmin, Ncol_sp, erf_drop):
        with open('inputs/bcs.inp', 'w') as f:
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
        with open('inputs/tech_params.inp', 'w') as f:
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
        with open('inputs/term_ind.inp', 'w') as f:
            f.write('#<flag>: <boolean>  #<comment>\n')
            f.write('lyacool:            {:d}      #\n'.format(lyacool))
            f.write('tidalforce:         {:d}      #\n'.format(tidalforce))
            f.write('bolo_heat_cool:     {:.5f} #Can also be used to ramp in bolometric heating/cooling \n'.format(bolo_heat_cool))
            if integrate_out == False:    
                f.write('integrate_outward:  {:d}      #\n'.format(0))
            else:
                f.write('integrate_outward:  {:d}      #\n'.format(1))
                