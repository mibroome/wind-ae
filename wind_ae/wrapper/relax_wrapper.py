#!/usr/bin/env python3
"""
relax_wrapper.py:
    Contains the class for wrapping the relaxation code. Handles loading, saving, running, and converging boundary conditions for the relaxation simulations. 
    The goal is find the solution of interest to the end-user and automate the choices one must make to be intelligently picked and consistent, i.e., the boundary conditions and ramping thru parameter space.
    Handling of the solution is done by the other major class wind_solution found in wind_ae/wrapper/wrapper_utils/windsoln.py.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE, STDOUT
from scipy.interpolate import CubicSpline 
import scipy.interpolate as si 
from scipy.optimize import fsolve
from os.path import exists
import math
import sys
from IPython import display
from scipy.interpolate import CubicSpline
from scipy.special import exp1

from wind_ae.wrapper.wrapper_utils import constants as const
from wind_ae.wrapper.wrapper_utils.system import system
from wind_ae.wrapper.wrapper_utils.physics import physics
from wind_ae.wrapper.wrapper_utils.inputs import input_handler
from wind_ae.wrapper.wrapper_utils.windsoln import wind_solution
from wind_ae.wrapper.wrapper_utils.plots import four_panel_plot, six_panel_plot, energy_plot
from wind_ae.wrapper.wrapper_utils.plots import quick_plot
from wind_ae.wrapper.wrapper_utils.metals import metal_class
from wind_ae.wrapper.wrapper_utils.spectrum import spectrum
import wind_ae.McAstro.atoms.atomic_species as McAtom
import pandas as pd
import time
import importlib.resources as pkg_resources


class wind_simulation:
    path = str(pkg_resources.files('wind_ae'))+'/'
    def __init__(self, csv_file=path+'inputs/guess.inp', name='Init. Planet',
                 expedite=True):
        self.path = str(pkg_resources.files('wind_ae'))+'/'
        self.last_print_rastered = False
        self.strat_kappa = 1e-2
        self.inputs = input_handler()
        self.first_print = True
        self.skip = False
        self.static_bcs = False
        self.df = pd.read_csv(pkg_resources.files('wind_ae.wrapper.wrapper_utils').joinpath('dere_table29.dat'),
                              sep=r'\s+', names=list(range(45)))  
        self.df = self.df.rename(columns={0:'Z', 1:'Ion', 2:'NS', 3:'I', 4:'Tmin'})
        self.clear = 100*' '+'\n'
        self.try_turning_off = True
#         self.width_factor = 1

    def _raster_print(self,msg,end='', pad=130):
        print('\r'+msg + ' '*(pad - len(msg)), end=end)
        self.last_print_rastered = True
        return
    def _normal_print(self,msg):
        if self.last_print_rastered:
            print('\n'+msg+self.clear)
        else:
            print(msg+self.clear)
        self.last_print_rastered = False
        return
        
    # def load_nonexpedited(self, csv_file=None):
    #     if csv_file is None:
    #         csv_file = self.last_load
    #     self.windsoln = wind_solution(file=csv_file,calc_postfacto=True)
    #     return


    def load_uservars(self, csv_file=None):
        """
        Loads all user variables for plotting and analysis. Does not rewrite input parameter files, so can be run concurrently with simulations.

        Args:
            csv_file (str, optional): Path to the csv file containing the planet solution. Defaults to None.

        Returns:
            None

        """
        if csv_file is None:
            csv_file = self.last_load
        self.windsoln = wind_solution(file=csv_file, calc_postfacto=True)
        return
    

    def generate_rate_coeffs(self):
        """
        Generates rate coefficients for secondary ionization and populates array in src/rate_coeffs.h.
        Rate coefficients are generated from interpolation over the Dere (2007) table.

        Returns:
            None

        """
        new_spec = np.genfromtxt(self.path+'inputs/spectrum.inp',skip_header=9,delimiter=',')
        E_wl = new_spec[:,0]
        f = open(self.path+'/inputs/spectrum.inp','r')
        ff = f.readlines()
        nspecies = int(ff[1].split(':')[1])
        species_new = [sp[:-2] for sp in ff[9].split(',')[2:][1::2]]
        species_new[-1] = species_new[-1][:-1]
        f.close()
        def Interpolater(j,m):
            species_ionized = species_new[m]
            species_name_spaced = McAtom.formatting_species_list([species_ionized])[0]
            ma = McAtom.atomic_species(species_name_spaced)
            Z = ma.Z
            Ne = ma.Ne
            species_ionizing = species_new[j]
            species_name_spaced = McAtom.formatting_species_list([species_ionizing])[0]
            ma = McAtom.atomic_species(species_name_spaced)
            ion_pot_m = ma.ion_pot

            Ion = Z-Ne+1
            dere = self.df[(self.df['Z']==Z)&(self.df['Ion']==Ion)].values[0][5:]
            dere = dere[np.isnan(dere)==False]
            split = len(dere)//2

            x = dere[:split]
            rho = dere[split:]
            I = self.df[(self.df['Z']==Z)&(self.df['Ion']==Ion)]['I'].values[0]
            k = 8.67e-5  # eV/K
            f = 2
            T = (I/k)*(np.exp(np.log(f)/(1-x))-f)
            E = k*T
            t = (k*T)/I
            R = t**(-1/2)*I**(-3/2)*rho*exp1(1/t)
            # spline = interp.CubicSpline(E[(E>40)&(E<ulim)]*const.eV,R[(E>40)&(E<ulim)])
            spline = CubicSpline(E*const.eV,R)

            E0 = E_wl - ion_pot_m #energy left after ionizing species m
            return spline(E0)

        # with open(self.path+'inputs/add_params.inp', 'w') as f:
        f = open(self.path+'src/rate_coeffs.h','w')
        f.write('/*Array of rate coefficients (cm3/s) from Dere (2007)*/\n')
        f.write('/*of species m being secondarily ionized (Col 1) by a photoelectron released after impacting species j (Col 2)*/\n')
        f.write('/*as a function of E0 = E - I[j] (ergs) from the spectrum.*/\n')

        nspecies = len(species_new)
        f.write("double R[%d][%d] = {"%(nspecies**2+nspecies,len(E_wl)))
        # for species in sim.windsoln.species:
        for j in range(nspecies):
            for m in range(nspecies):
                # f.write("{%s,%s," %(sim.windsoln.species[j],sim.windsoln.species[m]))
                Rs = Interpolater(j,m,)
                f.write("{"+(','.join('%.4e' %rs for rs in Rs))+"},\n")
        f.write('};\n')
        f.close()
        return
    
    def load_planet(self, csv_file, calc_postfacto=True, name='Loaded Planet',
                    print_atmo=True, print_warnings=True):
        """
        Loads a planet solution file into inputs/guess.inp as a new guess and updates input parameter folders.
        Ramping to a new planet with ramp_to() will rewrite parameter folders for the desired new planet, but the guess will remain the same.

        Args:
            csv_file (str): Path to the csv file containing the planet solution.
            calc_postfacto (bool, optional): If True, calculates post-facto variables (e.g., photoionization heating, etc.). Defaults to True.
            name (str, optional): Name of the planet, used for identification. Defaults to 'Loaded Planet'.
            print_atmo (bool, optional): If True, prints the atmosphere composition and mass fractions. Defaults to True.
            print_warnings (bool, optional): If True, prints warnings from the wind_solution analysis. Defaults to True.

        Returns:
            None
        """
        # load header and data
#         self.failed_deeper_bcs_ramp = False
        planet = wind_solution(file=csv_file, 
                               calc_postfacto=calc_postfacto,print_warnings=print_warnings)
        for j in range(planet.nspecies): #had a problem with adding spaces to e.g., 'He I'
            planet.species_list[j] = (planet.species_list[j]).replace(' ','')
            
        #Changing Nspecies, M (# pts in relaxation region), RHOSCALE (convergence condition) and remaking C code    
        f = open(self.path+'src/defs.h', 'r') 
        for line in f.readlines():
            splitline = line.split()
            if len(splitline) >= 2:
                line_var = splitline[0]+' '+splitline[1]
            if line_var == '#define NSPECIES':
                nspecies_def = int(line.split()[2])
            elif line_var == '#define M':
                og_length = int(line.split()[2]) 
            elif line_var == '#define RHOSCALE':
                og_RHOSCALE = float(line.split()[2])
            elif line_var == '#define N_ADD_PARAMS':
                og_N_ADD_PARAMS = float(line.split()[2])
        f.close()
        
        #rewriting
        new_length = len(planet.soln['q'][planet.soln['q']<=1])
        new_RHOSCALE = 10**np.floor(np.log10(planet.rho_rmin*0.01))
        nspecies_new = int(planet.nspecies)
        new_N_ADD_PARAMS = planet.n_add_prms
        
        f = open(self.path+'src/defs.h','w')
        h = (open(self.path+'src/defs-master.h','r')).readlines()
        for idx,hline in enumerate(h):
            splitline = hline.split()
            if len(splitline) >= 2:
                line_var = splitline[0]+' '+splitline[1]
            if line_var == '#define NSPECIES':
                index1 = idx
            elif line_var == '#define M':
                index2 = idx
            elif line_var == '#define RHOSCALE':
                index3 = idx
            elif line_var == '#define N_ADD_PARAMS':
                index4 = idx
        h[index1] = '#define NSPECIES %d\n' %nspecies_new
        h[index2] = '#define M %d            /* number of points */\n' %new_length
        h[index3] = '#define RHOSCALE %.1f\n' %new_RHOSCALE
        h[index4] = '#define N_ADD_PARAMS %d\n' %new_N_ADD_PARAMS
        f.writelines(h)
        f.close() 

        remake = False
        if nspecies_def != nspecies_new:
            self._raster_print(f'Nspecies has changed from %d to %d. Remaking C code...' %(nspecies_def, nspecies_new))
            remake = True
        #For regridded solutions, the number of points in relaxtion region may change
        if og_length != new_length:
            self._normal_print("\nNumber of points in relaxation region has changed from %d to %d. Remaking C code..." 
                  %(og_length, new_length))
            remake=True
        #For solutions with very high rho at the lower boundary, the rho convergence condition should be raised
        if new_RHOSCALE != og_RHOSCALE:
            self._normal_print("\nRHOSCALE (convergence condition) has changed from %d to %d. Remaking C code..." 
                  %(og_RHOSCALE, new_RHOSCALE)) 
            remake=True
        #Future users can add additional parameters
        if new_N_ADD_PARAMS != og_N_ADD_PARAMS:
            self._normal_print("N_ADD_PARAMS (num. of additional params added by user) has changed from %d to %d. Remaking C code..." 
                  %(og_N_ADD_PARAMS, new_N_ADD_PARAMS)) 
            remake=True
        if remake == True:
            sub = Popen('make',cwd=self.path, stdout=PIPE, stderr=PIPE) 
            output, error_output = sub.communicate()
        
        # write input parameters
        self.inputs.write_planet_params(*planet.planet_tuple)
        self.inputs.write_physics_params(*planet.physics_tuple)
        self.inputs.write_bcs(*planet.bcs_tuple)
        self.inputs.write_tech(*planet.tech_tuple)
        #should load with solution in guess
        self.inputs.write_flags(*planet.flags_tuple,integrate_out=planet.flags_tuple[3])
        self.inputs.write_additional_params(*planet.add_tuple)
        self.inputs.write_spectrum(*planet.spectrum_tuple)
        self.generate_rate_coeffs()
        # load planet as guess.inp if not already
        if csv_file != 'inputs/guess.inp':
            sub = Popen(["cp", csv_file, self.path+'inputs/guess.inp'],
                        stdout=PIPE, stderr=PIPE)
            sub.wait()
            output, error_output = sub.communicate()
            if error_output:
                self._normal_print(f'ERROR: Failed loading {csv_file:s} as guess.inp\n')
                self._normal_print(error_output)
                return
        # If all is successful update wind_simulation object
        sub = Popen(["cp", 'inputs/guess.inp', 'saves/windsoln.csv'],cwd=self.path, stdout=PIPE, stderr=PIPE)
        sub.wait()
        output, error_output = sub.communicate()
        if error_output:
            self._normal_print(f'ERROR: Failed copying {csv_file:s} into saves/windsoln.csv\n')
            self._normal_print(error_output)
            return
        
        self.guess = planet
        self.windsoln = planet
        self.system = system(*planet.planet_tuple, name=name)
        self.mu = planet.calc_mu()[0]
        # self.atmosphere = atmosphere(self.system,
        #                              self.guess.T_rmin*self.guess.scales[3],
        #                              self.mu,
        #                              kappa_opt=self.strat_kappa)
        # self.atmosphere.windbase_radius(self.guess.rho_rmin
        #                                 *self.guess.scales[1])
        self.physics = physics(*planet.physics_tuple)
        if print_atmo==True:
            print("Atmosphere Composition")
            print('  Species:   '+(',        '.join('%s' %sp.replace(' ','') for sp in planet.species_list)))
            print('  Mass frac: '+(', '.join('%.2e' %hx for hx in planet.HX)))
            print('')
        
        self.ramp_class = "system"
        self.last_load = csv_file
        return


    def load_spectrum(self, generate=False, wl_norm=1e-7, print_warnings=False):
        """
        Loads a high resolution spectrum and, if generate=True, smooths it, keeping the flux at ionization edges for the species in the current windsoln accurate, and conserving flux locally in the vicinity of spectral peaks. The spectrum specified in the windsoln is the one loaded and smoothed.
        One can change the spectrum using ramp_to_user_spectrum() to ramp to a user-input spectrum or change the wavelength range and integrated flux via ramp_spectrum(). To change the total flux but maintain spectrum shape and range, use ramp_flux().

        Args:
            generate (bool, optional): If True, generates a new smoothed spectrum stored in inputs/spectrum.inp. Defaults to False.
            wl_norm (float, optional): Converts nm to cm. Should not need to be changed. Defaults to 1e-7.
            print_warnings (bool, optional): If True, prints warnings from McAstro/planets/insolation/glq_spectrum.csv. Defaults to False.

        Returns:
            None
        """
        if self.windsoln.spec_src_file != 'scaled-solar':
            self.spectrum = spectrum(lisird=False,spectrum_file=self.windsoln.spec_src_file,
                                     wl_norm=wl_norm,print_warning=print_warnings)
        else:
            self.spectrum = spectrum(date=self.windsoln.spec_date, wl_norm=wl_norm,
                                     print_warning=print_warnings)
        for name in self.windsoln.species_list:
            self.spectrum.add_species(name)
        soln_resolved = self.windsoln.spec_resolved/wl_norm
        self.spectrum.set_resolved(*soln_resolved)
        soln_normalized = self.windsoln.spec_normalized/wl_norm
        self.spectrum.set_normalized(*soln_normalized)
        soln_window = self.windsoln.spec_window/wl_norm
        self.spectrum.set_window(*soln_window, kind=self.windsoln.spec_kind)
        for j in range(self.windsoln.nspecies):
            self.windsoln.species_list[j] = (self.windsoln.species_list[j]).replace(' ','')
        if generate:
            self.spectrum.generate(kind=self.windsoln.spec_kind, savefile=self.path+'inputs/spectrum.inp')
            self.generate_rate_coeffs()
        return

    
    def save_planet(self, filename=None, filepath='', overwrite=False, polish=False):
        """
        Copies wind_ae/saves/windsoln.csv to the desired folder.
        Note: It is not necessary to define a folder. All files will be saved in the saves/ repo by default.

        Args:
            filename (str, optional): Name and/or path+name. .csv not required. Defaults to None.
            filepath (str, optional): Path to the folder where the file will be saved. Defaults to ''.
            overwrite (bool, optional): If True, will overwrite existing file with the same name. Defaults to False.
            polish (bool, optional): If True, will polish the boundary conditions before saving. Defaults to False.

        Returns:
            None
        """
        if polish == True:
            self.polish_bcs()
        folder = filepath
        name = filename
        if folder != '':
            if len(folder) > 1 and folder[-1] != '/':
                folder += '/'
            os.makedirs(folder, exist_ok=True)
        if name is not None:
            if name[-4:] != '.csv':
                name = folder+name+'.csv'
            else:
                name = folder+name
        else:
            Mp, Rp, Mstar, semimajor, Ftot, Lstar = self.windsoln.planet_tuple
            name = (self.path+("saves/Mp:{:g}_Rp:{:g}_Mstar:{:g}_a:{:g}_Ftot:{:g}_Nsp:{:d}_Spect:{:s}.csv"
                             .format(Mp, Rp, Mstar, semimajor, Ftot, 
                                     self.windsoln.nspecies, self.windsoln.spectrum_tuple[3])))
        if not overwrite and os.path.isfile(name):
            self._normal_print('File already exists.\n'
                  '  To overwrite use save_planet(overwrite=True).')
            return
        print("Saving %s" % name)
        sub = Popen(["cp", self.path+'saves/windsoln.csv', name],
                    stdout=PIPE, stderr=PIPE)
        sub.wait()
        output, error_output = sub.communicate()
        if error_output:
            self._normal_print(error_output)
        return
    
    
    
    def easy_output_file(self, outputs=['v', 'T', 'rho'],
                         output_file=path+'saves/simplified_outputs/output.dat',
                         comments='', overwrite=False):
        """
        Writes desired solution variables to a csv file that can be easily shared.

        Args:
            outputs (list of str, optional): Desired solution columns. 'r' not required. Defaults to ['v', 'T', 'rho'].
            output_file (str, optional): Output file name. Can include path. Defaults to 'output.dat'.
            comments (str, optional): Additional comments to include in the header. Defaults to ''.
            overwrite (bool, optional): If True, will overwrite existing file. Defaults to False.

        Returns:
            None
        """
        if not overwrite and os.path.isfile(output_file):
            self._normal_print('File already exists.\n'
                  '  To overwrite use set overwrite=True.')
            return
        f = open(output_file,'w')
        self.windsoln.add_user_vars()
        s = self.windsoln
    #     if len(comments) == 0:
        header = f'#Planet: {s.Mp/const.Mjupiter:.4f} Mjup, {s.Rp/const.Rjupiter:.4f} Rjup, {s.semimajor/const.au:.3f} au \n'
        header += f'#Star: {s.Mstar/const.Msun:.3f} Msun, {s.Lstar/const.Lsun:.3f} Lsun \n'
        header += f'#Spectrum: {s.Ftot:.0f} ergs/s/cm2, {s.spectrum_tuple[3]}, [{s.spec_normalized[0]/1e-7:.3f},{s.spec_normalized[1]/1e-7:.3f}] nm\n'
        comments = header + '#'+comments+'\n'    
        f.write(comments)
        col_line='r,r/Rp'
        for var in outputs:
            if var != 'r':
                col_line += f',{var}'
        col_line+='\n'
        f.write(col_line)

        for i in range(len(s.soln['r'])):
            line = '%.1f,%.3f' %(s.soln['r'][i],s.soln_norm['r'][i])    
            for var in outputs:
                if var != 'r':
                    line += ',%.5e' %s.soln[var][i]
            line+='\n'
            f.write(line)
        f.close() 
        
        return


    def run_wind(self, expedite=False, calc_postfacto=True, verbose=False, retry=False):
        """
        Runs the C relaxation code. If successful, returns 0, saves solution in saves/windsoln.csv, and copies this solution to be the next guess in inputs/guess.inp.

        Args:
            expedite (bool, optional): If True, will also run outward integration. Will return 4 if relaxation code has solved, but outward integration fails (does not affect profile or mass loss accuracy). Defaults to False.
            calc_postfacto (bool, optional): If True, performs post-facto calculations (e.g., heating/cooling, number density, etc.). Defaults to True.
            verbose (bool, optional): If True, prints verbose output. Defaults to False.
            retry (bool, optional): If True, retries on analysis errors. Defaults to False.

        Returns:
            int: Status code.
                0 - Both relaxation and outward integration have succeeded (if expedite=False), or relaxation has succeeded (if expedite=True).
                4 - Relaxation has succeeded, but outward integration has failed.
                1 - Relaxation has failed.
                2 - Error in copying saves/windsoln.csv to inputs/guess.inp (rare).
                3 - Analysis error when executing wind_solution() to load current solution into windsoln object.
        """
        # Run wind relaxation code
        if expedite == True: #do not integrate out, else use default of loaded solution
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)
        sub = Popen('./bin/relaxed_ae',cwd=self.path, stdout=PIPE, stderr=PIPE) 
        output, error_output = sub.communicate()
        if error_output:
            if error_output[0:55] == b'Numerical Recipes run-time error...\nstep size underflow':
                return 4
            else:
                return 1
            if verbose:
                self._normal_print(error_output)
        else:
            # If successful then update our guess to new solution
            sub = Popen(["cp", 'saves/windsoln.csv', 'inputs/guess.inp'],cwd=self.path,
                        stdout=PIPE, stderr=PIPE)
            output, error_output = sub.communicate()
            if error_output:
                self._normal_print(error_output)
                return 2
            #if calc_postfacto==True, takes time to run add_user_vars()
            self.windsoln = wind_solution(calc_postfacto=calc_postfacto, print_warnings=False)
            if self.windsoln.error:
                # If we ran into analysis errors and not runtime errors, then
                # we probaby need to converge the bcs.
                if not retry:
                    # Check retry to avoid infinite recursion
                    self._normal_print("Retry triggered")
                    self.converge_mol_atomic_transition(_called_in_ramp_bcs=True)
                    self.converge_Ncol_sp(expedite=False)
                    return self.run_wind(retry=True)
                else:
                    self._normal_print("Failed on retry to fix analysis errors")
                    return 3
            # saves/windsoln.csv and inputs/guess.inp are the same
            self.guess = self.windsoln
            return 0
                     

# Ramping Functions
    def ramp_to(self, system=None, intermediate_converge_bcs=False, final_polish=False,
                integrate_out=True, static_bcs=False, make_plot=False):
        """
        Ramps the current wind solution to a new system. Wrapper for ramp_Ftot(), ramp_grav(), ramp_star().

        Args:
            system (system, optional): The target system to ramp towards. Syntax: system = system(Mp, Rp, Mstar, a, Ftot, Lstar) in cgs. Defaults to None.
            intermediate_converge_bcs (bool, optional): If True, converges boundary conditions during and after each ramp function. Increases runtime, but may improve likelihood of convergence for sensitive solutions or solutions distant from the initial guess. Defaults to False.
            final_polish (bool, optional): If True, converges boundary conditions (i.e. polish_bcs()) after all ramping is done. Defaults to False.
            integrate_out (bool, optional): If True, integrates final solution out past sonic point to Coriolis radius. Defaults to True.
            static_bcs (bool, optional): If True, keeps current boundary conditions static during ramping. Defaults to False.
            make_plot (bool, optional): If True, generates velocity, temperature, density, and ionization fraction plots of steps in the ramping process. Defaults to False.

        Returns:
            int: Status code. 0 for success, other values for failure modes.
        """
        if (system is None):# and physics is None):
            self._normal_print("Please provide a system to ramp towards. system=system(Mp,Rp,Mstar,a,Ftot,Lstar).")
            return 0
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        fail = 0
        result = 0
        if system is not None:
            self.ramp_class = "system"
            result = self.ramp_var("Ftot", system.value("Ftot"),
                                  converge_bcs=intermediate_converge_bcs, make_plot=make_plot,
                                  expedite=True,static_bcs=static_bcs,
                                  integrate_out=False)
            fail = result
            
            if (fail != 0) and (fail != 5):
                return fail  
            
            # Ramps Mp and Rp simultaneously, ~constant surface gravity
            result = self.ramp_grav(system, converge_bcs=intermediate_converge_bcs,
                                   make_plot=make_plot, expedite=True,static_bcs=static_bcs,
                                   integrate_out=False)
            fail = result 
            
            if (fail != 0) and (fail != 5):
                return fail            
                
            # Ramps Mstar and semimajor simultaneously, ~constant Hill radius
            result = self.ramp_star(system, converge_bcs=intermediate_converge_bcs,
                                   make_plot=make_plot, expedite=True,static_bcs=static_bcs,
                                   integrate_out=False)
            fail = result

        if final_polish == True:
            self.polish_bcs(integrate_out,static_bcs)            
        elif (integrate_out == True) & (final_polish == False):
            self.converge_Rmax()

        return fail

    
    def ramp_Ftot(self, F_goal, converge_bcs=False, expedite=False, integrate_out=True,
                  static_bcs=False, make_plot=True):
        """
        Ramps the current wind solution to a new Ftot in ergs/s/cm2.

        Args:
            F_goal (float): The target Ftot value to ramp towards.
            converge_bcs (bool, optional): If True, converges boundary conditions after ramping Ftot. Defaults to False.
            expedite (bool, optional): If True, expedites the ramping process, but may reduce likelihood of convergence for sensitive solutions. Defaults to False.
            integrate_out (bool, optional): If True, integrates out past the sonic point. Defaults to True.
            static_bcs (bool, optional): If True, uses static boundary conditions. Defaults to False.
            make_plot (bool, optional): If True, produces a four panel plot of velocity, temperature, density, and ionization fraction for steps in the ramping process. Defaults to True.

        Returns:
            int: Status code from ramp_var.
        """
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        return self.ramp_var("Ftot",F_goal,converge_bcs=converge_bcs,make_plot=make_plot,
                        expedite=expedite,integrate_out=integrate_out,static_bcs=static_bcs)
    
    
    def ramp_var(self, var, var_end, var_class=None, delta=0.02,
                 delta_additive=False, converge_bcs=False, make_plot=True,
                 expedite=False, integrate_out=True, static_bcs=False):
        """
        Ramps a variable in the system or physics class to a target value with adaptive stepsizes.

        Args:
            var (str): Name of the variable to ramp. Options:
                var_class='system': var = 'Ftot', 'Mp', 'Rp', 'Mstar', 'semimajor'
                var_class='physics': var = 'HX', 'species_list', 'atomic_masses', 'Z', 'ne', 'molec_adjust'
            var_end (float): Target value for the variable.
            var_class (str, optional): Class of the variable ('system' or 'physics'). Defaults to None.
            delta (float, optional): Default stepsize. Defaults to 0.02.
            delta_additive (bool, optional): If True, uses additive steps. Defaults to False.
            converge_bcs (bool, optional): If True, converges boundary conditions during ramping. Defaults to False.
            make_plot (bool, optional): If True, plots ramping progress. Defaults to True.
            expedite (bool, optional): If True, expedites ramping. Defaults to False.
            integrate_out (bool, optional): If True, integrates out after ramping. Defaults to True.
            static_bcs (bool, optional): If True, uses static boundary conditions. Defaults to False.

        Returns:
            int: Status code. 0 for success, other values for failure modes.
        """
        if var_class is None:
            var_class = self.ramp_class
        if var_class == "system":
            var_val = self.system.value(var)
        elif var_class == "physics":
            var_val = self.physics.value(var)
        else:
            self._normal_print("Unrecognized variable class for ramping: %s"
                  .format(var_class))
            return -100
        
        if (var == "Mp") or (var=='Mstar') or (var == "atomic_mass"):
            var_unit = 'g'
        elif (var == "Rp") or (var == "semimajor"):
            var_unit = 'cm'
        elif var == "Ftot":
            var_unit = 'ergs/s/cm2'
        elif var == "Lstar":
            var_unit = 'ergs/s'
        elif var == "molec_adjust":
            var_unit = '[m_H]'
        else:
            var_unit = ''
        print(f'\rRamping {str(var)} from {float(var_val):.3e} to {float(var_end):.3e} {str(var_unit)}.',
              end='                                                     \n')
        self.last_print_rastered = True
        if var_val == var_end or abs(var_val-var_end)/var_end < 1e-10:
            print(f'  {var:s} already done.',
                  end='                                                 \n')
            return 0
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        # Make sure inputfile matches the planet's parameters and hasn't changed
        self.inputs.write_planet_params(*self.windsoln.planet_tuple)
        self.inputs.write_physics_params(*self.windsoln.physics_tuple)
        if expedite==True:
            # Turn off integrating inward and outward
            expedite_flag_tuple = self.windsoln.flags_tuple
            original_flag_tuple = [x for x in expedite_flag_tuple]
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool
            expedite_flag_tuple[3] = 0
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out=False)
        # Check if ramping up or down
        flip = 1
        if var_val > var_end:
            flip *= -1
        delta *= flip
        if make_plot:
            # Setup plots for real time analysis of ramping
            fig, ax = plt.subplots(2, 2)
            # Plot current state before ramping
            max_sonic = self.windsoln.R_sp
            four_panel_plot(self.windsoln, ax,
                            label=(f'{var_val:.2e}'),
                            sub_sonic=True, past_rmin=True)
            fig.tight_layout()
            fig.canvas.draw()
            # Use percentage change to show change
            # Every n steps to prevent clutter
            plot_cntr = 0
            plot_every_pc = 0.10  # percent change between plots
            plot_every_n = 5  # successful steps between plots
        # Use deepcopy of current system to ramp to final system
        temp = None
        if var_class == "system":
            temp = self.system.deepcopy()
        elif var_class == "physics":
            temp = self.physics.deepcopy()
        failed = 0
        prcnt_chng = 1.0
        conv_cntr = 0
        conv_every_pc = 0.25
        conv_every_n = 10
        while (var_val != var_end):
            
            # Ramp temp system and write to planet parameters
            val_temp = var_val
            if delta_additive:
                val_temp += delta
            else:
                val_temp *= 1+delta
            if (flip*val_temp > flip*var_end):
                val_temp = var_end
            temp.assign(var, val_temp)
            if var_class == "system":
                self.inputs.write_planet_params(*temp.system_tuple())
            elif var_class == "physics":
                self.inputs.write_physics_params(*temp.physics_tuple())
            print("\r  Trying: {:s}:{:.6e}, delta:{:.4g}                   "
                  .format(var, val_temp, flip*(val_temp-var_val)/var_val),
                   end="                                                   ")
            self.last_print_rastered = True
            # See if can relax to partially ramped system
#             if converge_bcs:
#                 self.converge_mol_atomic_transition() #TEST - can afford to run each step
            result = self.run_wind(expedite=True,calc_postfacto=False)
            if result != 0:
                failed += 1
                delta /= 2.
                # Try converging. Can be expensive, so try only when smaller steps isn't working.
#                 if failed == 4:
#                     self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
#                                  10, #an arbitrary value that works well for ramping
#                                  self.windsoln.bcs_tuple[-1])
                if failed == 5:
#                     self.isotherm_start(run_wind=False) #updates where bolometric heat/cool dominate
                    if static_bcs:
                        # Converge boundary conditions on partial solution
                        print("\n   ...Intermediate ramping BCs activated. Static_bcs=True, so will not ramp BCs (may affect ability to converge).")
                    else:
                        # self._raster_print("   ...Intermediate ramping BCs activated.\n")
                        print("\n   ...Intermediate ramping BCs activated.")
                    self.ramp_base_bcs(static_bcs=static_bcs,tolerance=0.05)
                    self.converge_mol_atomic_transition()
                    self.converge_Ncol_sp(expedite=True,quiet=True)
            elif result == 2:
                self._normal_print("\nFailed to copy. Not good.")
                return 2
            else:
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge_bcs and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    
                    # converge boundary conditions on partial solution
                    self.ramp_base_bcs(intermediate=True,tolerance=0.1,static_bcs=static_bcs) #only converge if >10% diff in goal BC               
                    self.converge_mol_atomic_transition()
#                     self.converge_Ncol_sp() #MAYBE
                    conv_cntr = 0
                print("\r  Success %s:%.6e, delta:%.4g" %
                      (var, val_temp, flip*(val_temp-var_val)/var_val),
                       end="                                               ")
                self.last_print_rastered = True
                # update our system to partially ramped system
                var_val = val_temp
                if var_class == "system":
                    self.system.assign(var, var_val)
                elif var_class == "physics":
                    self.physics.assign(var, var_val)
                if make_plot:
                    # check if we should plot partial progress
                    if (prcnt_chng >= plot_every_pc and
                        plot_cntr >= plot_every_n):
                        max_sonic = max(max_sonic, self.windsoln.R_sp)
                        four_panel_plot(self.windsoln, ax,
                                        label=(f'{var_val:.2e}'),
                                       sub_sonic=True, past_rmin=True)
                        fig.tight_layout()
                        fig.canvas.draw()
                        plot_cntr = 0
                        prcnt_chng = 1.0
                    else:
                        prcnt_chng *= 1.+delta
                # Double delta if not ramping down faster than or equal to 50%
                if not (flip == -1 and delta >= 0.5):
                    delta *= 2
            if failed > 10:
                # Reset planet_params.inp to last sucessful state
                if var_class == "system":
                    self.inputs.write_planet_params(*self.system.system_tuple())
                elif var_class == "physics":
                    self.inputs.write_physics_params(*self.physics.physics_tuple())
                if expedite:
                    # return original flags
                    self.inputs.write_flags(*original_flag_tuple,integrate_out)
                print(f"\nERROR: Failing to converge {var}.")
                return 101
        print("\r  Final: ",end="\n")
        self.last_print_rastered = False
        # If did not plot last state of system do now
        if make_plot:
            if prcnt_chng != 1.0:
                four_panel_plot(self.windsoln, ax,
                                label=(f'{var_val:.2e}'),
                                sub_sonic=True, past_rmin=True)
            fig.tight_layout(pad=0.3)
            fig.suptitle(var)
            fig.subplots_adjust(bottom=0.3, top=0.9)
            fig.legend(bbox_to_anchor=(0.5, 0.0), loc='lower center', ncol=2)
            fig.canvas.draw()
            fig.canvas.flush_events()
        # Integrate inwards and outwards
        if expedite: #b/c expediting turned off outward integration before
            expedite_flag_tuple = self.windsoln.flags_tuple
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool 
            expedite_flag_tuple[3] = 1
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out) #user defined integrate out will override
        if integrate_out == True:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out) 
        result = self.run_wind(verbose=True,calc_postfacto=False)

        if result == 0:
            if converge_bcs:
                polish = self.polish_bcs(integrate_out,static_bcs)
                if polish == 0:
                    return 0
                if polish == 5:
                    self._normal_print("\nRamping successful, but polishing BCs to self consistency failed.")
                    return 5
            else:
                return 0
        elif result == 4:
            self._normal_print("\nRamping successful, but unable to integrate past sonic point. Option: Try polish_bcs().")
            return 4
        else:
            self._normal_print(f'\nERROR: Failure at the end game... result={result}')
            return 1
    

    def ramp_grav(self, system, delta=0.02, converge_bcs=False, make_plot=True,
                  expedite=False, integrate_out=True, static_bcs=False):
        """
        Ramps planet mass and radius along lines of constant surface gravity with adaptive stepsizes.

        Args:
            system (system): System object with Mp, Rp, Mstar, semimajor, Ftot, Lstar in cgs units.
            delta (float, optional): Default stepsize. Defaults to 0.02.
            converge_bcs (bool, optional): If True, runs polish_bcs(), which can be costly but may improve likelihood of convergence. Defaults to False.
            make_plot (bool, optional): If True, plots temperature, velocity, density, and ionization fraction profiles at steps throughout the ramp. Defaults to True.
            expedite (bool, optional): If True, skips integrating out to improve speed. Defaults to False.
            integrate_out (bool, optional): If True, requires integrate out to R_cori at the end of the ramp. Defaults to True.
            static_bcs (bool, optional): If True, skips running ramp_base_bcs() when ramper is having difficulty converging. Useful when not wanting to use default BCs. Defaults to False.

        Returns:
            int: Status code.
                0 - Successfully ramped (and integrated out and/or converged BCs)
                101 - Failed to ramp to desired Mp, Rp
                1 - (Unusual) Successfully ramped to Mp, Rp then encountered trouble solving. Last working solution should be sufficient.
                4 - Successfully ramped to Mp, Rp, but unable to integrate out. Try running converge_Rmax() or integrate_out().
                5 - Successfully ramped to Mp, Rp, but unable to polish all BCs. See printout for more details.
        """
        var_Mp = self.system.value('Mp')
        var_Rp = self.system.value('Rp')
        var_sg = var_Mp/var_Rp**2
        srt_Mp = var_Mp
        srt_Rp = var_Rp
        srt_sg = var_sg
        end_Mp = system.value('Mp')
        end_Rp = system.value('Rp')
        end_sg = end_Mp/end_Rp**2
        prct_Mp = abs((end_Mp - srt_Mp)/end_Mp)
        prct_Rp = abs((end_Rp - srt_Rp)/end_Rp)
        print("\rRamping {:s} from {:.3e} to {:.3e} g AND "
              "{:s} from {:.3e} to {:.3e} cm."
              .format('Mp', srt_Mp, end_Mp, 'Rp', srt_Rp, end_Rp), end='\n')
        self.last_print_rastered = True

        if srt_Mp != end_Mp:
            slope = (end_sg/srt_sg-1.)/(end_Mp/srt_Mp-1.)
        else:
            slope = 0
        if (var_Mp == end_Mp or abs(var_Mp-end_Mp)/end_Mp < 1e-10):
            print("  {:s} already done.".format('Mp'), end='\n')
        if (var_Rp == end_Rp or abs(var_Rp-end_Rp)/end_Rp < 1e-10):
            print("  {:s} already done.".format('Rp'), end='\n')
        if ((var_Mp == end_Mp or abs(var_Mp-end_Mp)/end_Mp < 1e-10) and
            (var_Rp == end_Rp or abs(var_Rp-end_Rp)/end_Rp < 1e-10)):
            return 0
        # Make sure inputfile matches the planet's parameters and hasn't changed
        self.inputs.write_planet_params(*self.windsoln.planet_tuple)
        self.inputs.write_physics_params(*self.windsoln.physics_tuple)
        if expedite==True:
            # Turn off integrating inward and outward
            expedite_flag_tuple = self.windsoln.flags_tuple
            original_flag_tuple = [x for x in expedite_flag_tuple]
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool
            expedite_flag_tuple[3] = 0
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out=False)
        # Check if ramping up or down
        M_flip = 1
        if srt_Mp > end_Mp:
            M_flip *= -1
        R_flip = 1
        if srt_Rp > end_Rp:
            R_flip *= -1
        M_delta = delta*M_flip
        R_delta = delta*R_flip
        if make_plot:
            # Setup plots for real time analysis of ramping
            fig, ax = plt.subplots(2, 2)
            # Plot current state before ramping
            max_sonic = self.windsoln.R_sp
            four_panel_plot(self.windsoln, ax,
                            label=(f'({var_Mp:.2e}, {var_Rp:.2e})'),
                            sub_sonic=True, past_rmin=True)
            fig.tight_layout()
            fig.canvas.draw()
            # Use percentage change to show change
            # Every n steps to prevent clutter
            plot_cntr = 0
            plot_every_pc = 0.10  # percent change between plots
            plot_every_n = 5  # successful steps between plots
        # Use deepcopy of current system to ramp to final system
        temp = self.system.deepcopy()
        failed = 0
        prcnt_chng = 1.0
        conv_cntr = 0
        conv_every_pc = 0.25
        conv_every_n = 10
        while (var_Mp != end_Mp) or (var_Rp != end_Rp):
            self.ramp_base_bcs(intermediate=True,tolerance=0.1,static_bcs=static_bcs) #only converge if >10% diff in goal BC  
            if (var_Mp != end_Mp) and (var_Rp != end_Rp):
                # Ramp Mp and Rp along linear rate of surface gravity change
                temp_Mp = var_Mp*(1.+M_delta)
                if (M_flip*temp_Mp >= M_flip*end_Mp):
                    temp_Mp = end_Mp
                    temp_Rp = end_Rp
                temp_sg = slope*(temp_Mp/srt_Mp-1.)+1.
                temp_sg *= srt_sg
                temp_Rp = (temp_Mp/temp_sg)**(0.5)
                temp.assign('Mp', temp_Mp)
                temp.assign('Rp', temp_Rp)
                print("\r  Trying: {:s}:{:.6e} & {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mp', temp_Mp, 'Rp', temp_Rp,
                              M_flip*(temp_Mp-var_Mp)/var_Mp),
                       end="                                               ")
                self.last_print_rastered = True
            elif var_Mp != end_Mp:
                # If we never needed to ramp Rp, ramp Mp linearly
                temp_Mp = var_Mp*(1.+M_delta)
                if (M_flip*temp_Mp >= M_flip*end_Mp):
                    temp_Mp = end_Mp
                temp.assign('Mp', temp_Mp)
                print("\r  Trying: {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mp', temp_Mp, M_flip*(temp_Mp-var_Mp)/var_Mp),
                       end="                                               ")
                self.last_print_rastered = True
            elif var_Rp != end_Rp:
                # If we never needed to ramp Mp, ramp Rp linearly
                temp_Rp = var_Rp*(1.+R_delta)
                if (R_flip*temp_Rp >= R_flip*end_Rp):
                    temp_Rp = end_Rp
                temp.assign('Rp', temp_Rp)
                print("\r  Trying: {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Rp', temp_Rp, R_flip*(temp_Rp-var_Rp)/var_Rp),
                       end="                                               ")
                self.last_print_rastered = True
            else:
                print("\nERROR: should be impossible surf_grav ramp condition")
                return
            # Write ramped variables to input file and try updating relaxation
            self.inputs.write_planet_params(*temp.system_tuple())
            result = self.run_wind(expedite=True,calc_postfacto=False)
            if result != 0:
                failed += 1
                if var_Mp != end_Mp:
                    M_delta /= 2.
                else:
                    R_delta /= 2.
                # Try converging. Can be expensive, so only try every once in a while
                if failed == 4:
                    if static_bcs:
                        print("\n  ...Intermediate ramping BCs activated. Static_bcs=True, so will not ramp BCs (may affect ability to converge).")
                    else:
                        # self._raster_print(" ...Intermediate ramping BCs activated.\n")
                        print("\n  ...Intermediate ramping BCs activated.")
                    self.converge_mol_atomic_transition()
                    self.ramp_base_bcs(static_bcs=static_bcs,tolerance=0.05) 
                    self.converge_mol_atomic_transition()
#                     self.converge_Ncol_sp()
            elif result == 2:
                print("\nERROR: Failed to copy windsoln to guess, not great.")
                return 2
            else: # Successfully relaxed with updated parameters
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge_bcs and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    # converge boundary conditions on partial solution
                    self.ramp_base_bcs(intermediate=True,tolerance=0.1,static_bcs=static_bcs) #only converge if >10% diff in goal BC
                    self.converge_mol_atomic_transition()
                    conv_cntr = 0
                # update our system to partially ramped system
                if (var_Mp != end_Mp) and (var_Rp != end_Rp):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mp', temp_Mp, 'Rp', temp_Rp,
                            M_flip*(temp_Mp-var_Mp)/var_Mp),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Mp = temp_Mp
                    var_Rp = temp_Rp
                    self.system.assign('Mp', var_Mp)
                    self.system.assign('Rp', var_Rp)
                elif var_Mp != end_Mp:
                    print("\r  Success %s:%.6e, M_delta:%.4g"
                          %('Mp', temp_Mp, M_flip*(temp_Mp-var_Mp)/var_Mp),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Mp = temp_Mp
                    self.system.assign('Mp', var_Mp)
                elif var_Rp != end_Rp:
                    print("\r  Success %s:%.6e, R_delta:%.4g"
                          %('Rp', temp_Rp, R_flip*(temp_Rp-var_Rp)/var_Rp),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Rp = temp_Rp
                    self.system.assign('Rp', var_Rp)
                if make_plot:
                    # check if we should plot partial progress
                    if (prcnt_chng >= plot_every_pc and
                        plot_cntr >= plot_every_n):
                        max_sonic = max(max_sonic, self.windsoln.R_sp)
                        four_panel_plot(self.windsoln, ax,
                                        label=(f'({var_Mp:.2e}, {var_Rp:.2e})'),
                                        sub_sonic=True, past_rmin=True)
                        fig.tight_layout()
                        fig.canvas.draw()
                        plot_cntr = 0
                        prcnt_chng = 1.0
                    else:
                        if srt_Mp != end_Mp:
                            prcnt_chng *= 1.+M_delta
                        else:
                            prcnt_chng *= 1.+R_delta
                # Double delta if not ramping down faster than or equal to 50%
                if not (M_flip == -1 and M_delta >= 0.5):
                    M_delta *= 2
                if not (R_flip == -1 and R_delta >= 0.5):
                    R_delta *= 2
            if failed > 10:
                # Reset planet_params.inp to last successful state
                self.inputs.write_planet_params(*self.system.system_tuple())
                if expedite:
                    # return original flags
                    self.inputs.write_flags(*original_flag_tuple,integrate_out)
                self._normal_print("\nERROR: Failing to converge gravity ramp.")
                return 101
        print("\r  Final: ",end="\n")
        self.last_print_rastered = False
        # If did not plot last state of system do now
        if make_plot:
            if prcnt_chng != 1.0:
                four_panel_plot(self.windsoln, ax,
                                label=(f'({var_Mp:.2e}, {var_Rp:.2e})'),
                                sub_sonic=True, past_rmin=True)
            fig.tight_layout(pad=0.3)
            fig.suptitle('(Mp, Rp)')
            fig.subplots_adjust(bottom=0.3, top=0.9)
            fig.legend(bbox_to_anchor=(0.5, 0.0), loc='lower center', ncol=2)
            fig.canvas.draw()
            fig.canvas.flush_events()
        # Integrate inwards and outwards
        if expedite: #b/c expediting turned off outward integration before
            expedite_flag_tuple = self.windsoln.flags_tuple
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool 
            expedite_flag_tuple[3] = 1
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out) #user defined integrate out will override
        result = self.run_wind(verbose=True,calc_postfacto=False)

        if result == 0:
            if converge_bcs:
                polish = self.polish_bcs(integrate_out,static_bcs)
                if polish == 0:
                    return 0
                if polish == 5:
                    self._normal_print("\nRamping successful, but polishing BCs to self consistency failed.")
                    return 5
            else:
                return 0
        elif result == 4:
            self._normal_print("\nRamping successful, but unable to integrate past sonic point. Option: Try polish_bcs().")
            return 4
        else:
            self._normal_print(f'\nERROR: Failure at the end game... result={result}')
            return 1
        

    def ramp_star(self, system, delta=0.02, converge_bcs=True, make_plot=True,
                  expedite=False, integrate_out=True, static_bcs=False):
        """
        Ramps stellar mass (Mstar), semimajor axis (semimajor), and stellar bolometric luminosity (Lstar).
        If Mstar and semimajor axis change, ramps linearly along Hill radius rate of change.
        If semimajor axis and Lstar change, ramps linearly along optical flux (essentially along skin temperature since T_skin propto Fopt) rate of change.

        Args:
            system (system): System object that takes Mp, Rp, Mstar, semimajor, Ftot, Lstar in cgs units.
            delta (float, optional): Default stepsize. Defaults to 0.02. (Unlikely to need to change)
            converge_bcs (bool, optional): If True, runs polish_bcs(), which can be costly but may improve likelihood of convergence. Defaults to True.
            make_plot (bool, optional): If True, plots temperature, velocity, density, and ionization fraction profiles at steps throughout the ramp. Defaults to True.
            expedite (bool, optional): If True, skips integrating out to improve speed. Defaults to False.
            integrate_out (bool, optional): If True, requires integrate out to R_cori at the end of the ramp. Defaults to True.
                Note: R_cori will still need to be computed by running converge_Rmax().
            static_bcs (bool, optional): If True, skips running ramp_base_bcs() when ramper is having difficulty converging. Useful when not wanting to use default BCs. Defaults to False.

        Returns:
            int: Status code.
                0 - Successfully ramped (and integrated out and/or converged BCs)
                101 - Failed to ramp to desired Mstar, Lstar, semimajor
                1 - Successfully ramped to Mstar, Lstar, semimajor then encountered trouble solving. Last working solution should be sufficient (self.path+'saves/windsoln.csv').
                4 - Successfully ramped to Mstar, Lstar, semimajor but unable to integrate out. Try running converge_Rmax() or integrate_out().
                5 - Successfully ramped to Mstar, Lstar, semimajor but unable to polish all BCs. See printout for more details.
        """
        start = time.time()
        count = 0
        var_Mstar = self.system.value('Mstar')
        var_adist = self.system.value('semimajor')
        var_Lstar = self.system.value('Lstar')
        var_R_hill = var_adist/var_Mstar**(1./3.)
        var_Fopt  = var_Lstar/var_adist**2 #want to also keep this const
        srt_Mstar = var_Mstar
        srt_adist = var_adist
        srt_Lstar = var_Lstar
        srt_Fopt  = var_Fopt
        srt_R_hill = var_R_hill
        end_Mstar = system.value('Mstar')
        end_adist = system.value('semimajor')
        end_Lstar = system.value('Lstar')
        end_Fopt  = end_Lstar/end_adist**2
        end_R_hill = end_adist/end_Mstar**(1./3.)
        print("\rRamping {:s} from {:.3e} to {:.3e} AND "
              "{:s} from {:.3e} to {:.3e} AND {:s} from {:.3e} to {:.3e}."
              .format('Mstar', srt_Mstar, end_Mstar,
                      'semimajor', srt_adist, end_adist,
                     'Lstar', srt_Lstar, end_Lstar), end='\n')
        self.last_print_rastered = True
        if srt_Mstar != end_Mstar:
            slope = (end_R_hill/srt_R_hill-1.)/(end_Mstar/srt_Mstar-1.)
        else:
            slope = 0
        if (srt_Lstar != end_Lstar) and (srt_adist != end_adist):
            slope_F = (end_Fopt/srt_Fopt-1.)/(end_adist/srt_adist-1.)
        else:
            slope_F = 0
        #checking if already done
        if (var_Mstar == end_Mstar or
            abs(var_Mstar-end_Mstar)/end_Mstar < 1e-10):
            print("  {:s} already done.".format('Mstar'), end='\n')
        if (var_adist == end_adist or
            abs(var_adist-end_adist)/end_adist < 1e-10):
            print("  {:s} already done.".format('semimajor'), end='\n')
        if (var_Lstar == end_Lstar or
            abs(var_Lstar-end_Lstar)/end_Lstar < 1e-10):
            print("  {:s} already done.".format('Lstar'), end='\n')
        if ((var_Mstar == end_Mstar or
             abs(var_Mstar-end_Mstar)/end_Mstar < 1e-10) and
            (var_adist == end_adist or
             abs(var_adist-end_adist)/end_adist < 1e-10) and
            (var_Lstar == end_Lstar or
             abs(var_Lstar-end_Lstar)/end_Lstar < 1e-10)):
            return 0
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        # Make sure inputfile matches the planet's parameters and hasn't changed
        self.inputs.write_planet_params(*self.windsoln.planet_tuple)
        self.inputs.write_physics_params(*self.windsoln.physics_tuple)
        if expedite == True:
            # Turn off integrating inward and outward
            expedite_flag_tuple = self.windsoln.flags_tuple
            original_flag_tuple = [x for x in expedite_flag_tuple]
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool
            expedite_flag_tuple[3] = 0
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out=False)
        # Check if ramping up or down
        M_flip = 1
        if srt_Mstar > end_Mstar:
            M_flip *= -1
        R_flip = 1
        if srt_adist > end_adist:
            R_flip *= -1
        F_flip = 1
        if srt_Fopt > end_Fopt:
            F_flip *= -1
        L_flip = 1
        if srt_Lstar > end_Lstar:
            L_flip *= -1   
        M_delta = delta*M_flip
        R_delta = delta*R_flip
        F_delta = delta*F_flip
        L_delta = delta*L_flip
        if make_plot:
            # Setup plots for real time analysis of ramping
            fig, ax = plt.subplots(2, 2)
            # Plot current state before ramping
            max_sonic = self.windsoln.R_sp
            four_panel_plot(self.windsoln, ax,
                            label=(f'({var_Mstar:.2e}, {var_adist:.2e})'),
                            sub_sonic=True, past_rmin=True)
            fig.tight_layout()
            fig.canvas.draw()
            # Use percentage change to show change
            # Every n steps to prevent clutter
            plot_cntr = 0
            plot_every_pc = 0.10  # percent change between plots
            plot_every_n = 5  # successful steps between plots
        # Use deepcopy of current system to ramp to final system
        temp = self.system.deepcopy()
        failed = 0
        prcnt_chng = 1.0
        conv_cntr = 0
        conv_every_pc = 0.25
        conv_every_n = 10
        while (var_Mstar != end_Mstar) or (var_adist != end_adist) or (var_Lstar != end_Lstar):
            #If all variables need to be ramped
            if (var_Mstar != end_Mstar) and (var_adist != end_adist) and (var_Lstar != end_Lstar):
                # Ramp Mstar and adist along linear rate of Hill radius change
                temp_Mstar = var_Mstar*(1.+M_delta)
                if (M_flip*temp_Mstar >= M_flip*end_Mstar):
                    temp_Mstar = end_Mstar
                    temp_adist = end_adist
                    temp_Lstar = end_Lstar
                temp_R_hill = slope*(temp_Mstar/srt_Mstar-1.)+1.
                temp_R_hill *= srt_R_hill
                temp_adist = temp_R_hill*(temp_Mstar)**(1./3.)
                temp.assign('Mstar', temp_Mstar)
                temp.assign('semimajor', temp_adist)
                temp_Fopt = slope_F*(temp_adist/srt_adist-1.)+1.
                temp_Fopt *= srt_Fopt
                temp_Lstar = temp_Fopt*(temp_adist**2)
                temp.assign('Lstar', temp_Lstar)
                print("\r  Trying: {:s}:{:.6e} & {:s}:{:.6e} &  {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mstar', temp_Mstar, 'semimajor', temp_adist,'Lstar', temp_Lstar,
                              M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                       end="                                                                          ")
                self.last_print_rastered = True
            #If bolometric luminosity is constant, but M* and a are not
            elif (var_Mstar != end_Mstar) and (var_adist != end_adist) and (abs(var_Lstar-end_Lstar)/end_Lstar < 1e-10):
                # Ramp Mstar and adist along linear rate of Hill radius change
                temp_Mstar = var_Mstar*(1.+M_delta)
                if (M_flip*temp_Mstar >= M_flip*end_Mstar):
                    temp_Mstar = end_Mstar
                    temp_adist = end_adist
                temp_R_hill = slope*(temp_Mstar/srt_Mstar-1.)+1.
                temp_R_hill *= srt_R_hill
                temp_adist = temp_R_hill*(temp_Mstar)**(1./3.)
                temp.assign('Mstar', temp_Mstar)
                temp.assign('semimajor', temp_adist)
                print("\r  Trying: {:s}:{:.6e} & {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mstar', temp_Mstar, 'semimajor', temp_adist,
                              M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                       end="                                                    ")
                self.last_print_rastered = True
            #If M* is const, but L* and a are not
            elif (abs(var_Mstar-end_Mstar)/end_Mstar < 1e-10) and (var_adist != end_adist) and (var_Lstar != end_Lstar):
                # Ramp Mstar and adist along linear rate of Hill radius change
                temp_adist = var_adist*(1.+R_delta)
                if (R_flip*temp_adist >= R_flip*end_adist):
                    temp_adist = end_adist
                    temp_Lstar = end_Lstar
                temp_Fopt = slope_F*(temp_adist/srt_adist-1.)+1.
                temp_Fopt *= srt_Fopt
                temp_Lstar = temp_Fopt*(temp_adist**2) #will this be going in the right direction? Yes, should be
                temp.assign('semimajor', temp_adist)
                temp.assign('Lstar', temp_Lstar)
                print("\r  Trying: {:s}:{:.6e} & {:s}:{:.6e}, a_delta:{:.4g}"
                      .format('Lstar', temp_Lstar, 'semimajor', temp_adist,
                              R_flip*(temp_adist-var_adist)/var_adist),
                       end="                                                                            ") 
                self.last_print_rastered = True
            #If semimajor is const, but M* and L* are not. These values are independent, but this just covers this case
            elif (var_Mstar!=end_Mstar) and (abs(var_adist-end_adist)/end_adist < 1e-10) and (var_Lstar != end_Lstar):
                # If we never needed to ramp adist, ramp Lstar linearly
                temp_Lstar = var_Lstar*(1.+L_delta)
                if (L_flip*temp_Lstar >= L_flip*end_Lstar):
                    temp_Lstar = end_Lstar
                temp.assign('Lstar', temp_Lstar)
                # If we never needed to ramp adist, ramp Mstar linearly
                temp_Mstar = var_Mstar*(1.+M_delta)
                if (M_flip*temp_Mstar >= M_flip*end_Mstar):
                    temp_Mstar = end_Mstar
                temp.assign('Mstar', temp_Mstar)
                print("\r  Trying: {:s}:{:.6e}, L_delta:{:.4g}; {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Lstar', temp_Lstar,L_flip*(temp_Lstar-var_Lstar)/var_Lstar,
                              'Mstar', temp_Mstar,
                              M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                       end="                                                                ") 
                self.last_print_rastered = True
           #If only M* needs to be ramped
            elif (var_Mstar != end_Mstar) and (abs(var_adist-end_adist)/end_adist < 1e-10) and (abs(var_Lstar-end_Lstar)/end_Lstar < 1e-10):
                # If we never needed to ramp adist, ramp Mstar linearly
                temp_Mstar = var_Mstar*(1.+M_delta)
                if (M_flip*temp_Mstar >= M_flip*end_Mstar):
                    temp_Mstar = end_Mstar
                temp.assign('Mstar', temp_Mstar)
                print("\r  Trying: {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mstar', temp_Mstar,
                              M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                       end="                                               ")
                self.last_print_rastered = True
            #if only semimajor needs to be ramped
            elif (abs(var_Mstar-end_Mstar)/end_Mstar < 1e-10) and (var_adist != end_adist) and (abs(var_Lstar-end_Lstar)/end_Lstar < 1e-10):
                # If we never needed to ramp Mstar, ramp adist linearly
                temp_adist = var_adist*(1.+R_delta)
                if (R_flip*temp_adist >= R_flip*end_adist):
                    temp_adist = end_adist
                temp.assign('semimajor', temp_adist)
                print("\r  Trying: {:s}:{:.6e}, a_delta:{:.4g}"
                      .format('semimajor', temp_adist,
                              R_flip*(temp_adist-var_adist)/var_adist),
                       end="                                               ")
                self.last_print_rastered = True
                self.ramp_base_bcs(intermediate=True,tolerance=0.01,static_bcs=static_bcs)
            #if only L* needs to be ramped
            elif (abs(var_Mstar-end_Mstar)/end_Mstar < 1e-10) and (abs(var_adist-end_adist)/end_adist < 1e-10) and (var_Lstar != end_Lstar):
                # If we never needed to ramp adist, ramp Lstar linearly
                temp_Lstar = var_Lstar*(1.+L_delta)
                if (L_flip*temp_Lstar >= L_flip*end_Lstar):
                    temp_Lstar = end_Lstar
                temp.assign('Lstar', temp_Lstar)
                print("\r  Trying: {:s}:{:.6e}, L_delta:{:.4g}"
                      .format('Lstar', temp_Lstar,
                              L_flip*(temp_Lstar-var_Lstar)/var_Lstar),
                       end="                                               ")
                self.last_print_rastered = True
            else:
                self._normal_print("\n ERROR: should be impossible surf_grav ramp condition")
                return 1
            # Write ramped variables to input file and try updating relaxation
            self.inputs.write_planet_params(*temp.system_tuple())
            result = self.run_wind(expedite=True,verbose=True,calc_postfacto=False)
            fail=0
            while result == 4:
                fail+=1
                print(f'   Temporarily increasing Ncol_sp for numerical integration reasons. Failure {fail:d}',
                                   end='                                                                          ')
                self.windsoln.bcs_tuple[5][:] = self.windsoln.Ncol_sp*1.2 
                self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                if fail > 10:
                    self._normal_print("\nERROR: Failed to integrate outwards too many times. ")
                    return 1
                result = self.run_wind(expedite=True, calc_postfacto=False)

            if result == 1:  #too many iterations during relaxation, try smaller stepsize
                failed += 1
                if var_Mstar != end_Mstar:
                    M_delta /= 2.
                if var_adist != end_adist:
                    R_delta /= 2.
                if var_Lstar != end_Lstar:
                    L_delta /= 2.
                # Try converging. Can be expensive, so only do if really stuck.
                if failed == 4:
                    count += 1
                    if not static_bcs:
                        # Converge boundary conditions on partial solution
                        print("\n  ...Intermediate ramping BCs activated. Static_bcs=True, so will not ramp BCs (may affect ability to converge).")
                    else:
                        # self._raster_print("  ...Intermediate ramping BCs activated (Instance {count:.0f}).\n")
                        print("\n  ...Intermediate ramping BCs activated.")
                    self.converge_Ncol_sp(expedite=True,quiet=True)
                    self.ramp_base_bcs(static_bcs=static_bcs,tolerance=0.05) 
                    self.converge_mol_atomic_transition()
            elif result == 2:
                self._normal_print("\nERROR: Failed to copy windsoln to guess. Not great.")
                return 2
            elif result == 0: # Successfully relaxed with updated parameters
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge_bcs and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    # converge boundary conditions on partial solution
                    # First update atmosphere to update Rmin location
                    self.converge_mol_atomic_transition()
                    self.ramp_base_bcs(intermediate=True,tolerance=0.05,static_bcs=static_bcs) #only converge if >5% diff in goal BC
                    conv_cntr = 0
                # update our system to partially ramped system
                if (var_Mstar != end_Mstar) and (var_adist != end_adist) and (var_Lstar != end_Lstar):
                    print("\r  Success %s:%.6e & %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar, 'semimajor', temp_adist, 'Lstar', temp_Lstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                           end="                                                       ")
                    self.last_print_rastered = True
                    var_Mstar = temp_Mstar
                    var_adist = temp_adist
                    var_Lstar = temp_Lstar
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('Lstar', var_Lstar)
                    self.system.assign('semimajor', var_adist)
                elif (var_Mstar != end_Mstar) and (var_adist != end_adist):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar, 'semimajor', temp_adist,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Mstar = temp_Mstar
                    var_adist = temp_adist
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('semimajor', var_adist)
                elif (var_Lstar != end_Lstar) and (var_adist != end_adist):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Lstar', temp_Lstar, 'semimajor', temp_adist,
                            L_flip*(temp_Lstar-var_Lstar)/var_Lstar),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Lstar = temp_Lstar
                    var_adist = temp_adist
                    self.system.assign('Lstar', var_Lstar)
                    self.system.assign('semimajor', var_adist)
                elif (var_Lstar != end_Lstar) and (var_Mstar != end_Mstar):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Lstar', temp_Lstar, 'Mstar', temp_Mstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                           end="                                                ")
                    self.last_print_rastered = True
                    var_Mstar = temp_Mstar
                    var_Lstar = temp_Lstar
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('Lstar', var_Lstar)
                elif var_Mstar != end_Mstar:
                    print("\r  Success %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Mstar = temp_Mstar
                    self.system.assign('Mstar', var_Mstar)
                elif var_adist != end_adist:
                    print("\r  Success %s:%.6e, R_delta:%.4g"
                          %('semimajor', temp_adist,
                            R_flip*(temp_adist-var_adist)/var_adist),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_adist = temp_adist
                    self.system.assign('semimajor', var_adist)
                elif var_Lstar != end_Lstar:
                    print("\r  Success %s:%.6e, L_delta:%.4g"
                          %('Lstar', temp_Lstar,
                            L_flip*(temp_Lstar-var_Lstar)/var_Lstar),
                           end="                                           ")
                    self.last_print_rastered = True
                    var_Lstar = temp_Lstar
                    self.system.assign('Lstar', var_Lstar)
                if make_plot:
                    # check if we should plot partial progress
                    if (prcnt_chng >= plot_every_pc and
                        plot_cntr >= plot_every_n):
                        max_sonic = max(max_sonic, self.windsoln.R_sp)
                        four_panel_plot(self.windsoln, ax,
                                        label=(f'({var_Mstar:.2e}, '
                                               f'{var_adist:.2e})'),
                                        sub_sonic=True, past_rmin=True)
                        fig.tight_layout()
                        fig.canvas.draw()
                        plot_cntr = 0
                        prcnt_chng = 1.0
                    else:
                        if srt_Mstar != end_Mstar:
                            prcnt_chng *= 1.+M_delta
                        else:
                            prcnt_chng *= 1.+R_delta
                # Double delta if not ramping down faster than or equal to 50%
                if not (M_flip == -1 and M_delta >= 0.5):
                    M_delta *= 2
                if not (R_flip == -1 and R_delta >= 0.5):
                    R_delta *= 2
                if not (L_flip == -1 and L_delta >= 0.5):
                    L_delta *= 2
            if failed > 10:
                # Reset planet_params.inp to last sucessful state
                self.inputs.write_planet_params(*self.system.system_tuple())
                if expedite:
                    # return original flags
                    self.inputs.write_flags(*original_flag_tuple,integrate_out)
                self._normal_print("\nERROR: Failing to converge on star ramp.")
                return 101
        print("\r  Final: ",end="\n") 
        self.last_print_rastered = False
        # If did not plot last state of system do now
        if make_plot:
            if prcnt_chng != 1.0:
                four_panel_plot(self.windsoln, ax,
                                label=(f'({var_Mstar:.2e}, {var_adist:.2e})'),
                                sub_sonic=True, past_rmin=True)
            fig.tight_layout(pad=0.3)
            fig.suptitle('(Mstar, semimajor)')
            fig.subplots_adjust(bottom=0.3, top=0.9)
            fig.legend(bbox_to_anchor=(0.5, 0.0), loc='lower center', ncol=2)
            fig.canvas.draw()
            fig.canvas.flush_events()
        # Integrate inwards and outwards
        if expedite: #b/c expediting turned off outward integration before
            expedite_flag_tuple = self.windsoln.flags_tuple
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool 
            expedite_flag_tuple[3] = 1
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out) #user defined integrate out will override
        result = self.run_wind(calc_postfacto=False)

        if result == 0:
            if converge_bcs:
                polish = self.polish_bcs(integrate_out,static_bcs)
                if polish == 0:
                    return 0
                if polish == 5:
                    self._normal_print("\nRamping successful, but polishing BCs to self consistency failed.")
                    return 5
            else:
                return 0
        elif result == 4:
            self._normal_print("\nRamping successful, but unable to integrate past sonic point. Option: Try polish_bcs().")
            return 4
        else:
            self._normal_print(f'\nERROR: Failure at the end game... result={result}')
            return 1
    
    
#Metals functions    
    def calc_metallicity(self, metals_list=None, Z=1):
        """
        Returns an array of the Lodders (2009) Z x solar metallicity MASS fractions for the metals in metals_list.
        Mass fractions must sum to 1, so the difference from 1 is added to the hydrogen mass fraction (always at index 0 in the array).

        Args:
            metals_list (list of str): List of metals, e.g., ['H I', 'He I', 'Mg II'].
            Z (float): Multiple of solar metallicity.

        Returns:
            np.ndarray: Array of mass fractions for each species in metals_list.
        """
        if metals_list is None:
            metals_list = self.windsoln.species_list
        self.metals = metal_class(self.windsoln)
        return self.metals.metallicity(metals_list,Z)
        
    
    def add_metals(self, desired_species_list, Z=1, custom_mfs=[], Ncol_sp_tot=0, 
                   integrate_out=True, converge_Ncol_sp=True):
        """
        Adds NEW metals to existing solutions. To ramp mass fractions or metallicity of species already in the solution, use ramp_metallicity(). Skips adding if species is already present. This function will converge to a self-consistent Ncol_sp for each species at the final step if converge_Ncol_sp=True.

        Args:
            desired_species_list (list of str): List of all species desired in the solution. Give element + ionization state (e.g., 'HI', 'h1', 'h I' for neutral hydrogen).
            Z (float, optional): Metallicity in units of solar metallicity. Defaults to 1.
            custom_mfs (list or array, optional): User-specified mass fractions. Defaults to [].
            Ncol_sp_tot (int, optional): Total sonic point column density across all species. Defaults to 0.
            integrate_out (bool, optional): If True, integrates past sonic point to Coriolis radius. Defaults to True.
            converge_Ncol_sp (bool, optional): If True, converges Ncol_sp for each species. Defaults to True.

        Returns:
            int: Status code. 0 for success, other values for failure modes.
        """
        self.metals = metal_class(self.windsoln)
        
        desired_species = desired_species_list
        McAtom.formatting_species_list(desired_species)
        #check whether these species are in the Verner list in McAstro
        skips = np.array([], dtype=int)
        for i in range(len(desired_species)):
            check_species = desired_species[i]
            try:
                McAtom.atomic_species(check_species).mass.iloc[0]
            except IndexError:
                self._normal_print("%s is not currently available in Wind-AE. Skipping." %check_species)
                skips = np.append(skips,int(i))
        desired_species = np.delete(desired_species,skips)
        unspaced_desired_list = [sp.replace(' ','') for sp in desired_species]

        #If user defined custom mass fractions, use those values
        unspaced_current_list = [sp.replace(' ','') for sp in self.windsoln.species_list]
        new_species = np.setdiff1d(unspaced_desired_list,
                                   unspaced_current_list,
                                  assume_unique=True)
        
        if len(custom_mfs) > 0: 
            if len(custom_mfs) != (len(unspaced_current_list)+len(new_species)):
                sys.exit("ERROR: Mass fraction and species list must be the same length.")
            if np.round(np.sum(custom_mfs),5) != 1:
                self._normal_print('WARNING: Total Mass Fraction must sum to 1. sum(ZX) = %.3f' %np.sum(custom_mfs))
            species_mf_dict = dict(zip(unspaced_desired_list,custom_mfs))


        
        self.metals = metal_class(self.windsoln)
        self.metals.add_species_to_guess(new_species)
        if self.raise_Ncol_sp(by_factor=20) != 0:
            self._normal_print("\nAttempted to raise Ncol_sp to increase speed of convergence, but failed.")
        if self.run_wind(expedite=True)==0:
            self._raster_print("Metals added to guess, now ramping mass fractions.")
            #setting mass fractions to multiple of solar Z or to custom values
            if len(custom_mfs)==0:
                goal_mfs = self.metals.metallicity(self.windsoln.species_list,1)
            else:
                goal_mfs=[]
                for sp in self.windsoln.species_list:
                    goal_mfs = np.append(goal_mfs,species_mf_dict[sp])
            
            #Writing goal mass fractions to input files
            start_mfs = np.copy(self.windsoln.HX)
            #Attempting to jump in 1 step to goal mass fractions (usually works)
            ratio = goal_mfs/self.windsoln.atomic_masses
            goal_nfs = ratio/sum(ratio)
            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                 0.5*goal_nfs, #an arbitrary value that works well for ramping
                                 self.windsoln.bcs_tuple[-1])
            self.inputs.write_physics_params(goal_mfs,
                                            *self.windsoln.physics_tuple[1:])
            result = self.run_wind(expedite=True,calc_postfacto=False)
            if result != 0:
                self._normal_print("Failed to ramp mass fracs directly to goal. Taking 10% steps. ")
                percent = 0.1
                delta_mfs = (goal_mfs - start_mfs)*0.1
                step_mfs = self.windsoln.HX + delta_mfs
                step_mfs[0] = 1 - sum(step_mfs[1:])
                self.inputs.write_physics_params(step_mfs,
                                                *self.windsoln.physics_tuple[1:])
                # self._normal_print('------%.5e-------'%abs(np.mean((goal_mfs - self.windsoln.HX)/goal_mfs)))
                while abs(np.mean((goal_mfs - self.windsoln.HX)/goal_mfs)) > 1e-4:
                    failed = 0
                    self._raster_print(f'Avg. fractional difference from goal: {abs(np.mean((goal_mfs - self.windsoln.HX)/goal_mfs)):.3f}')
                    while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                        failed += 1
                        delta_mfs /= 2
                        percent /= 2

                        self._raster_print(f'   Fail {failed}: Trying smaller stepsize (x{percent:.3f})')
                        if failed > 8:
                            self._normal_print(f"\nStruggling to substep towards goal mass fracs: \nCurrent: {self.windsoln.HX} \nGoal: {goal_mfs}")
                            self._normal_print("Hint: Try polish_bcs(converge_Rmax=False) then running sim.ramp_metallicity() again.")
                            return 1
                        
                    step_mfs = self.windsoln.HX + delta_mfs
                    step_mfs[0] = 1 - sum(step_mfs[1:])
                    self.inputs.write_physics_params(step_mfs,
                                                    *self.windsoln.physics_tuple[1:])
                    # self._normal_print('------%.5e-------'%abs(np.mean((goal_mfs - self.windsoln.HX)/goal_mfs)))
                    
            if Z>1:
                self._normal_print(f"Metals successfully added in at Z=1. Ramping to Z={Z}.")
                return self.ramp_metallicity(Z=Z,integrate_out=integrate_out,
                                            converge_Ncol_sp=converge_Ncol_sp)
            elif integrate_out==True:
                if converge_Ncol_sp==True:
                    self._normal_print("Successfully ramped mass fracs. Integrating out & converging Ncol_sp...")
                    if self.converge_Rmax(final_converge_Ncol=False) != 0:
                        return 4
                    else:
                        return 0
                if converge_Ncol_sp==False:
                    self.inputs.write_flags(*self.windsoln.flags_tuple,
                                            integrate_out=True)
                    self._normal_print("Successfully ramped mass fractions. Integrating out, but not converging Ncol_sp.")
                    if self.run_wind(calc_postfacto=False) != 0:
                        return 4
                    else:
                        return 0
            else:
                if converge_Ncol_sp==True:
                    self._normal_print("Successfully ramped mass fracs. Converging Ncol_sp, but not integrating out...")
                    return self.converge_Ncol_sp()
                else:
                    self._raster_print("Successfully ramped mass fractions.")
                    return 0
        else:
            self._normal_print("\nFailed to add new species. That shouldn't have happened. Try polishing solution before adding metals.")
            return 1

        
    
    def remove_metals(self, remove_species_list, run_wind=True):
        """
        Removes a metal from the list. Note: Helium will take a significant amount of time to remove, so it is advisable to ramp directly to your desired solution from an existing pure-H solution.

        Args:
            remove_species_list (list of str): List of species to remove.
            run_wind (bool, optional): If True, runs wind relaxation after removal. Defaults to True.

        Returns:
            None
        """
        self.metals = metal_class(self.windsoln)
        try:
            self.metals.remove_species_from_guess(remove_species_list)
#             self.windsoln.nspecies -= len(remove_species_list)
            self.load_planet(self.path+'inputs/guess.inp',calc_postfacto=False,
                             print_atmo=False,print_warnings=False)
        except ValueError:
            self._normal_print('\nOne or more of the species you are attempting to remove is not present in simulation.')
            return 1
        if run_wind == True:
            if self.run_wind(expedite=True,calc_postfacto=False) != 0:
                warning = '''WARNING: Running sim with reduced number of species was unsuccessful.
                Consider reloading original solution and manually ramping down to 0 the mass 
                fractions of the species you want to remove. Recall that total mass fraction 
                must sum to 1.
                '''
                self._normal_print(warning)
            else:
                self._normal_print("\n"+','.join('%s' %sp for sp in remove_species_list)+" removed and new windsoln generated.")         
        return

    def ramp_metallicity(self, Z=1, custom_mfs=[], static_bcs=False,
                        integrate_out=True, converge_Ncol_sp=True, polish_bcs=False):
        """
        Ramps up the metallicity of the species present in the simulation. Can ramp in multiples of solar Z or set custom mass fractions.

        Args:
            Z (float, optional): Multiples of Lodders (2008) solar metallicity. Defaults to 1.
            custom_mfs (list of float, optional): Custom mass fractions. If empty, defaults to Z. If not empty, uses custom_mfs. Defaults to [].
            static_bcs (bool, optional): If True, does not ramp lower boundary conditions (may affect convergence). Defaults to False.
            integrate_out (bool, optional): If True, integrates outwards past the sonic point to the Coriolis radius. Defaults to True.
            converge_Ncol_sp (bool, optional): If True, attempts to converge the Ncol_sp parameter to self-consistency after solving (if integrate_out=False). Defaults to True.
            polish_bcs (bool, optional): If True, polishes all boundary conditions after solving. Defaults to False.

        Returns:
            int: Status code. 0 for success, other values for failure modes.
        """
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        goal_Z = Z
        if goal_Z >= 100:
            self._normal_print("Note: For high metallicities, consider increasing mean molecular weight (molec_adjust).")
            for sp in self.windsoln.species_list:
                if (sp[:2] == 'Fe') or (sp[:2] == 'Ca'):
                    self._normal_print("WARNING: At high Z, Fe and Ca line cooling are significant, but have not yet been implemented.")
            self._normal_print("\n")
            
        if len(custom_mfs) != 0:
            if len(custom_mfs) != len(self.windsoln.species_list):
                sys.exit("ERROR: Mass fraction and species list must be the same length.")
            if np.round(np.sum(custom_mfs),5) != 1:
                self._normal_print('WARNING: Total Mass Fraction must sum to 1. sum(ZX) = %.3f \n' %np.sum(custom_mfs))
            self._normal_print("NOTE: Goal mass fractions will override any goal Z metallicity provided.")
            goal_mass_fracs = custom_mfs
            
            current_HX = np.copy(self.windsoln.HX)
            ratio = goal_mass_fracs/self.windsoln.atomic_masses
            goal_nfs = ratio/sum(ratio)
            if goal_nfs[0] < 0.5:
                self._normal_print("NOTE: Secondary ionization calcs assume H is the dominant species, therefore may not be accurate.")
            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                 0.5*goal_nfs, #an arbitrary value that works well for ramping
                                 self.windsoln.bcs_tuple[-1])
            self.inputs.write_physics_params(goal_mass_fracs,*self.windsoln.physics_tuple[1:])
            self._raster_print(f' Goal {goal_mass_fracs}: Attempting {goal_mass_fracs}')
            if self.run_wind(expedite=True,calc_postfacto=False) == 0:
                self._normal_print(f'Mass fractions successfully ramped to: {self.windsoln.HX}')
                if integrate_out == True:
                    print("  Integrating out and attempting to converge Ncol_sp...")
                    return self.converge_Rmax()
                elif (integrate_out == False) & (converge_Ncol_sp==True):
                    print("  Attempting to converge Ncol_sp...")
                    return self.converge_Ncol_sp(expedite=False,quiet=True)
                elif converge_Ncol_sp == False:
                    print("Note: Ncol_sp has not been converged to self-consistency.")
                    return 0
            else:
                self._raster_print(' Taking smaller steps...')
                percent = 0.2
                delta = np.copy(goal_mass_fracs - self.windsoln.HX )
                while abs(sum((goal_mass_fracs - self.windsoln.HX)/self.windsoln.HX)) > 1e-6:               
                    step_mfs = self.windsoln.HX+percent*delta
                    self._raster_print(f' Goal {goal_mass_fracs}: Attempting {step_mfs}')
                    ratio = (step_mfs/self.windsoln.atomic_masses)
                    step_nfs = ratio/sum(ratio)
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                         0.5*step_nfs, #an arbitrary value that works well for ramping
                                         self.windsoln.bcs_tuple[-1]) 
                    # self.run_wind(expedite=True,calc_postfacto=False)
                    if sum(abs(percent*delta)) > sum(abs(goal_mass_fracs - self.windsoln.HX)):
                        self.inputs.write_physics_params(self.windsoln.HX+(goal_mass_fracs - self.windsoln.HX),
                                                         *self.windsoln.physics_tuple[1:])
                    else:
                        self.inputs.write_physics_params(step_mfs,*self.windsoln.physics_tuple[1:])
                 
                    fail = 0
                    while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                        fail+=1
                        if fail==3:
                            print("  ...Intermediate BC ramping activated.",
                                  self.converge_mol_atomic_transition(),self.ramp_base_bcs(intermediate=True,static_bcs=static_bcs,tolerance=0.05))     
                            self.raise_Ncol_sp(by_factor=10)                    
                        step = self.windsoln.HX+percent*delta/(fail+1)
                        self._raster_print(f'  Fail {fail:d}: Attempting {step}')
                        self.inputs.write_physics_params(step,
                                                         *self.windsoln.physics_tuple[1:]) 
                        if fail>10:
                            self._normal_print(f"Failed to ramp mass fractions. Current: {self.windsoln.HX}."
                                  f" Goal: {custom_mfs}")
                            return 1
                self._normal_print(f'Mass fractions successfully ramped to: {self.windsoln.HX}')
                
        else: 
            #Just getting a grid of metallicity to print what the current Z is 
            grid = np.zeros(2000)
            for i in range(2000):
                grid[i] = abs(self.calc_metallicity(self.windsoln.species_list,Z=i+1)[0]-
                              self.windsoln.HX[0])
            start_Z = np.where(grid==min(grid))[0][0]+1
            print("Starting metallicity: %d xSolar"%start_Z)
            #Ramping in metallicity space 
            current_Z = start_Z
            step_Z = current_Z
            while (1-current_Z/goal_Z) > 1e-5:
                if abs(goal_Z - current_Z) > 10:
                    if goal_Z > current_Z:
                        step_Z += 10
                    if goal_Z < current_Z:
                        step_Z -= 10
                elif abs(goal_Z - current_Z) < 10:
                    step_Z = goal_Z
                self._raster_print(f'  Attempting to ramp Z from {current_Z:.1f} to {step_Z:.1f}')
                step_mfs = self.metals.metallicity(self.windsoln.species_list,Z=step_Z)
                ratio = step_mfs/self.windsoln.atomic_masses
                step_nfs = ratio/sum(ratio)
                self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                     0.8*step_nfs, #an arbitrary value that works well for ramping
                                     self.windsoln.bcs_tuple[-1])
                self.run_wind(expedite=True,calc_postfacto=False)
                
                self.inputs.write_physics_params(step_mfs,self.windsoln.species_list,
                                                 self.windsoln.molec_adjust)
                fail = 1
                delta = goal_Z - current_Z
                while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                    step_Z = current_Z + delta/(2**fail)
                    self._raster_print(f' Failed {fail}: Attempting to ramp Z from {current_Z:.1f} to {step_Z:.1f}')
                    if fail>5:
                        self._normal_print('Failed at Z = %.1f'%current_Z)
                        sys.exit(1)
                        
                    step_mfs = self.metals.metallicity(self.windsoln.species_list,Z=step_Z)
                    self.inputs.write_physics_params(step_mfs,self.windsoln.species_list,
                                                     self.windsoln.molec_adjust)                    
                    fail+=1
                current_Z = step_Z    
            print(f'Success! Ramped to goal metallicity, Z = {current_Z:.0f} x Solar')

            if integrate_out == True:
                print("  Integrating out and attempting to converge Ncol_sp...")
                return self.converge_Rmax()
            elif (integrate_out == False) & (converge_Ncol_sp==True):
                self._normal_print("  Attempting to converge Ncol_sp...")
                return self.converge_Ncol_sp(expedite=False,quiet=True)
            elif converge_Ncol_sp == False:
                self._normal_print("  Note: Ncol_sp has not been converged to self-consistency.")
                return 0



#Polishing Boundary Conditions functions (many of these these enforce self-consistency, and 
#are not neccessary for precision, but are for maximal accuracy (within the inherent uncertainty 
#in the model))
    def polish_bcs(self, converge_Rmax=True, static_bcs=None, user_override_press=False,base_press=1,bolo_on=True):
        """
        Polishes upper and lower boundary conditions to self-consistency.

        Args:
            converge_Rmax (bool, optional): If True, converges Rmax (more costly). Defaults to True.
            static_bcs (bool, optional): If True, skips running ramp_base_bcs() and maintains the current boundary conditions. Defaults to None (i.e., the value of self.static_bcs).
            user_override_press (bool, optional): If True, allows user to override pressure boundary conditions. Defaults to False.
            base_press (float, optional): Base pressure in microbars (barye) to use for boundary conditions. Defaults to 1.

        Returns:
            int: 0 if successfully polished all BCs, 5 if failed at polishing one or more BCs. See printout for more details. If "Molec. to atomic transition" has failed, plot energy_plot() to assess whether the current solution is satisfactory.
        """
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        #Checking that base bcs (Rmin, rho, T) have been converged
        #--------------
        print('Polishing up boundary conditions...')
        self.ramp_base_bcs(static_bcs=static_bcs,polish=True,base_press=base_press,user_override_press=user_override_press) 
        goal_bcs = self.base_bcs(user_override_press=user_override_press,base_press=base_press)
        curr_bcs = np.array(self.windsoln.bcs_tuple[:4])   
        avg_diff = sum(np.array((abs(curr_bcs-goal_bcs)/goal_bcs)[[0,2,3]]))

        #Checking that the molecular to atomic transition is occuring at the correct radius
        #--------------
        isotherm = self.converge_mol_atomic_transition(polish=True,bolo_on=bolo_on)
        width,idx = self.erf_velocity(return_idx=True,polish=True)[-2:]
        
        #Does bolometric heating/cooling impede too far into wind? 
        #If so, make it a sharper drop off
        self.windsoln.add_user_vars(expedite=True) #do postfacto calcs to get 'heat_ion' & 'boloheat'
        while len(np.where(self.windsoln.soln['heat_ion'][idx:]<self.windsoln.soln['boloheat'][idx:])[0]) > 20:
            width /= 2
            self._raster_print(f"  ...Shortening transition from molecular to atomic region. Erf width = {width:.2f} Hsc")
            if self.converge_mol_atomic_transition(polish=True,width_factor=width) != 0:
                self._normal_print(f"Attempting to shorten erfc transition region between molecular and atomic regions failed at {width:.3f}Hsc. Check energy_plot() to ensure bolometric heating/cooling doesn't impede on photoionization heated region.")
                isotherm=1
                break
            idx = self.erf_velocity(return_idx=True,polish=True)[-1]
            self.windsoln.add_user_vars(expedite=True) #annoying but have to repopulate for some reason

        #If bolo heat/cooling drops off too early, it can induce an unphysical and numerically
        #unstable spike in advective cooling. This extends the range of bolo to smooth spike.
        # while any(-self.windsoln.soln['heat_advect'][100:idx+10]>
        #           self.windsoln.soln['heat_ion'][100:idx+10]):
        #     width+=1
        #     if width > 10:
        #         self._normal_print("Warning: 20 scaleheights is an unlikely width for the error function",
        #               " transitioning between molecular and atomic regions.",
        #               "\nCheck energy_plot(). Stopping here.")
        #         isotherm=1
        #         break
        #     self._raster_print(f"  ...Smoothing transition from molecular to atomic region. Erf width = {width:.2f} Hsc")
        #     isotherm = self.converge_mol_atomic_transition(polish=True,width_factor=width/2)
        #     if isotherm != 0:
        #         self._normal_print("Failed to smooth transition from molecular to atomic region."+self.clear,
        #               "Unphysical kinks in wind profile may be present at base of wind. ",
        #               "\n       Mass loss rate relatively unaffected.")
        #         isotherm=1
        #         break
        #     self.windsoln.add_user_vars(expedite=True)

       #Converging Rmax (sets Rmax=r_cori) and/or Ncol_sp
        rcori_result = 'Success'
        if converge_Rmax==True:
            Rmax = self.converge_Rmax(final_converge_Ncol=True)
            if Rmax == 1: #occurs when relax fails
                Ncol = 0
            if Rmax==0:
                Ncol=0 
            rcori_str = f'\n   Setting Rmax to R_coriolis:         %s'
        else:
            Ncol = self.converge_Ncol_sp()
            rcori_str = ''
                    
        warn = ''
        bc_result = 'Success'
        iso_result = 'Success'
        ncol_result = 'Success'
        warn2=f'''-\nPolishing boundary conditions:\n   Base boundary conditions:           %s\n   Molec. to atomic transition:        %s\n   Self-consist. Ncol at sonic point:  %s'''+rcori_str
        fails = 0
        if avg_diff > 1e-2:
            fails += 1
            warn += ' base boundary conditions (ramp_base_bcs),'
            if static_bcs==False and self.windsoln.bolo_heat_cool != 0:
                bc_result = 'Failed'
            elif static_bcs==True:
                bc_result = 'Held const.'
        if isotherm == 1 or (self.windsoln.bolo_heat_cool > 0.0 and self.windsoln.bolo_heat_cool < 1.0):
            fails += 1
            warn += ' bolometric heating/cooling & molecular layer (converge_mol_atomic_transition),'
            iso_result = 'Failed'
        else:
            self.converge_mol_atomic_transition(polish=True) #run again to ensure polished
        if Ncol != 0:
            fails += 1
            warn += ' Ncol at sonic point (converge_Ncol_sp),'
            ncol_result = 'Failed'
        try:
            if Rmax != 0:
                fails += 1
                warn += ' setting Rmax to R_coriolis (converge_Rmax)'
                rcori_result = 'Failed'
        except NameError:
            warn += ''
       
        if (len(warn)>0) and (warn[-1] == ','):
            warn = warn[:-1]
        
        if fails == 0:
            #enforces outward integration
            self.integrate_out(quiet=True)
            self._normal_print("  All boundary conditions succesfully polished to self-consistency.")
            return 0
        else:
            if converge_Rmax == True:
                self._normal_print(warn2 %(bc_result,iso_result,ncol_result,rcori_result))
            else:
                self._normal_print(warn2 %(bc_result,iso_result,ncol_result))
            return 5
        

        
    def converge_mol_atomic_transition(self,width_factor=0.,_called_in_ramp_bcs=False,polish=False, bolo_on=True):
        """
        _Formerly run_isotherm()_
        
        Ramps the complementary error function that governs the transition from the molecular to atomic regions in the atmosphere. Both the mean molecular weight and bolometric heating/cooling are transitioned using the same error function. 

        If the bolometric heating/cooling flag is off (sim.windsoln.bolo_heat_cool=0.0) there will be no bolometric heating/cooling, but the mean molecular weight will still transition from molecular to atomic unless sim.windsoln.molec_adjust = 0.0.
        To turn off the bolometric heating/cooling run sim.turn_off_bolo().

        Criterion for transition: when photoionization heating begins to dominate over the PdV cooling and a wind will launch.
        Below the wind, for an average planet, molecules have not photodissociated so mu should be mean molecular weight instead of mean atomic weight and the erf enforces this. Molecular opacities mean that bolometric heating and cooling dominate the energy budget and create an isotherm in the molecular region below the wind. In the optically thin atomic wind, bolometric heating and cooling are negligible, so the erf also enforces the drop off of bolometric heating and cooling.

        If a user needs to call this function for some reason, width_factor is the only argument they should need to change. For numerical or physical reasons it may be necessary for the transition to occur over more than 1 scaleheight. In this case, increase width_factor.

        Args:
            width_factor (float, optional): Factor to adjust the width in scaleheights of the transition region. Defaults to 0.
            _called_in_ramp_bcs (bool, optional): Set True if called within ramp_base_bcs to avoid recursion issues. Defaults to False.
            polish (bool, optional): If True, computes full postfacto heat_ion & cool_PdV for more accurate transition point. Defaults to False.
            bolo_on (bool, optional): If polish=True, bolo_on=True turns on bolometric heating/cooling if appropriate. Defaults to True.

        Returns:
            int or tuple: Index of transition if return_idx is True, otherwise status code.
        """
        #if the molecular layer has been turned off and not polishing
        if (self.windsoln.bolo_heat_cool == 0.0) and (polish==False):
            return 0
        #if we have stipulated it should be skipped
        if (polish==False) and (self.skip==True):
            return 0
        
        #loading last working solution to do this
        self.load_planet(self.path+'saves/windsoln.csv',calc_postfacto=False,
                         print_atmo=False,print_warnings=False)
    
        
        #if there is an unphysical spike in heat_advect at base of wind, widen erfc
        #to-do: probably should just do this in polish, not every time since costly
        if polish:
            #generate the goal erfc parameters
            v_drop,rate,width,idx = self.erf_velocity(_called_in_ramp_bcs=_called_in_ramp_bcs,
                                                        polish=polish,
                                                        width_factor=width_factor,
                                                        return_idx=True) 
            self.windsoln.add_user_vars(expedite=True)
            while any(-self.windsoln.soln['heat_advect'][100:idx+10]>
                    self.windsoln.soln['heat_ion'][100:idx+10]):
                width+=1
                if width > 10:
                    self._normal_print("Warning: 20 scaleheights is an unlikely width for the error function transitioning between molecular and atomic regions. \nCheck energy_plot(). Stopping here.")
                    #reset to last working solution
                    self.load_planet(self.path+'saves/windsoln.csv',calc_postfacto=False,
                                    print_atmo=False,print_warnings=False)
                    break
                self._raster_print(f"  ...Smoothing transition from molecular to atomic region to prevent numerical instabilities. Trying erf width = {width:.2f} Hsc")
                v_drop,rate = self.erf_velocity(width_factor=width,_called_in_ramp_bcs=_called_in_ramp_bcs,polish=polish)
                #for now, keep current radial erfc drop off location and just change width
                isotherm = self.ramp_molecular_erfc(v_drop=self.windsoln.erf_drop[0], rate=rate, polish=polish)
                if isotherm != 0:
                    # if polish == True:
                    self._normal_print("Failed to smooth transition from molecular to atomic region. Unphysical kinks in wind profile may be present at base of wind. \n       Mass loss rate relatively unaffected.")
                    self.load_planet(self.path+'saves/windsoln.csv',calc_postfacto=False,
                                    print_atmo=False,print_warnings=False)
                    break
                self.windsoln.add_user_vars(expedite=True)    
        
        # Now ramp the erfc to the goal values as usual
        v_drop,rate,width,idx = self.erf_velocity(_called_in_ramp_bcs=_called_in_ramp_bcs,
                                polish=polish,
                                width_factor=width_factor,
                                return_idx=True) 
        result = self.ramp_molecular_erfc(v_drop, rate, polish=polish) 
        if result != 0:
            return 1            

        #on final polishing run, check if bolometric heating/
        #cooling should be added back in
        if polish:
            self.windsoln.add_user_vars(expedite=True)
            heat = self.windsoln.soln['heat_ion']
            cool = self.windsoln.soln['cool_PdV']
#                 heat,cool = self.quick_calc_heat_cool()
            if len(np.where(-cool[:20]<heat[:20])[0]) > 3:
                out_str = f"  NOTE: Photoionization heating dominates down to base of sim ({self.windsoln.Rmin:.3f} Rp).\n       | Molecular layer still turned off. Max error in dM/dt ~ 10%. \n       | See documentation for workaround."
                print(out_str)
                return 0
            
            if (self.windsoln.bolo_heat_cool == 0) and (bolo_on==True):
                return self.turn_on_bolo()
            
            return 0
            # if (self.windsoln.molec_adjust <= 0) and (molec_adjust_on == True):
            #     # print("triggered in molec_on")
            #     return self.turn_on_molec_adjust(_called_indep=False)
        else:
            return 0
        
    def ramp_molecular_erfc(self,v_drop,rate,polish=False):
        '''Ramps the molecular-to-atomic transition complementary error function to the desired radial drop off location and rate. These values can be obtained from sim.erf_velocity().

        Should rarely need to be run by user. Manually changing the erfc parameters is only necessary if sim.converge_mol_atomic_transition() fails to converge the erfc parameters automatically.

        Args:
            v_drop (float): Velocity at the radial location of the erfc drop in units of cm/s
            rate (float): Rate of the erfc drop off

        Returns:
            int: Status code. 0 for success, 1 for failure.
        '''
        diffs = np.array([abs(v_drop-self.windsoln.erf_drop[0])/v_drop,
                         abs(rate-self.windsoln.erf_drop[1])/rate])
        if any(diffs) > 1e-5:
            current_erfs =  np.array([self.windsoln.erf_drop[0],self.windsoln.erf_drop[1]])
            goal_erfs    = np.array([v_drop,rate])
            try_smoothing = True
            fail = 0
            smaller_delta=np.array([0,0])
            # if (self.windsoln.bolo_heat_cool == 1):
            if (self.skip == False) or (polish == True):
                while np.mean(abs(goal_erfs - current_erfs)/goal_erfs) > 1e-5:
                    self._raster_print(f"   ...Ramping molecular-to-atomic transition erfc: Current {current_erfs[0]:.3e}, {current_erfs[1]:.3e}. Goal {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}.")
                    delta = goal_erfs - current_erfs
                    if fail > 0: #take smaller stepsizes to goal
                        if abs(sum(delta)) < abs(sum(smaller_delta)):
                            smaller_delta = delta
                        self._raster_print(f"   ...Proceeding with smaller stepsizes: Current {current_erfs[0]:.3e}, {current_erfs[1]:.3e}. Trying: {current_erfs[0]+smaller_delta[0]:.3e}, {current_erfs[1]+smaller_delta[1]:.3e}. Goal {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}.")
                        self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                current_erfs+smaller_delta)
                    else:
                        self._raster_print(f"   ...Trying: {current_erfs[0]+delta[0]:.3e}, {current_erfs[1]+delta[1]:.3e}. Current {current_erfs[0]:.3e}, {current_erfs[1]:.3e}. Goal {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}.")
                        self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                current_erfs+delta)
                    
                    fail = 0
                    while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                        fail += 1
                        limit = 10
                        if self.windsoln.nspecies > 5:
                            limit = 4
                        if fail >= limit:
                            current_erfs = self.windsoln.erf_drop
                            self._normal_print(f"  Too many attempts to ramp error function. Skipping subsequent attempts until polishing. Goal: {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}. Current: {current_erfs[0]:.3e}, {current_erfs[1]:.3e}.")
                            self.skip = True
                            return 1
                        elif (fail <= 2) & (try_smoothing==True): #try smoothing out erf
                            step_erfs = np.copy(current_erfs)
                            step_erfs[1] *= 5*fail
                            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                    step_erfs)
                            self._raster_print(f"  Erfc Fail {fail:d}: Smoothing transition: {step_erfs[1]}")
                            if fail == 2:
                                try_smoothing=False #if this doesn't work, don't keep trying

                        elif (fail<limit) : #try taking a smaller step
                            smaller_delta = delta/(fail+1)
                            self._raster_print(f"   ..Fail {fail}: Current {current_erfs[0]:.3e}, {current_erfs[1]:.3e}. Trying {current_erfs[0]+smaller_delta[0]:.3e}, {current_erfs[1]+smaller_delta[1]:.3e}. Goal {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}."+self.clear)
                            step_erfs = current_erfs+smaller_delta
                            # self._raster_print(f" Erfc Fail {fail:d}: goal {goal_erfs[0]:.3e} {goal_erfs[1]:.3e}, step {step_erfs[0]:.3e} {step_erfs[1]:.3e}")
                            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                    step_erfs)
                    current_erfs = self.windsoln.erf_drop
                # self._raster_print(f"  Successfully ramped molecular-to-atomic transition erfc to {self.windsoln.erf_drop[0]:.3e}, {self.windsoln.erf_drop[1]:.3e}."+self.clear)
                return 0
            else:
                # self._raster_print(f"  Skipping ramping of molecular-to-atomic transition erfc until polishing.")
                # self.skip = True
                return 1

        
    def ramp_base_bcs(self,user_override_press=False, base_press=1,static_bcs=None,
                      Kappa_opt=4e-3,Kappa_IR=1e-2,molec_adjust=None,
                      adiabat=False, polish=False,
                      intermediate=False, tolerance=0.1):
        """
        Ramps Rmin, mass density at Rmin (rho_rmin), and temperature at Rmin (T_rmin) to the physical values computed in base_bcs().
        Default is the values at base_press = 1 microbar. If the loaded solution base BC is at higher pressure, base_press will stay that value in subsequent iterations.

        Args:
            base_press (int or float, optional): Desired pressure at base of simulation in microbars. Defaults to 1. Default will be overridden if loaded solution has different base press, unless user_override_press=True.
            user_override_press (bool, optional): To lower or raise pressure of base, set True. Defaults to False.
            static_bcs (bool, optional): If True, skips running ramp_base_bcs() and maintains the current boundary conditions. Defaults to None (i.e., the value of self.static_bcs).
            Kappa_opt (float, optional): Optical opacity. Sets bolometric heating/cooling in molecular region below wind. Defaults to 4e-3.
            Kappa_IR (float, optional): IR opacity. Defaults to 1e-2.
            molec_adjust (float, optional): Dimensionless mean molecular weight. Mu of H2 is 2.3*mH. Defaults to self.windsoln.molec_adjust if =None. If molec_adjust <= 0, turns off molecular layer.
            adiabat (bool, optional): If True, computes base BCs assuming atmosphere is adiabatic below wind. Defaults to False.
            polish (bool, optional): If True, skips ramping T when self.windsoln.bolo_heat_cool=0 to avoid costly iteration. Defaults to False.
            intermediate (bool, optional): If True, returns early if average difference is below tolerance. Defaults to False.
            tolerance (float, optional): Tolerance for intermediate ramping. Defaults to 0.1.

        Returns:
            int: 0 if successfully ramped base boundary conditions, 1 if failed.
        """
        if static_bcs is None:
            static_bcs = self.static_bcs
        else:
            self.static_bcs = static_bcs
        if static_bcs==True:
            return
        rho_scale = self.windsoln.scales_dict['rho']
        T_scale = self.windsoln.scales_dict['T']
        #Unless user wants to override the base pressure, take the pressure to be
        #the base pressure saved to solution
        if user_override_press == False:
            P = self.windsoln.soln['rho'][0]*const.kB*self.windsoln.soln['T'][0]/(self.windsoln.calc_mu()[0])
            rounded_base = np.round(P/10)*10
            if rounded_base < 1:
                rounded_base = 1
            if rounded_base != base_press:
                base_press = rounded_base #=P would prob be more accurate in case someone sets it to 2
        
        if molec_adjust is not None:
            self._raster_print(f"  Setting molec_adjust to {molec_adjust:.2f}")
            self.ramp_molec_adjust(molec_adjust,_called_indep=False,converge_final=False)
        else:
            molec_adjust = self.windsoln.molec_adjust

        goal_bcs = self.base_bcs(None,
                                 Kappa_opt,Kappa_IR,adiabat,base_press,
                                 user_override_press)
        goal_bcs = np.array(goal_bcs)
        if self.windsoln.bolo_heat_cool == 0: #an approximation
            goal_bcs[3] = np.average(self.windsoln.soln['T'][:10])/self.windsoln.scales_dict['T']

        curr_bcs = np.array(self.windsoln.bcs_tuple[:4])        
        avg_diff = np.array((abs(curr_bcs-goal_bcs)/goal_bcs)[[0,2,3]])
        if (intermediate==True) and all(avg_diff < tolerance):
            return
        
        if sum(abs((goal_bcs - curr_bcs)/goal_bcs)[[0,2,3]]) > 1e-2:
            self._raster_print(f"  Attempting to ramp Rmin:{curr_bcs[0]:.2f}->{goal_bcs[0]:.2f}Rp, rho:{curr_bcs[2]*rho_scale:.3e}->{goal_bcs[2]*rho_scale:.3e}g/cm3, T:{curr_bcs[3]*T_scale:.0f}->{goal_bcs[3]*T_scale:.0f}K")
            self.inputs.write_bcs(*goal_bcs,*self.windsoln.bcs_tuple[4:])

            if self.run_wind(expedite=True,calc_postfacto=False) == 1: 
                self._raster_print(f"\n  ..Initial jump failed. Ramping variables individually."+self.clear)
                if abs((self.windsoln.Rmin - goal_bcs[0])/goal_bcs[0]) > 4e-3:
                    self.ramp_Rmin(goal_bcs[0],_called_in_polish=True)
                if abs((self.windsoln.T_rmin  - goal_bcs[3])/goal_bcs[3]) > 9e-2:
                    # if (polish == True) & (self.windsoln.bolo_heat_cool==0):
                    #     print(f"   Skipping ramping T to avoid costly iterations. Current: {self.windsoln.T_rmin*T_scale:.0f}K. Estimated goal: {goal_bcs[3]*T_scale:.0f}K")
                    # else:
                    self.ramp_T_rmin(goal_bcs[3],_called_in_polish=True)
                if abs((self.windsoln.rho_rmin -  goal_bcs[2])/goal_bcs[2]) > 1e-4:
                    self.ramp_rho_rmin(goal_bcs[2],_called_in_polish=True)

        curr_bcs = np.array(self.windsoln.bcs_tuple[0:4])
        if sum(abs((goal_bcs - curr_bcs)/goal_bcs)[[0,2,3]]) > tolerance*1.01: 
            if (polish == True) & (self.windsoln.bolo_heat_cool==0):
                # self._normal_print("Successfully ramped base boundary conditions.")    
                self._raster_print("  Successfully ramped base boundary conditions.")          
                return 0
            else:
                self._raster_print("  Failed to ramp bcs. This can occur when attempting to set base at too high or too low of pressure (e.g., for very high or very low gravity planets).")
                self._normal_print(f"  Rmin: Current {curr_bcs[0]:.2f}, Goal {goal_bcs[0]:.2f}Rp; Rho: Curr {curr_bcs[2]*rho_scale:.3e}, Goal {goal_bcs[2]*rho_scale:.3e}g/cm3; T: Curr {curr_bcs[3]*T_scale:.0f}, Goal {goal_bcs[3]*T_scale:.0f}K")
                return 1
        else:
            self._raster_print("  Successfully ramped base boundary conditions.")          
            # self._normal_print("Successfully ramped base boundary conditions.")           
            return 0

    
    
    def base_bcs(self,molec_adjust=None,Kappa_opt=4e-3,
                 Kappa_IR=1e-2,adiabat=False,
                 base_press=1,user_override_press=False): 
        """
        Sets base of simulation below XUV tau=1 surface and assumes bolometric heating & cooling of molecules dominates there.

        Computes the density and temperature at either R_ubar (microbar pressure radius) or R_IR (vertical tau=1 surface to IR radiation). Between this radius and the nanobar radius where wind is launched, the temperature profile is an isotherm at the skin temperature. If R_IR is below Rp, it takes R_mbar to be the base of the simulation.

        Assumes that Rp is the slant path optical tau = 1 surface. Computes skin temperature by balancing bolometric heating and cooling. Computes effective temperature from stellar bolometric luminosity. For highly irradiated planets (most planets for which atmospheric escape is significant), it is an isotherm at T_eff=L_star/(16*pi*sigma_sb*a^2) between R_IR and Rp or an isotherm at T_skin between R_mbar and Rp. If adiabat=True, computes the R_IR assuming an adiabat between T_eff and T_skin.

        For high metallicity atmospheres, users should increase the molecular adjustment factor (molec_adjust). Likewise, they should change Kappa_opt and Kappa_IR.

        Args:
            molec_adjust (float, optional): Dimensionless mean molecular weight. Mu of H2 is 2.3*mH. Defaults to self.windsoln.molec_adjust if None. If molec_adjust == 0, turns off molecular layer.
            Kappa_opt (float, optional): Optical opacity. Defaults to 4e-3.
            Kappa_IR (float, optional): IR opacity. Defaults to 1e-2.
            adiabat (bool, optional): If True, computes IR BCs assuming an adiabat. Defaults to False.
            base_press (int or float, optional): Pressure at base of simulation in microbars. Defaults to 1.
            user_override_press (bool, optional): If True, overrides base pressure. Defaults to False.

        Returns:
            tuple: (R, Rmax, rho, T)
                R (float): Radius at base of simulation in units of Rp
                Rmax (float): Set by original boundary conditions, units of Rp
                rho (float): Density in units of RHO0 from defs.h (1e-15 ~nBAR*MH/K/T0)
                T (float): Effective temperature in units of T0 (1e4 K)
        If R_IR > Rp:
            Returns: R_IR, Rmax, rho_IR, T_eff
        Else:
            Returns: Rp, Rmax, rho_Rp, T_skin
        """
#         #computing current pressure at base. If current base pressure > 1, set current as base pressure
#         #this is because some high flux planets may need lower lower BC to capture all heating
        rho_scale = self.windsoln.scales_dict['rho']
        T_scale = self.windsoln.scales_dict['T']
        # mu = self.windsoln.calc_mu()[0]
        if (molec_adjust is None):
            molec_adjust = self.windsoln.molec_adjust
            mu = self.windsoln.calc_mu()[0]
        else:
            mu = molec_adjust*const.mH
        if user_override_press == False:
            P = self.windsoln.soln['rho'][0]*const.kB*self.windsoln.soln['T'][0]/mu
            rounded_base = np.round(P/10)*10
            if rounded_base < 1:
                rounded_base = 1
            if rounded_base != base_press:
                base_press = rounded_base #=P would prob be more accurate in case someone sets it to 2
        
        Rp = self.windsoln.Rp
        # if 0 < molec_adjust <= 1:
        #     self._normal_print("WARNING: Molecular adjustment factor should be >1 to account for increased mean weight due to molecules, instead of atoms, lower in atmosphere (below the wind).")
        F_opt = self.windsoln.Lstar/(4*np.pi*self.windsoln.semimajor**2)
        T_skin = (F_opt*(Kappa_opt+Kappa_IR/4)/(2*const.sig_SB*Kappa_IR))**0.25
        T_eff = (F_opt/(4*const.sig_SB))**0.25
        #rho_Rp derived from optical slant path tau=1 geometry
        rho_Rp = np.sqrt(mu*const.G*self.windsoln.Mp / (8 * Rp**3 *const.kB*T_skin)) / Kappa_opt
        P_Rp = rho_Rp*const.kB*T_skin / mu #in barye = 1e-6 bar

        #take solution base as base_press to ease computational burden of going down to Rp
        cs2 = (const.kB*T_skin/mu)
        R_mbar = (((cs2/(const.G*self.windsoln.Mp))* np.log(base_press/P_Rp) 
                   + 1/self.windsoln.Rp )**(-1))
        rho_mbar = base_press * mu / (const.kB*T_skin)

        if adiabat == True:
            G = const.G
            kB = const.kB
            denom = (7/2)*(kB*T_eff/mu) * ( 1 - (2*Kappa_opt/Kappa_IR + 1/2)**0.25 ) + G*self.windsoln.Mp/Rp
            R_IR = G*self.windsoln.Mp / denom
            Kappa_term =  (2*Kappa_opt/Kappa_IR+1/2)**(-6/8)/Kappa_opt
            rho_IR = np.sqrt( mu*G*self.windsoln.Mp / (8*Rp**3 * kB*T_eff) ) * Kappa_term
        else:
            Hsc = const.kB*T_eff*(Rp**2) /(mu*const.G*self.windsoln.Mp) #scaleheight at Rp (approximately Hsc at R_IR)
            rho_IR = 1/(Kappa_IR*Hsc)
            # R_IR = rho_Rp*np.exp(-Rp/Hsc)
            if Rp/Hsc > 690:
                R_IR = 0
            else:
                R_IR = Rp**2 /( Hsc * np.log(rho_IR/rho_Rp*np.exp(Rp/Hsc)) )
        if R_IR > Rp:
            self._normal_print('Using vertical tau=1 for IR photons, R_IR.')
            return R_IR/Rp,self.windsoln.Rmax,rho_IR/rho_scale,T_eff/T_scale
        else:
            if self.windsoln.bolo_heat_cool == 0: #an approximation
                T_skin = np.average(self.windsoln.soln['T'][:10])

            return R_mbar/Rp,self.windsoln.Rmax,rho_mbar/rho_scale,T_skin/T_scale

        
        
    def converge_Rmax(self,final_converge_Ncol=True):
        """
        Self-consistently sets Rmax to be the Coriolis length. Does not worry about converging other boundaries. Cannot be expedited as the Coriolis length must be calculated.

        Args:
            final_converge_Ncol (bool, optional): If True, converges Ncol_sp self-consistently again after Rmax has been set to Rcori. Defaults to True.

        Returns:
            int: Status code. 0 for success, other values for failure modes.
        """
        try: # Check if r_cori has be calculated
            self.windsoln.R_cori
        except AttributeError:
            if self.windsoln.integrate_outward == 0:
                self.inputs.write_flags(*self.windsoln.flags_tuple,
                                        integrate_out=True)
                out = self.run_wind(calc_postfacto=False)
                attempt = 0
                while out == 4:
                    attempt += 1
                    self.windsoln.bcs_tuple[-2] = self.windsoln.Ncol_sp*10*attempt
                    self._raster_print(f"  ..Raising Ncol to allow outward integration. Attempt {attempt:d}.")
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                    if attempt > 6:
                        self._normal_print("Failed to integrate outwards even after raising Ncol. Suspicious.")
                        self.converge_Ncol_sp()
                        return 1
                    out = self.run_wind(calc_postfacto=False)
                else:
                    if out == 1:
                        self._normal_print("Relaxation error.")
                        return 1
                    elif out == 0:
                        self._raster_print("Successfully integrated outwards. Now converging Rmax...")
                        self.windsoln.calc_Coriolis()
            else:
                self.windsoln.calc_Coriolis()
    
        # First while statement extends domain far past true Coriolis length
        while (self.windsoln.Rmax < self.windsoln.R_cori):
            self._raster_print("  ..Rmax {:.4e}, r_Cori {:.4e}"
                  .format(self.windsoln.Rmax, self.windsoln.R_cori))
            bcs_tuple = self.windsoln.bcs_tuple
            # Just go far enough to ensure true R_cori is within domain
            bcs_tuple[1] = min(self.windsoln.R_cori+10,
                               1.5*self.windsoln.R_cori)
            self.inputs.write_bcs(*bcs_tuple)
            failed = 0
            if self.run_wind(calc_postfacto=False) != 0:
                self._normal_print("\nERROR: Failed to integrate outwards to r_Cori: {:.3e} "
                      "reverting back to {:.3e}."
                      .format(self.windsoln.R_cori, self.windsoln.Rmax))
                bcs_tuple[1] = self.windsoln.Rmax
                self.inputs.write_bcs(*bcs_tuple)
                return 1
            self.windsoln.calc_Coriolis()
            
        # Second while statement whittles down to true Coriolis length
        while (abs(1.-self.windsoln.Rmax/self.windsoln.R_cori) > 1e-2):
            self._raster_print("  ..Rmax {:.4e}, r_Cori {:.4e}"
                  .format(self.windsoln.Rmax, self.windsoln.R_cori))
            delta_Rmax = self.windsoln.Rmax-self.windsoln.R_cori
            bcs_tuple = self.windsoln.bcs_tuple
            bcs_tuple[1] = self.windsoln.Rmax-0.99*delta_Rmax
            self.inputs.write_bcs(*bcs_tuple)
            failed = self.run_wind(calc_postfacto=False)
            if failed:
                self._normal_print('\nERROR: Failed to even get relax code to run')
                return 1
            self.windsoln.calc_Coriolis()
        if final_converge_Ncol == True:
            result = self.converge_Ncol_sp(expedite=False,integrate_out=True)
            if result == 0:
                self._normal_print(f"Successfully converged Rmax to {self.windsoln.Rmax:.6f} (Ncol also converged).")
            if result == 3:
                self._normal_print(f"Successfully converged Rmax to {self.windsoln.Rmax:.6f} (To preserve outward integration, Ncol was slightly increased => Mdot error < ~10%).")
        else:
            self._normal_print(f"Successfully converged Rmax to {self.windsoln.Rmax:.6f}. It is recommended to run converge_Ncol_sp(integrate_out=True) as Ncol will have updated.")
        self.windsoln.add_user_vars()
        # self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
        self.integrate_out(quiet=True)
        return 0
    

    def converge_Ncol_sp(self,expedite=False,integrate_out=False,quiet=False):
        """
        Description:
            Self-consistently converges Ncol at the sonic point, such that the 
            column density boundary condition for each species, Ncol_sp,
            matches the number density at the sonic point.

        Args:
            expedite - bool; default = False. True speeds solutions by not integrating
                                              past Sonic Point.
            integrate_out - bool; default = False. Forces outward integration by raising Ncol_sp
                                                   incrementally. Final Ncol_sp may be slightly higher 
                                                   than self-consistent.
            quiet - bool; default = False. Prints warnings about Ncol_sp being only estimates if
                                            integrate_out = False.
        """
        scale = self.windsoln.scales_dict['Ncol_HI'] #scales the same across species
        #self-consistent Ncol finder
        Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False)
        Ncol_goal[Ncol_goal==0] = 1e-10*Ncol_goal[0]
        Ncol_current = np.copy(self.windsoln.Ncol_sp)
        
        if expedite == True:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)
        else:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
                
            
        #needs large convergence condition, otherwise will iterate indefinitely
        niter = 0 #a check
        while np.abs(np.mean((Ncol_goal-Ncol_current)/Ncol_goal)) > 0.08:
            niter += 1
            if niter > 10:
                self._normal_print(f'Too many iterations (current average diff {avg:.2e}).')
                return 1
            avg = np.abs(np.mean((Ncol_goal-Ncol_current)/Ncol_goal))
            self._raster_print(f'  ...Ramping Ncol Iter {niter:d}: Current average diff {avg:.2e}')            
            self.windsoln.bcs_tuple[5][:] = Ncol_goal  
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            fail = 1
            result = self.run_wind(expedite=expedite,calc_postfacto=False)
                    
            while result != 0:
                attempt = 0
                if result == 4:
                    if integrate_out == False:
                        self._raster_print(" ...Integration error. Turning off outward integration.")
                        self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)
                        result = self.run_wind(expedite=expedite,calc_postfacto=False) 
                    else:
                        self._raster_print("  ...Integration error. To preserve outward integration, raising Ncol_sp slightly and returning.")
                        self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
                        self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                             Ncol_goal*1.1,
                                             self.windsoln.bcs_tuple[:-1])
                        fail = 0 
                        result = self.run_wind(expedite=False,calc_postfacto=False)
                        while result==4:
                            fail+=1
                            self._raster_print(f"  ...Attempting to integrate out. Ncol_sp x{(1+(fail+1)/10):.1f}")
                            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                                                  Ncol_goal*1+(fail+1)/10,
                                                  self.windsoln.bcs_tuple[:-1])
                            result = self.run_wind(expedite=False,calc_postfacto=False)
                        if result == 0:
                            if fail == 0:
                                return 0  
                            if fail > 0:
                                self._raster_print(f"  Successfully integrated out w/ Ncol_sp = {(1+(fail)/10):.1f} x goal Ncol_sp. Error in Mdot < ~10%.\n")
                                return 3
                else:
                    delta = (Ncol_goal - Ncol_current)/(2*fail)
                    Ncol_step = Ncol_current + delta
                    self.windsoln.bcs_tuple[5][:] = Ncol_step 
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                    self._raster_print(f'  ...Fail {fail:d}: Attempting Ncol_sp = {Ncol_step*scale}.')

                    if fail > 10:
                        self._normal_print(f"\nFailed at Ncol_sp = {Ncol_current} & Goal = {Ncol_goal}.") #summed neutral number density
                        return 2
                    result = self.run_wind(expedite=expedite,calc_postfacto=False)
                    fail += 1
            
            Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False) 
            Ncol_current = np.copy(self.windsoln.Ncol_sp)

        Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False) 
        avg = np.abs(np.mean((Ncol_goal-Ncol_current)/Ncol_goal))
        if quiet == False:
            self._raster_print(f' Success! Average Ncol difference: {avg:.2e}')

        #Final convergence of all species at once (MAY NOT WANT THIS BECAUSE WILL MOVE GOAL again)
        if niter != 0:
            if expedite == True:
                if quiet == False:
                    self._raster_print('  ...Attempting Final Ncol Convergence. (Note: This is an estimate. Cannot converge precisely without outward integration).')  
                    integrate_out=0
            if expedite == False:
                if quiet == False:
                    self._raster_print('  ...Attempting Final Ncol Convergence.')  

            self.windsoln.bcs_tuple[5][:] = Ncol_goal
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=integrate_out)
            result = self.run_wind(expedite=expedite,calc_postfacto=False)
            if result == 4:
                if quiet == False:
                    self._normal_print(f'Failure at last Ncol convergence. Last working solution should be fine for most purposes. Average diff from self-consistent value: {avg:.2e}')
                return 4
            elif result == 1:
                self._normal_print('Failure at the endgame.')
                return 1            
        return 0

    
    
    def self_consistent_Ncol(self,method=1,warning=True):
        """
        Computes the self-consistent column density sonic point boundary condition from the neutral number density of a given species.

        Args:
            method (int, optional): Method for calculation. Defaults to 1.
            warning (bool, optional): If True, prints warnings. Defaults to True.

        Returns:
            tuple: (currents, goals)
                currents (np.ndarray): Array of current Ncol_sp values.
                goals (np.ndarray): Array of self-consistent Ncol_sp values.
        """
        #sonic_index = index of sonic point radius
        scale = self.windsoln.scales_dict['Ncol_HI'] #scales the same across species
        goals = np.zeros(self.windsoln.nspecies)
        currents = np.zeros(self.windsoln.nspecies)
        sonic_index = np.searchsorted(self.windsoln.soln_norm['r'],
                                      self.windsoln.R_sp)
        for j,sp_name in enumerate(self.windsoln.species_list):
            sp_name = sp_name.replace(' ','') 
            n_neutral = self.windsoln.HX[j]*self.windsoln.soln['rho'] 
            self.windsoln.species_list
            n_neutral *= self.windsoln.soln['Ys_'+sp_name]
            n_neutral /= self.windsoln.atomic_masses[j] #num dens of neutral atoms of species
            nn = n_neutral*np.diff(self.windsoln.soln['r'],prepend=self.windsoln.Rp)
            #goal self-consistent boundary cond.
            goals[j] = np.flip(np.cumsum(np.flip(nn)))[sonic_index]/scale
            #boundary condition in input
            if method == 2:
                currents[j] = self.windsoln.soln['Ncol_'+sp_name][sonic_index]/scale
        if method == 1: 
            currents = self.windsoln.Ncol_sp
        if self.windsoln.flags_tuple[-1] == 0:
            if warning:
                self._normal_print("This is an estimate of Ncol. Cannot converge precisely without outward integration.") 
            return currents, goals#*8 #to account for the fact that we can't cumsum n out to Rmax
        if self.windsoln.flags_tuple[-1] == 1:
            return currents, goals 
        
        
        
#     def converge_mol_atomic_transition(self,width_factor=0.,return_idx=False,_called_in_ramp_bcs=False,polish=False):
#         """
#         Ramps the complementary error function that governs the transition from the molecular to atomic regions in the atmosphere.

#         Criterion for transition: when photoionization heating begins to dominate over the PdV cooling and a wind will launch.
#         Below the wind, for an average planet, molecules have not photodissociated so mu should be mean molecular weight instead of mean atomic weight and the erf enforces this. Molecular opacities mean that bolometric heating and cooling dominate the energy budget and create an isotherm in the molecular region below the wind. In the optically thin atomic wind, bolometric heating and cooling are negligible, so the erf also enforces the drop off of bolometric heating and cooling.

#         If a user needs to call this function for some reason, width_factor is the only argument they should need to change. For numerical or physical reasons it may be necessary for the transition to occur over more than 1 scaleheight. In this case, increase width_factor.

#         Args:
#             width_factor (float, optional): Factor to adjust the width in scaleheights of the transition region. Defaults to 0.
#             return_idx (bool, optional): If True, returns radial index of transition point. Defaults to False.
#             _called_in_ramp_bcs (bool, optional): Set True if called within ramp_base_bcs to avoid recursion issues. Defaults to False.
#             polish (bool, optional): If True, computes full postfacto heat_ion & cool_PdV for more accurate transition point. Defaults to False.

#         Returns:
#             int or tuple: Index of transition if return_idx is True, otherwise status code.
#         """
#         if (self.windsoln.bolo_heat_cool == 0) and (polish==False):
#             return
#         if (polish==False) and (self.skip==True):
#             return
        
#         #loading last working solution to do this
#         self.load_planet(self.path+'saves/windsoln.csv',calc_postfacto=False,
#                          print_atmo=False,print_warnings=False)
#         v_drop,rate = self.erf_velocity(_called_in_ramp_bcs=_called_in_ramp_bcs,
#                                         polish=polish,
#                                         width_factor=width_factor) 
       
#         diffs = np.array([abs(v_drop-self.windsoln.erf_drop[0])/v_drop,
#                          abs(rate-self.windsoln.erf_drop[1])/rate])
#         if any(diffs) > 1e-5:
#             current_erfs =  np.array([self.windsoln.erf_drop[0],self.windsoln.erf_drop[1]])
#             goal_erfs    = np.array([v_drop,rate])
#             try_smoothing = True
#             fail = 0
#             smaller_delta=np.array([0,0])
#             if (self.windsoln.bolo_heat_cool == 1):
#                 if (self.skip == False) or (polish == True):
#                     while np.mean(abs(goal_erfs - current_erfs)/goal_erfs) > 1e-5:
#                         delta = goal_erfs - current_erfs
#                         if fail > 0: #take smaller stepsizes to goal
#                             if abs(sum(delta)) < abs(sum(smaller_delta)):
#                                 smaller_delta = delta
#                             self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
#                                                   current_erfs+smaller_delta)
#                         else:
#                             self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
#                                                   current_erfs+delta)
                        
#                         fail = 0
#                         while self.run_wind(expedite=True,calc_postfacto=False) != 0:
#                             fail += 1
#                             limit = 10
#                             if self.windsoln.nspecies > 5:
#                                 limit = 4
#                             if fail >= limit:
#                                 current_erfs = self.windsoln.erf_drop
#                                 self._normal_print(f"Too many attempts to ramp error function. Skipping subsequent attempts until polishing. Goal: {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}.",
#                                       f" Current: {current_erfs[0]:.3e}, {current_erfs[1]:.3e}.")
#                                 self.skip = True
#                                 return 1#self.turn_off_bolo()
#                             elif (fail <= 2) & (try_smoothing==True): #try smoothing out erf
#                                 step_erfs = np.copy(current_erfs)
#                                 step_erfs[1] *= 5*fail
#                                 self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
#                                                       step_erfs)
#                                 self._raster_print(f" Erfc Fail {fail:d}: Smoothing transition: {step_erfs[1]}")
#                                 if fail == 2:
#                                     try_smoothing=False #if this doesn't work, don't keep trying

#                             elif (fail<limit) : #try taking a smaller step
#                                 smaller_delta = delta/(fail+1)
#                                 step_erfs = current_erfs+smaller_delta
#                                 self._raster_print(f" Erfc Fail {fail:d}: goal {goal_erfs[0]:.3e} {goal_erfs[1]:.3e}, step {step_erfs[0]:.3e} {step_erfs[1]:.3e}")
#                                 self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
#                                                       step_erfs)
#                         current_erfs = self.windsoln.erf_drop
#                     return 0
#                 else:
#                     return 1
                            
#         #on final polishing run of isotherm, check if bolometric heating/
#         #cooling should be added back in
#         if (self.windsoln.bolo_heat_cool == 0) and (polish==True):
#             self.windsoln.add_user_vars(expedite=True)
#             heat = self.windsoln.soln['heat_ion']
#             cool = self.windsoln.soln['cool_PdV']
# #                 heat,cool = self.quick_calc_heat_cool()
#             if len(np.where(-cool[:20]<heat[:20])[0]) > 3:
#                 out_str = f"  NOTE: Photoionization heating dominates down to base of sim ({self.windsoln.Rmin:.3f} Rp).\n       | Bolometric heat/cooling still turned off. Max error in dM/dt ~ 10%. \n       | See documentation for workaround."
#                 print(out_str)
#                 return 1
#             else: #if appropriate, ramp back in bolometric heating/cooling
#                 while self.windsoln.bolo_heat_cool < 1:
#                     flags = self.windsoln.flags_tuple
#                     delta = 1-self.windsoln.bolo_heat_cool
#                     flags[2] += delta
#                     self.inputs.write_flags(*flags)
#                     fail = 0
#                     while self.run_wind(calc_postfacto=False) != 0:
#                         fail+=1
#                         flags[2] = self.windsoln.bolo_heat_cool + 0.1/fail
#                         self._raster_print(f'Ramping back in bolometric heating & cooling. Trying factor of {flags[2]:.2f}')
#                         self.inputs.write_flags(*flags)
#                         if fail>10:
#                             self._normal_print("Warning: Bolometric heating/cooling failed to ramp back in.")
#                             energy_plot(self.windsoln)
#                             return 1
#                 self._normal_print("  Bolometric heating/cooling successfully ramped back in.")
#                 return 0
                            
# #             else:
# #                 return 1
#         else:
#             return 0
        
        
    def _quick_calc_heat_cool(self):
        ''' Inexpensively calculates the photoionization heating and PdV cooling as a function of radius. Used to compute the arguments in the complementary error function that enforces the molecular to atomic transition.

            Args: 
                None
            Returns:
                heat (array): photoionization heating rate as a function of radius
                cool (array): PdV cooling rate as a function of radius
        '''
        n_tot = np.zeros_like(self.windsoln.soln['rho'])
        
        unspaced_list = [sp.replace(' ','') for sp in McAtom.formatting_species_list(self.windsoln.species_list)]
        for j,species in enumerate(unspaced_list):
            n = self.windsoln.HX[j]*self.windsoln.soln['rho']
            n /= self.windsoln.atomic_masses[j]
            n_tot += n
            self.windsoln.soln['n_'+species] = n
            
        n_HII = (1-self.windsoln.soln['Ys_HI'])*self.windsoln.soln['n_HI']

        # Multifrequency calculations
        background_ioniz_frac = n_HII/n_tot #only based on H fraction
        background_ioniz_frac[background_ioniz_frac<0] = 1e-10

        #Re-written for speed, so most is reshaping 
        heating_rate = np.zeros((len(self.windsoln.soln), self.windsoln.nspecies))
        total_heating = np.zeros(len(self.windsoln.soln))

        Ncol_arr    = self.windsoln.soln.iloc[:,4+self.windsoln.nspecies:4+2*self.windsoln.nspecies]
        Ncol_arr[Ncol_arr<0] = 0      #less than 0 because unconverged soln returns negative Ncols. This feels sus.
        sigmas = self.windsoln.sim_spectrum.iloc[:,2:2+self.windsoln.nspecies]
        taus = np.dot(Ncol_arr,sigmas.T)
        wPhi = np.multiply(np.tile(self.windsoln.sim_spectrum['wPhi'],(len(taus),1)),np.exp(-taus))*self.windsoln.Ftot

        # Adapted from Mocassin (Shull & Steenberg 1985)
        # Accounts for secondary ionizations due to highly energetic (>100eV) incoming photons
        frac_in_heat = 0.9971 * (1 - pow(1-pow(background_ioniz_frac,0.2663),1.3163))
    #     self._normal_print(frac_in_heat)
        for s, species in enumerate(self.windsoln.species): #48s
            E_matrix = np.tile(self.windsoln.E_wl,(len(background_ioniz_frac),1))
            Ncol_matrix = np.tile(Ncol_arr.iloc[:,s],(self.windsoln.npts,1)).T
            sig_matrix = np.tile(sigmas.iloc[:,s],(len(background_ioniz_frac),1))

            taus_temp = taus
            taus_temp[taus==0] = np.inf
            f = np.nan_to_num(sig_matrix*Ncol_matrix / taus_temp) #frac of incoming photon energy that will interact with species s
            E0_matrix = E_matrix - self.windsoln.ion_pot[s]

            heatfrac_matrix = (np.tile(frac_in_heat,(len(self.windsoln.E_wl),1))).T
            heatfrac_matrix[E0_matrix<6.408707e-11] = 1 #zeroing where E_0 too low for secondary ionizations (100eV)

            heating_rate[:,s] = np.sum((heatfrac_matrix)*E0_matrix*sig_matrix*f*wPhi,axis=1)
            n_abs = self.windsoln.soln['n_'+species]
            heating_rate[:,s] *= n_abs #/rho temporary

        heat = np.sum(heating_rate,axis=1)

        P = self.windsoln.soln['rho']*const.kB*self.windsoln.soln['T']/(self.windsoln.molec_adjust*const.mH)
        cool = P*self.windsoln.soln['v']/self.windsoln.soln['rho']
        cool *= np.gradient(self.windsoln.soln['rho'], self.windsoln.soln['r'])

        return heat, cool
    
    
    
    def turn_off_bolo(self):
        """
        Turns off bolometric heating and cooling AND the mean molecular weight adjustment in the molecular region. These values are currently coupled by the error function defined in converge_mol_atomic_transition().

        Returns:
            int: 0 if successful, 1 if failed to turn off bolometric heating/cooling.
        """
        self.load_planet(self.path+'saves/windsoln.csv',calc_postfacto=False,
                         print_atmo=False,print_warnings=False)
#         if failed_bolo_turn_off == True:
#             self._normal_print("Previously failed to turn off bolometric heating and cooling, so not trying here.")
#             return
        if self.windsoln.bolo_heat_cool == 1:
            self._raster_print('  ...Turning off bolometric heating/cooling.')
        #turning off bolo_heat_cool
        while self.windsoln.bolo_heat_cool > 0:
            flags = self.windsoln.flags_tuple
            bolo_flag = np.copy(flags[2])
            delta = -flags[2]
            flags[2] = 0 #turning off bolo heating and cooling
            self.inputs.write_flags(*flags)
            fail = 0
            while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                fail+=1
                if fail>5:
                    self._normal_print(f"WARNING: Failed to turn off bolometric heating/cooling. This is unusual. Current multiplicative factor: {self.windsoln.bolo_heat_cool}. Goal: 0.")
                    self.windsoln.flags_tuple[2] = 1
                    self.inputs.write_flags(*self.windsoln.flags_tuple)
                    return 1
                delta /= 2
                flags[2] = bolo_flag + delta
                self._raster_print(f'Turning off bolometric heating & cooling. Trying factor of {flags[2]:.2f}')
                self.inputs.write_flags(*flags) 
        self._raster_print('  ...Successfully turned off bolometric heating/cooling (and molecular layer)\n')      
        return 0
    
    def turn_on_bolo(self):
        """
        Turns on bolometric heating and cooling in the region below the wind.

        Returns:
            int: 0 if successful, 1 if failed to turn on bolometric heating/cooling.
        """
        while self.windsoln.bolo_heat_cool < 1:
            flags = self.windsoln.flags_tuple
            delta = 1-self.windsoln.bolo_heat_cool
            flags[2] += delta
            self.inputs.write_flags(*flags)
            fail = 0
            while self.run_wind(calc_postfacto=False) != 0:
                fail+=1
                flags[2] = self.windsoln.bolo_heat_cool + 0.1/fail
                self._raster_print(f'Ramping back in bolometric heating & cooling. Trying factor of {flags[2]:.2f}')
                self.inputs.write_flags(*flags)
                if fail>10:
                    self._normal_print("Warning: Bolometric heating/cooling failed to ramp back in.")
                    energy_plot(self.windsoln)
                    return 1
        
        self._normal_print("  Bolometric heating/cooling successfully ramped back in.")
        return 0   
        
    def erf_velocity(self,return_idx=False,_called_in_ramp_bcs=False, 
                     polish=False,width_factor=0, called_in_mol_layer=False):
        """
        Defines the drop-off radius of the complementary error function that governs the drop-off of bolometric heating/cooling and the mean molecular weight in the isothermal part of the wind as photoionization heating begins to dominate and the atmosphere becomes atomic and non-isothermal.

        Args:
            return_idx (bool, optional): If True, returns radial index where drop-off occurs. Defaults to False.
            _called_in_ramp_bcs (bool, optional): Avoids recursion of ramping bases if called in ramp_base_bcs. Defaults to False.
            polish (bool, optional): If True, computes full bolometric heating/cooling and photoionization heating curves for improved accuracy. Defaults to False.
            width_factor (float, optional): Sets width in radius space over which the erfc drops off as a factor of scaleheight at the transition index (Hsc). If 0, uses current width of transition in Hsc. Defaults to 0.

        Returns:
            float: Velocity at location of drop-off in units of cm/s.
            float: Rate of drop-off.
            tuple (optional): If return_idx=True, also returns (width_factor, drop_idx).
        """
                    
        #now turns off boloheat/cool if sim base not deep enough
        if polish == False:
            heat, cool = self._quick_calc_heat_cool()
        if polish == True:
            self.windsoln.add_user_vars(expedite=True)
            heat = self.windsoln.soln['heat_ion']
            cool = self.windsoln.soln['cool_PdV']
        v = self.windsoln.soln_norm['v']
        r = self.windsoln.soln_norm['r']

        try:
            #finds where photoion heating first starts to dominate over PdV cool
            drop_index = np.where(heat>-cool)[0][0]
        except IndexError:
            drop_index = 0 
        
        #Drop-off rate / gradient of erf
        def rate_calc(idx,width_factor):
            #Approximate pressure scaleheight at drop radius in units of Rp
            Hsc =  const.kB*self.windsoln.soln['T'][idx]*self.windsoln.soln['r'][idx]
            Hsc /= (self.windsoln.calc_mu()[0]*const.G*self.windsoln.Mp)
            #approximates mu as molecular value
            #Compute the desired rate of erfc drop-off as 
            if idx >= 10:
                slope = (v[idx+10] - v[idx-10])/(r[idx+10] - r[idx-10])
            else:
                slope = (v[idx+10] - v[idx])/(r[idx+10] - r[idx])
            width = width_factor*Hsc 
            rate = slope*width
            return rate
            
        #if user does not specify a width_factor, will use the current approximate value
        current_rate = self.windsoln.erf_drop[1]
        Hsc =  const.kB*self.windsoln.soln['T'][drop_index]*self.windsoln.soln['r'][drop_index]
        Hsc /= (self.windsoln.calc_mu()[drop_index]*const.G*self.windsoln.Mp)
        if drop_index >= 10:
            slope = (v[drop_index+10] - v[drop_index-10])/(r[drop_index+10] - r[drop_index-10])
        else:
            slope = (v[drop_index+10] - v[drop_index])/(r[drop_index+10] - r[drop_index])
#         self._normal_print('drop index',drop_index)
        current_width = current_rate / (slope*Hsc)
        if width_factor == 0:
            width_factor = current_width
        if width_factor >= 15:
            width_factor = 15
        #If this is being called inside of the ramping function, do this 
        #to avoid recursion
        if (drop_index<=10) and (_called_in_ramp_bcs==True):
            if return_idx==False:
                return v[0],rate_calc(0,width_factor)
            else:
                return v[0],rate_calc(0,width_factor),width_factor,0
        if (drop_index<=10) and (called_in_mol_layer==True):
            if return_idx==False:
                return v[0],rate_calc(0,width_factor)
            else:
                return v[0],rate_calc(0,width_factor),width_factor,0
        elif (drop_index<=10):# and (self.failed_deeper_bcs_ramp==True):
            if (self.windsoln.molec_adjust > 0) and (self.try_turning_off):
                # print(self.windsoln.molec_adjust)
                self._raster_print("\nCurrent simulation Rmin nearly or fully inside of wind. Molecular layer will be turned off.\n")
                self.turn_off_bolo()
                self.try_turning_off=False
                if return_idx==False:
                    return v[0],rate_calc(0,width_factor)
                else:
                    return v[0],rate_calc(0,width_factor),width_factor,0
        rate = rate_calc(drop_index,width_factor)
        #If none of the above conditions are triggered, simply return v
        if return_idx == True:
            return v[drop_index],rate,width_factor,drop_index
        else:
            return v[drop_index],rate   
  
        
        
    def raise_Ncol_sp(self,to_total=0.8,by_factor=0,expedite=False):
        """
        Increase per-species column density at the sonic point. This is distinct from converging Ncol_sp. Too low of Ncol_sp may cause numerical instability, so if ramping is stuck, try raising Ncol_sp. Can also diagnose by plotting six_panel_plot—if Ncol_sp drops off sharply after the sonic point for several species, this may be the source of numerical instability.

        Args:
            to_total (float, optional): If by_factor=0, raises sum(Ncol_sp) until it equals to_total. Defaults to 0.8.
            by_factor (float, optional): If not 0, multiplies Ncol_sp by this factor. Defaults to 0.
            expedite (bool, optional): If True, does not integrate out. Defaults to False.

        Returns:
            int: Status code from run_wind.
        """
        if by_factor == 0:
            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                     to_total*(self.windsoln.Ncol_sp/sum(self.windsoln.Ncol_sp)), 
                     self.windsoln.bcs_tuple[-1])
            return self.run_wind(expedite=expedite)
        else:
            self.inputs.write_bcs(*self.windsoln.bcs_tuple[:-2],
                     by_factor*self.windsoln.Ncol_sp, 
                     self.windsoln.bcs_tuple[-1])
            return self.run_wind(expedite=expedite)


    def turn_off_tidal_grav(self):
        """
        Turns off the tidal gravity.

        Returns:
            None
        """
        flag = self.windsoln.flags_tuple[1]
        stepsize=0.1
        while np.round(flag,4) > 0:
            flag+= -stepsize    
            self._raster_print(f'  Trying {flag:.4f} x tidal grav term')
            self.inputs.write_flags(1,flag,1,0,integrate_out=False)
            if self.run_wind() == 0:
                continue
            else:
                flag+=stepsize #reset to last working
                stepsize/=10 
                self._raster_print(f'  Failed. Trying smaller stepsize {stepsize:.4f}.')
        self._normal_print('  Success! Tidal gravity turned off.')
        return

        
    def ramp_T_rmin(self, goal_T,integrate_out=False,_called_in_polish=False):
        """
        Ramps normalized temperature at Rmin.

        Args:
            goal_T (float): Goal temperature at Rmin in units of 1e4 K.
            integrate_out (bool, optional): If True, integrates out after ramping. Defaults to False.

        Returns:
            int: 0 if successful, 1 if failed to converge.
        """
        T_scale = self.windsoln.scales_dict['T']
        if goal_T > 1.1:
            self._normal_print("WARNING: T should be in units of 1e4 K. Returning...")
            return 
        while (abs(1.-self.windsoln.T_rmin/goal_T) > 1e-10):
            bcs_tuple = self.windsoln.bcs_tuple
            bcs_tuple[3] = goal_T
            self.inputs.write_bcs(*bcs_tuple)
            failed = 0
            self._raster_print(f'..T_rmin {self.windsoln.T_rmin*T_scale:.0f}, goal {goal_T*T_scale:.0f}')
            while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                failed += 1
                delta_T = (bcs_tuple[3]-self.windsoln.T_rmin)/2
                bcs_tuple[3] = self.windsoln.T_rmin+delta_T
                self.inputs.write_bcs(*bcs_tuple)
                self._raster_print(f'..T_rmin {self.windsoln.T_rmin*T_scale:.0f}, '
                      f'try {bcs_tuple[3]*T_scale:.0f}')
                if failed > 15:
                    self._normal_print("\nStruggling to substep towards new T_rmin")
                    return 1
        self._raster_print(f"   Successfully converged T_rmin to {self.windsoln.T_rmin*T_scale:.0f}")
        if _called_in_polish==False:
            self.converge_mol_atomic_transition(_called_in_ramp_bcs=True)        
        if integrate_out==True:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
            return self.run_wind()
        else:
            return 0
    
    
    
    def ramp_Rmin(self,goal_Rmin,integrate_out=False,_called_in_polish=False):
        """
        Ramps normalized Rmin, where Rmin is the base of the simulation and should ideally be a radius "below the wind", i.e., in the molecular region below ~1 microbar.

        Args:
            goal_Rmin (float): The desired Rmin value in normalized units of Rp.
            integrate_out (bool, optional): If True, integrates out to Coriolis radius after ramping. Defaults to False.

        Returns:
            int: 0 if successful, 1 if failed to converge.
        """
        if goal_Rmin > 30:
            self._normal_print("WARNING: Rmin should be in units of Rp. Returning...")
            return
        current_Rmin = self.windsoln.Rmin
#         goal_Rmin    = self.base_bcs()[0]
        delta = 0.02 
        if goal_Rmin < current_Rmin:
            delta *= -1
        while (abs(1.-self.windsoln.Rmin/goal_Rmin) > 5e-3):
            bcs_tuple = self.windsoln.bcs_tuple
            if abs(goal_Rmin - current_Rmin) < abs(delta):
                delta = goal_Rmin - current_Rmin
            step_Rmin = current_Rmin + delta
            bcs_tuple[0] = step_Rmin
            self.inputs.write_planet_params(*self.windsoln.planet_tuple)
            self.inputs.write_bcs(*bcs_tuple)
            failed = 0
            self._raster_print(f'..Goal Rmin: {goal_Rmin:.4f}Rp. Current: {self.windsoln.Rmin:.4f}Rp. Trying {step_Rmin:.4f}')
            fail = 0
            result = self.run_wind(expedite=True,calc_postfacto=False)
            while result != 0:
                fail+=1
                if fail==1:
                    self.converge_mol_atomic_transition(_called_in_ramp_bcs=True) #TEST
                if result == 4:
                    self.inputs.write_flags(*self.windsoln.flags_tuple,
                                            integrate_out=False)
                    self._normal_print("Turning off outward integration temporarily.")
                else:
                    bcs_tuple = self.windsoln.bcs_tuple
                    step_Rmin = current_Rmin + delta/(fail+1)
                    bcs_tuple[0] = step_Rmin
                    self.inputs.write_planet_params(*self.windsoln.planet_tuple)
                    self._raster_print(f'..Fail {fail:d}: Rmin {current_Rmin:.5g}Rp, trying {step_Rmin:.5g}Rp')
                    self.inputs.write_bcs(*bcs_tuple)
                if fail > 10:
                    self._normal_print(f"   Failed at {current_Rmin} Rp. Goal {goal_Rmin} Rp.")
                    return 1
                result = self.run_wind(expedite=True,calc_postfacto=False)
            current_Rmin = self.windsoln.Rmin
        self._raster_print(f"   Successfully converged Rmin to {self.windsoln.Rmin:.6f} Rp.\n")
        # self.converge_mol_atomic_transition(_called_in_ramp_bcs=True) #TEST
        if _called_in_polish==False:
            self.converge_mol_atomic_transition(_called_in_ramp_bcs=True)
        if integrate_out==True:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
            return self.run_wind()
        else:
            return 0
    
    
    
    def ramp_rho_rmin(self,goal_rho,integrate_out=False,_called_in_polish=False):
        """
        Ramps to normalized mass density at Rmin.

        Args:
            goal_rho (float): The desired rho value in normalized units of RHO0 (found in src/defs.h, default 1e-15).
            integrate_out (bool, optional): If True, integrates out to Coriolis radius after ramping. Defaults to False.

        Returns:
            int: 0 if successful, 1 if failed to converge.
        """
        rho_scale = self.windsoln.scales_dict['rho']
        if goal_rho < 10:
            self._normal_print("WARNING: Rho should be in units of RHO0 =  %.0e ~nBAR*MH/K/T0"%rho_scale)
            return
        failed = 0
        OOM_old = np.copy(np.floor(np.log10(self.windsoln.rho_rmin)))
        while (abs(1.-self.windsoln.rho_rmin/goal_rho) > 1e-10):
            # OOM_current = np.floor(np.log10(self.windsoln.rho_rmin))
            OOM_current = np.floor(np.log10(goal_rho)) #new: 11/11/25
            if OOM_old != OOM_current:
                self._normal_print(f'..Order of mag of rho_rmin has changed, C source code does not update dynamically. \nIf run takes a long time or fails to converge, halt run, and restart run to update rho_rmin convergence scale.')            
                # new_RHOSCALE = 10**np.floor(np.log10(self.windsoln.rho_rmin*0.001))
                new_RHOSCALE = 10**np.floor(np.log10(goal_rho*0.001))
                h = (open(self.path+'src/defs.h','r')).readlines()
                f = open(self.path+'src/defs.h','w')
                for idx,hline in enumerate(h):
                    splitline = hline.split()
                    if len(splitline) >= 2:
                        line_var = splitline[0]+' '+splitline[1]
                    if line_var == '#define RHOSCALE':
                        index3 = idx
                h[index3] = '#define RHOSCALE %.1f\n' %new_RHOSCALE
                f.writelines(h)
                f.close() 
                sub = Popen('make',cwd=self.path, stdout=PIPE, stderr=PIPE) 
                output, error_output = sub.communicate() #FIX (put output check)
                print(output,error_output)
#                 self._normal_print(error_output)
                OOM_old = np.copy(OOM_current)
                # OOM_old = np.copy(np.floor(np.log10(self.windsoln.rho_rmin)))

            if failed > 2:
                bcs_tuple = self.windsoln.bcs_tuple
                bcs_tuple[2] = self.windsoln.rho_rmin+(bcs_tuple[2]-self.windsoln.rho_rmin)/(10**failed) #here
                self.inputs.write_bcs(*bcs_tuple)
                self._raster_print(f'  ..Proceeding with smaller stepsize: rho_rmin {self.windsoln.rho_rmin:.5g}, goal {goal_rho:.5g}')
                fail2 = 0
                while self.run_wind(expedite=True,calc_postfacto=False) == 1: 
                    fail2 += 1
                    delta_rho = (bcs_tuple[2]-self.windsoln.rho_rmin)/10**(failed+fail2)
                    bcs_tuple[2] = self.windsoln.rho_rmin+delta_rho
                    self.inputs.write_bcs(*bcs_tuple)
                    self._raster_print(f'  ..Attempt {failed:d}: rho_rmin {self.windsoln.rho_rmin*rho_scale:.3e} g/cm3, '
                          f'try {bcs_tuple[2]*rho_scale:.3e}')
                    if fail2 > 8:
                        self._normal_print("\n  Struggling to substep towards new rho_rmin")
                        return 1
            else:
                bcs_tuple = self.windsoln.bcs_tuple
                bcs_tuple[2] = goal_rho
                self.inputs.write_bcs(*bcs_tuple)
                failed = 0
                self._raster_print(f' ..Goal rho_rmin: {goal_rho*rho_scale:.3e} g/cm3. Current: {self.windsoln.rho_rmin*rho_scale:.3e}.')

                while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                    failed += 1
                    if failed == 1:
                        self.converge_mol_atomic_transition(_called_in_ramp_bcs=True) #TEST
                    delta_rho = (bcs_tuple[2]-self.windsoln.rho_rmin)/10
                    bcs_tuple[2] = self.windsoln.rho_rmin+delta_rho
                    self.inputs.write_bcs(*bcs_tuple)
                    self._raster_print(f'  ..Attempt {failed:d}: rho_rmin {self.windsoln.rho_rmin*rho_scale:.3e} g/cm3, '
                          f'try {bcs_tuple[2]*rho_scale:.3e}')
                    if failed > 8:
                        self._normal_print("\n  Struggling to substep towards new rho_rmin")
                        return 1
        self._raster_print(f"   Successfully converged rho_rmin to {self.windsoln.rho_rmin*rho_scale:.5e} g/cm3")
        if _called_in_polish==False:
            self.converge_mol_atomic_transition(_called_in_ramp_bcs=True)
        if integrate_out==True:
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
            return self.run_wind()
        else:
            return 0
    
    
    
#Spectrum tools    
    def flux_norm(self,goal_flux_in_range,eV_range=[13.6,100],ramp=False,plot=True,
                 integrate_out=True,converge_bcs=True):
        """
        Computes the total flux (across the loaded spectrum range) necessary to achieve the goal_flux_in_range in the eV_range identified. If ramp=True, also ramps the solution to that Ftot.

        Args:
            goal_flux_in_range (float): Desired flux in ergs/s/cm2 in range given by eV_range.
            eV_range (list or array, optional): Range to normalize to in electron volts (13.6-100 is EUV). Defaults to [13.6, 100].
            ramp (bool, optional): If True, ramps the solution to the computed Ftot. Defaults to False.
            plot (bool, optional): If True, plots ramping progress. Defaults to True.
            integrate_out (bool, optional): If True, integrates out to Coriolis radius after ramping. Defaults to True.
            converge_bcs (bool, optional): If True, polishes boundary conditions after ramping. Defaults to True.

        Returns:
            int: 0 if ramp successful or already done, error codes (1, 4, etc.) if ramp unsuccessful.
            float: goal_total_flux if ramp=False (does not ramp).
        """
        euv = np.array(eV_range)*const.eV
        E = self.windsoln.E_wl
        flux_per_bin = self.windsoln.Ftot*self.windsoln.wPhi_wl*E
        current_euv_flux = sum(flux_per_bin[(E>euv[0])&(E<euv[1])])
        ratio = goal_flux_in_range/current_euv_flux
        goal_total_flux = ratio*self.windsoln.Ftot
        if ramp == True:
            if abs(goal_total_flux-self.windsoln.Ftot)/self.windsoln.Ftot > 0.01:
                return self.ramp_var("Ftot",goal_total_flux,make_plot=plot,
                                    converge_bcs=converge_bcs,integrate_out=integrate_out)
            else:
                self._normal_print("  Flux ramping already done.")
                return 0
        else:
            return goal_total_flux

        
        
    def ramp_spectrum(self,Fnorm=0.0,norm_spec_range=[],
                      goal_spec_range=[],units='eV',normalize=True,
                      kind='full',plot=False):
        """
        Ramps stellar spectrum wavelength/energy range to new wavelength/energy range.

        Args:
            Fnorm (float, optional): Flux in ergs/s/cm2 AT SEMIMAJOR AXIS OF PLANET. If 0.0, flux is normalized to the current value in norm_spec_range. Otherwise, ramps to given Fnorm. Defaults to 0.0.
            norm_spec_range (list or array, optional): Desired range over which to normalize the flux in units of 'units'. Defaults to [].
            goal_spec_range (list or array, optional): Custom upper and lower limits of spectrum in units of 'units'. Defaults to [].
            units (str, optional): Units of range values. Options are 'cm', 'nm', 'eV'. Defaults to 'eV'.
            normalize (bool, optional): If True, flux in new range will be normalized to Fnorm in norm_spec_range. Defaults to True.
            kind (str, optional): 'full' or 'mono' - spectrum frequency type. Defaults to 'full'.
            plot (bool, optional): If True, plots ramping progress. Defaults to False.

        Returns:
            int: 0 if ramping successful, 1 if ramping with smaller stepsize unsuccessful, 2 if need to provide Ftot or set norm_flux=True.

        Example:
            >>> ramp_spectrum(1000,norm_spec_range=[13.6,100],goal_spec_range=[13.6,2000],units='eV',normalize=True)
            Will generate a solution with XUV spectrum over 13.6-2000 eV and total integrated flux over 13.6-2000 eV such that the flux in the norm_spec_range (13.6-100 eV) is 1000 ergs/s/cm^2. 
        """
        self.load_spectrum()
        ons = np.copy(np.array([self.spectrum.norm_span[0], self.spectrum.norm_span[1]]))
        wl_norm = 1e-7
        if (kind == 'mono') and (self.windsoln.spec_kind == 'multi'):
            self._normal_print("WARNING: Trouble ramping from multifrequency to monofrequency solutions. Start from a monofrequency solution and ramp flux using ramp_var('Ftot',...).")
            return
        
        if self.windsoln.nspecies > 4:
            self._normal_print("WARNING: Ramping spectrum with large number of metals can be expensive and fail.\n  Suggestion: Ramp spectrum for H,He version of planet, then add metals.")
               
        if len(norm_spec_range)== 0:
            norm_spec_range = goal_spec_range
            
        def format_range(spec_range,units):
            #converts bounds from given units to 'nm' and returns energy bounds for easy normalization
#             if len(spec_range) == 0:
#                 spec_range = self.windsoln.spec_resolved/self.spectrum.wl_norm
#                 units = 'nm'
#                 E_top = const.hc/(spec_range[1]*wl_norm) #ergs
#                 E_bottom = const.hc/(spec_range[0]*wl_norm)
#                 goal_span = spec_range
#             else:
            spec_range = np.array(spec_range)

            if spec_range[0] > spec_range[1]:
                self._normal_print('WARNING: Spectrum range limits should be in ascending order')
                return 3
            if units == 'eV':
                E_top = spec_range[1]*const.eV #ergs
                E_bottom = spec_range[0]*const.eV
                goal_span = np.flip(const.hc/(spec_range*const.eV)/wl_norm)
            elif units == 'nm':
                E_top = const.hc/(spec_range[0]*wl_norm) #ergs
                E_bottom = const.hc/(spec_range[1]*wl_norm)
                goal_span = spec_range
            elif units == 'cm':
                E_top = const.hc/(spec_range[0]) #ergs
                E_bottom = const.hc/(spec_range[1])
                goal_span = spec_range / wl_norm
            else:
                sys.exit("Invalid units. Currently only equipped for 'cm','nm', 'eV'.")

            return goal_span, E_bottom, E_top
        
        goal_span, E_bot_goal, E_top_goal = format_range(goal_spec_range,units)
        norm_span, E_bot_norm, E_top_norm = format_range(norm_spec_range,units)
        curr_span = self.windsoln.spec_resolved/self.spectrum.wl_norm
        
        if normalize == True:
            if Fnorm == 0.0:
                if (norm_span[0]<curr_span[0]) or (norm_span[1]>curr_span[1]):
                    self._normal_print('When Fnorm=0.0, flux is normalized to the current value in norm_spec_range.')  
                    self._normal_print(f'However, bounds of normalization range exceed the current bounds of spectrum. (Norm: [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm. Current: [{curr_span[0]:.2f},{curr_span[1]:.2f}]nm.)')
                    self._normal_print('So, flux will be normalized to total flux in current spectrum range.')
                    norm_span[0] = max(norm_span[0],curr_span[0])
                    norm_span[1] = min(norm_span[1],curr_span[1])
                #set Fnorm to the current value in the normalization range
                E = self.windsoln.E_wl
                flux_per_bin = E*self.windsoln.Ftot*self.windsoln.wPhi_wl
                Fnorm = sum(flux_per_bin[(E<E_top_norm) & (E>E_bot_norm)])
                flux_str = f'Spectrum will be normalized such that sum(Flux[{norm_span[0]:.2f}, {norm_span[1]:.2f}]nm) = {Fnorm:.0f} ergs/s/cm2.'
                self._normal_print(flux_str)
        

        if kind != self.windsoln.spec_kind:
            self._normal_print("Warning: Ramper sometimes has difficulty changing from 'full' to 'mono'.")
            self._normal_print("Consider ramping from existing monofrequency solution.")
        
        if self.windsoln.spec_src_file == 'scaled-solar':
            spec = spectrum(date=self.windsoln.spec_date)
        else:
            spec = spectrum(lisird=False,spectrum_file=self.windsoln.spec_src_file)
        for sps in self.windsoln.species_list:
            spec.add_species(sps)

        self._normal_print(f'Goal: {goal_span} nm')
        spec.set_resolved(*goal_span) #
        spec.set_normalized(*goal_span) 
        spec.set_window(*goal_span,kind=kind)
        spec.generate(kind=kind, savefile=self.path+'inputs/spectrum.inp')
        self.generate_rate_coeffs()

        if self.run_wind(expedite=True,calc_postfacto=False) == 0:
            if normalize == True:
                self._raster_print(f'Ramped spectrum wavelength range, now normalizing spectrum. \n ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm') 
                ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                          const.hc/(norm_span[0]*wl_norm)/const.eV]
                final_result = self.flux_norm(Fnorm,ranges,
                                   ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
            else:
                final_result = 0
        else:
            #If it cannot make it in one leap, ramp spectrum range
            delta = np.copy(goal_span - self.windsoln.spec_resolved/self.spectrum.wl_norm)
            gen_fail = 0 
            last_spec = self.windsoln.spec_resolved
            while abs(sum((goal_span-self.windsoln.spec_resolved/self.spectrum.wl_norm)/goal_span)) > 1e-3:
                avg_diff = abs(sum((goal_span-self.windsoln.spec_resolved/self.spectrum.wl_norm)/goal_span))
                curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                step_span = self.windsoln.spec_resolved/self.spectrum.wl_norm + delta / 5
                if abs(sum((goal_span-step_span)/goal_span)) > avg_diff:
                    step_span = goal_span
                self._raster_print(f'  Fail 1: Current [{curr[0]:.1f},{curr[1]:.1f}]. Attempting [{step_span[0]:.1f},{step_span[1]:.1f}]')
                spec.set_resolved(*step_span) #
                spec.set_normalized(*step_span) 
                spec.set_window(*step_span,kind='full')
                spec.generate(kind='full', savefile=self.path+'inputs/spectrum.inp')
                self.generate_rate_coeffs()
                
                if np.array_equal(self.windsoln.spec_resolved,last_spec):
                    if gen_fail > 5:
                        self._normal_print("Spectrum ramping cannot proceed. Try ramping spectrum for an H,He atmosphere then add metals.")
                        return 1
                    gen_fail += 1
                    delta /= 3 
                last_spec = np.copy(self.windsoln.spec_resolved)

                fail = 1
                while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                    fail += 1
                    factor = 1/(5*2**(fail-1))
                    if fail > 6:
                        curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                        eV = const.hc/self.windsoln.spec_resolved / const.eV
                        self._normal_print(f'Failure: Too many attempts. Final range: [{curr[0]:.2f},{curr[1]:.2f}]nm ([{eV[1]:.2f},{eV[0]:.2f}]eV).')
                        self._normal_print('          To ramp manually, see tutorial.')
                        if normalize == True:
                            self._raster_print(f'Now normalizing spectrum. \n ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm') 
                            ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                                    const.hc/(norm_span[0]*wl_norm)/const.eV]
                            result =  self.flux_norm(Fnorm,ranges,
                                                ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
                            curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                            self._normal_print(f'Final Ftot at planet across [{curr[0]:.2f},{curr[1]:.2f}]nm = {self.windsoln.Ftot:.0f} ergs/s/cm2.')
                        return 1
                    curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                    step_span = curr + delta*factor
                    if abs(sum((goal_span-step_span)/goal_span)) > avg_diff:
                        step_span = goal_span
                    self._raster_print(f'  Fail {fail}: Current [{curr[0]:.1f},{curr[1]:.1f}]. Attempting {step_span}nm')
                    spec.set_resolved(*step_span) #
                    spec.set_normalized(*step_span) 
                    spec.set_window(*step_span,kind='full')
                    spec.generate(kind='full', savefile=self.path+'inputs/spectrum.inp') 
                    self.generate_rate_coeffs()

            else:
                if normalize == True:
                    self._raster_print('  Ramped spectrum wavelength range, now normalizing spectrum.'+self.clear+f'  ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm') 
                    ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                            const.hc/(norm_span[0]*wl_norm)/const.eV]
                    result =  self.flux_norm(Fnorm,ranges,
                                        ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
                    curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                    self._normal_print(f'Final Ftot at planet across [{curr[0]:.2e},{curr[1]:.2e}]nm = {self.windsoln.Ftot:.0f} ergs/s/cm2.')
                    final_result = result
                    # return result
                else:
                    final_result = 0
                    # return 0
        
        if plot == True:
            fig, ax = plt.subplots()
            max_mk = ((self.spectrum.data_norm['wl']>=
                   min(ons[0], self.spectrum.norm_span[0]))
                  &(self.spectrum.data_norm['wl']<=
                    max(ons[1], self.spectrum.norm_span[1])))
            mk = ((self.spectrum.data_norm['wl']>=self.spectrum.norm_span[0])
                  &(self.spectrum.data_norm['wl']<=self.spectrum.norm_span[1]))

            v3 = ax.axvline(self.spectrum.norm_span[0], zorder=-2, ls='-', c='r',
                            lw=3, label='Normalized')
            ax.axvline(self.spectrum.norm_span[1], zorder=-2, ls='-', c='r', lw=3)
            # Plot old spans (reduced alpha)
            ax.axvline(ons[0], zorder=0, ls='--', c='y', lw=3, alpha=0.5)
            ax.axvline(ons[1], zorder=0, ls='--', c='y', lw=3, alpha=0.5)
            # Plot spectrum
            s0, = ax.plot(self.spectrum.data_norm['wl'][max_mk],
                         self.spectrum.data_norm['phi_wl'][max_mk], lw=1,
                         label='Spectrum')
            s1, = ax.plot(self.spectrum.data_norm['wl'][mk],
                         self.spectrum.data_norm['phi_wl_smth'][mk], lw=2,
                         label='Smoothed')
            # Plot arrows if changed
            ax2 = ax.twinx()
            ax2.get_yaxis().set_visible(False)
            compare = [[ons, self.spectrum.norm_span]]
            sc = ['y', 'm', 'r']
            oldplot = False
            lines = [s0, s1, v3]
            for i, span in enumerate(compare):
                old = span[0]
                new = span[1]
                for o_s, n_s in zip(old, new):
                    if abs(o_s-n_s)/o_s > 1e-2:
                        ax2.arrow(o_s, 0.75+0.1*i, n_s-o_s, 0, alpha=0.3,
                                  color=sc[i], width=0.015, head_width=0.05,
                                  head_length=3, length_includes_head=True)
            ax.set_yscale('log')
            ax.set_title('Spectrum ramp')
            ax.set_xlabel(r'Wavelength ($\lambda$) [nm]')
            ax.set_ylabel('Normalized spectral photon irradiance\n'
                          r'$\left(\phi_{\lambda}\right)$ [s$^{-1}$ cm$^{-2}$]')
            fig.tight_layout(pad=0.3)
            fig.subplots_adjust(bottom=0.3, top=0.9)
            labels = [l.get_label() for l in lines]
            fig.legend(lines, labels, bbox_to_anchor=(0.5, 0.0), loc='lower center',
                       ncol=3)
            plt.show()

        return final_result


    def ramp_to_user_spectrum(self,spectrum_filename,species_list=[],
                              updated_F=0.0,norm_spec_range=[],goal_spec_range=[],
                              units='eV',normalize=True,plot=True,ramp_range=True,
                              hires_savgol_window=None):
        '''Ramping stellar spectrum wavelength/energy range to new wavelength/energy range.
        
        Args:
            spectrum_filename (str): Name of the formatted spectrum file saved in McAstro/stars/additional_spectra/.
            species_list (list of str, optional): Species list for spectrum binning. Defaults to [].
            updated_F (float, optional): Summed flux over norm_spec_range in ergs/s/cm² at planet's semimajor axis. Defaults to 0.0.
            norm_spec_range (list or array, optional): Desired range over which to normalize, in units of 'units'. Defaults to [].
            goal_spec_range (list or array, optional): Custom upper and lower limits of spectrum, in units of 'units'. Defaults to [].
            units (str, optional): Units for range values ('eV', 'cm', 'nm'). Defaults to 'eV'.
            normalize (bool, optional): If True, flux in new range will be normalized to Fnorm in norm_spec_range. Defaults to True.
            kind (str, optional): Spectrum frequency type ('full' or 'mono').
            plot (bool, optional): If True, plot the spectrum. Defaults to False.
            ramp_range (bool, optional): If False, will not ramp spec range or flux. Defaults to True.
            hires_savgol_window (int, optional): Window length for Savitzky–Golay smoothing if user spectrum has >10,000 bins. Must be odd. Defaults to 501 (Higher = more smoothing).
                Once you have identified an appropriate value, it will be automatically saved to the McAstro/stars/additional_spectra/spectrum_filename for future loads of that spectrum.

        Returns:
            None
        '''
        #generating a smoothed and binned version of the user-input code
        wl_norm=1e-7
        if hires_savgol_window is not None:
            if hires_savgol_window%2==0: #must be odd
                hires_savgol_window+=1
                
            f = open(self.path+'McAstro/stars/spectrum/additional_spectra/'+spectrum_filename,'r')
            h = f.readlines()
            h[2] = '%d\n' %hires_savgol_window
            f.close()
            g = open(self.path+'McAstro/stars/spectrum/additional_spectra/'+spectrum_filename,'w')
            g.writelines(h)
            g.close()
        
        spec = spectrum(lisird=False,spectrum_file=spectrum_filename,
                        wl_norm=wl_norm,print_warning=True)
        
        if len(species_list) == 0:
            species_list = self.windsoln.species_list
        spaced_list = McAtom.formatting_species_list(species_list)
        unspaced_list = [sp.replace(' ','') for sp in spaced_list]
        for species in unspaced_list: 
            spec.add_species(species)
        #ramp first to same range as existing solution for interpolation 
        spec.set_resolved(*self.windsoln.spec_resolved/wl_norm) 
        spec.set_normalized(*self.windsoln.spec_normalized/wl_norm) 
        spec.set_window(*self.windsoln.spec_window/wl_norm,kind='full')
        spec.generate(kind='full', savefile=self.path+'inputs/goal_spectrum.inp')
        spec.generate(kind='full', savefile=self.path+'inputs/spectrum.inp')
        
        old_E = self.windsoln.sim_spectrum['E']
        old_wPhi = self.windsoln.sim_spectrum['wPhi']
        
        self.inputs.write_flags(*self.windsoln.flags_tuple,
                                integrate_out=False)
        if self.run_wind() != 0:
            new = np.genfromtxt(self.path+'inputs/goal_spectrum.inp',
                                delimiter=',',skip_header=10)


            start_wPhi = si.Akima1DInterpolator(np.flip(old_E),
                                                np.flip(self.windsoln.sim_spectrum['wPhi']))
            goal_wPhi = si.Akima1DInterpolator(np.flip(new[:,0]),np.flip(new[:,1]))
            wPhi = goal_wPhi(old_E)
            wPhi[0] = wPhi[1] #deals with edge errors
            wPhi[-1] = wPhi[-2]

            percent=1
            delta_wPhi = np.copy((goal_wPhi(old_E) - self.windsoln.sim_spectrum['wPhi']))
            step_wPhi = np.array(self.windsoln.sim_spectrum['wPhi'] + delta_wPhi*percent)
            step_wPhi[0] = step_wPhi[1] #deals with errors at edges
            step_wPhi[-1] = step_wPhi[-2]

            self.inputs.write_spectrum(*self.windsoln.spectrum_tuple[:3],
                                      spectrum_filename,
                                      *self.windsoln.spectrum_tuple[4:10],
                                      step_wPhi,
                                      *self.windsoln.spectrum_tuple[11:])
            self.inputs.write_flags(*self.windsoln.flags_tuple,
                                    integrate_out=False)
            result = self.run_wind()

            avg = abs(np.average((wPhi - self.windsoln.sim_spectrum['wPhi'])))
            try:
                while abs(np.average((wPhi - self.windsoln.sim_spectrum['wPhi'])/self.windsoln.sim_spectrum['wPhi']))>1e-5:
                    spec.generate(kind='full', savefile=self.path+'inputs/spectrum.inp')
                    self.generate_rate_coeffs()
                    if self.run_wind() != 0:
                        fail = 0
                        percent = 0.1

                        current_wPhi = np.copy(self.windsoln.sim_spectrum['wPhi'])
                        step_wPhi = np.array(current_wPhi+delta_wPhi*(percent))
                        #if next step is larger than the remainder, set step to goal
                        if avg < abs(np.average(delta_wPhi*percent)):
                            step_wPhi = wPhi
                        if plot == False:
                            self._raster_print(f"  ..Trying a {percent*100:.0f}% step")
                        step_wPhi[0] = step_wPhi[1]
                        step_wPhi[-1] = step_wPhi[-2]
                        self.inputs.write_spectrum(*self.windsoln.spectrum_tuple[:3],
                                                  spectrum_filename,
                                                  *self.windsoln.spectrum_tuple[4:10],
                                                  step_wPhi,
                                                  *self.windsoln.spectrum_tuple[11:]) 
                        self.inputs.write_flags(*self.windsoln.flags_tuple,
                                                integrate_out=False)
                        while self.run_wind(expedite=True,calc_postfacto=False) != 0:
                            fail+=1
                            self._raster_print(f"  Fail {fail:d}: delta = {percent/(fail+1):.3f}")
                            step_wPhi = np.array(current_wPhi+delta_wPhi*(percent/(fail+1)))
                            step_wPhi[0] = step_wPhi[1]
                            step_wPhi[-1] = step_wPhi[-2]
                            self.inputs.write_spectrum(*self.windsoln.spectrum_tuple[:3],
                                                      spectrum_filename,
                                                      *self.windsoln.spectrum_tuple[4:10],
                                                      step_wPhi,
                                                      *self.windsoln.spectrum_tuple[11:]) 
                            self.inputs.write_flags(*self.windsoln.flags_tuple,
                                                    integrate_out=False)
                            if fail > 10:
                                self._normal_print("Failed to ramp to new stellar spectrum.")
                                return 1

                        self.converge_Ncol_sp(expedite=True,quiet=True)
                        self.converge_mol_atomic_transition()
                        self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)

                        avg = abs(np.average((wPhi - self.windsoln.sim_spectrum['wPhi'])))
                        if plot == True:
                            plt.clf()
                            plt.plot(const.hc/old_E/1e-7,start_wPhi(old_E),label='Original')
                            plt.plot(const.hc/old_E/1e-7,wPhi,label='Goal')
                            plt.plot(const.hc/old_E/1e-7,self.windsoln.sim_spectrum['wPhi'],
                                         label='Current',ls='--')
                            plt.xlabel('Wavelength (nm)')
                            plt.ylabel('Photon Density')
                            plt.legend()
                            plt.yscale('log')
                            plt.title(f" Success! Average difference now {avg:.2e}")
                            display.display(plt.gcf())
                            display.clear_output(wait=True)
                        else:
                            self._raster_print(f" Success! Average difference now {avg:.2e}")
            except ValueError:
                self._normal_print(f"Succeeded in full jump to {spectrum_filename}")
        self._normal_print("*** REMINDER: Goal spectrum (orange) oversmoothed? Lower the hires_savgol_window value. ***")
        self._normal_print("              (Note: lower values will take longer to run).")


        #Ramping to desired spectral range and flux
        self.converge_Ncol_sp(expedite=True,quiet=True)
        self.converge_mol_atomic_transition()
        self._normal_print(f'\nSuccess ramping to {spectrum_filename} spectrum shape!')
        if plot == True:
            plt.plot(const.hc/old_E/1e-7,old_wPhi,label='Original')
            plt.plot(const.hc/self.windsoln.sim_spectrum['E']/1e-7,self.windsoln.sim_spectrum['wPhi'],
                         label='Current')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel(r'Photon Density at 1 au')
            plt.legend()
            plt.yscale('log')
            plt.show()

        if ramp_range==True:
            self._normal_print('     Now ramping flux and spectral range.')
            self.load_spectrum()
            if len(goal_spec_range) == 0:
                goal_spec_range = self.windsoln.spec_resolved/self.spectrum.wl_norm
                units = 'nm'
            if updated_F == 0:
                f = open(self.path+"McAstro/stars/spectrum/additional_spectra/"+spectrum_filename,"r")
                updated_F = float(f.readlines()[1])/(self.windsoln.semimajor/const.au)**2
                norm_spec_range = [12.4,91.1]
                f.close()
            #Ramping the spectrum range and flux in normalized range
            result = self.ramp_spectrum(updated_F,norm_spec_range,
                                        goal_spec_range,units,normalize,
                                        kind='full',plot=plot)
            if (result == 0) or (result==4):
                self._normal_print("Success! Ramped to user-input stellar spectrum "+spectrum_filename)
                return 0
            else:
                self._normal_print("Ramped successfully to user-input spectrum, but ramping spectral range and/or flux was not successful.")
                return 2
        else:
            return 0

        
        
    def format_user_spectrum(self,wl,flux_1au,wl_units='cm',spectrum_name='',
                            comment='',overwrite=False):
        """
        Users may import a custom spectrum. This function formats the spectrum to be readable by the code.

        Args:
            wl (list or array): Wavelengths.
            flux_1au (list or array): Flux per wavelength bin at 1 AU from star in ergs/s/cm2. To scale from Earth, multiply by (distance to star in au / 1 au)^2. Not flux density (to convert from flux density, multiply by delta(wl)).
            wl_units (str, optional): Units for wavelength. Options are 'cm', 'nm', 'm', 'A' (angstroms). Defaults to 'cm'.
            spectrum_name (str, optional): Name of file. Will be saved in solution spectrum. Defaults to ''.
            comment (str, optional): Any comments (name of spectrum and other info). No need to add newline character. Defaults to ''.
            overwrite (bool, optional): If True, overwrites existing file of same name. Defaults to False.

        Returns:
            None
        """
        spectrum_name = spectrum_name.replace(' ','_')
        wl = np.array(wl)
        flux_1au = np.array(flux_1au)

        if len(wl) != len(flux_1au):
            self._normal_print("Wavelength and flux at 1 au arrays must be same length.")
            return    

        file = self.path+'McAstro/stars/spectrum/additional_spectra/'+spectrum_name

        if (exists(file)==True) and (overwrite == False):
            self._normal_print("To overwrite existing file, set overwrite=True.")
            return

        g = open(file,'w')
#         if len(comment) > 0:
        g.write('#'+comment+'\n')

        if wl_units == 'm':
            wl_cm = wl*1e2
        elif wl_units == 'nm':
            wl_cm = wl*1e-7
        elif wl_units == 'A':
            wl_cm = wl*1e-8
        elif wl_units == 'cm':
            wl_cm = wl
        else:
            self._normal_print("Invalid units: Options = 'A','nm','cm',or 'm'.")
            return
        if np.median(flux_1au) < 1e-8:
            self._normal_print("WARNING: Flux at 1 au from star should be in units of ergs/s/cm2 (check that your input is not the flux at Earth).")
            
        if wl_cm[-1] > 1e-4:
            self._normal_print(f"WARNING: Double check - Input wavelength array might not be in units of {wl_units:s}.")
        
        total_EUV_flux = sum(flux_1au[(wl>min(np.max(wl),1.2e-6))&(wl<max(np.min(wl),9.1e-6))])
        g.write('%.5f\n'%total_EUV_flux)
        g.write('0\n')
        
        #Converting to ergs/s/cm2/cm flux density
        delta_lam = np.diff(wl_cm,prepend=wl_cm[0]-(wl_cm[1]-wl_cm[0]))
        flux_dens_1au = flux_1au/delta_lam
        
        df = pd.DataFrame(np.column_stack((
            wl_cm,
            flux_dens_1au,
            np.zeros_like(wl_cm),
            const.c/wl_cm,
            flux_dens_1au*wl_cm**2/const.c,
        )), columns=['wl','F_wl','unc','nu','F_nu'])

        df.to_csv(g, header=True, index=False)
        g.close()

        self._normal_print('Spectrum saved in Wind-AE readable format at '+file)

        return
    
    #Assorted helpful functions 
    def energy_plot(self,ax=0,alpha=0.8,all_terms=False,
                    CII_line_cool=False,CIII_line_cool=False,OII_line_cool=False,OIII_line_cool=False,
                    legend=True,sub_sonic=True):
        """ Plots energy balance terms used in the energy equation (Broome et al. 2025)

        Args: 
            windsoln: The wind solution object containing the simulation data. ("sim.windsoln")
            ax: The axis to plot on (default is 0, which creates a new figure)
            alpha: Transparency level for the plot lines (default is 0.8). 
                    Useful when overplotting multiple on same axes
            all_terms: If True, plot terms not included in Wind-AE, e.g., free-free cooling, conduction (default is False)
            CII_line_cool: If True, include CII line cooling terms (default is False)
            CIII_line_cool: If True, include CIII line cooling terms (default is False)
            OII_line_cool: If True, include OII line cooling terms (default is False)
            OIII_line_cool: If True, include OIII line cooling terms (default is False)
            legend: If True, display the legend (default is True)
            sub_sonic: If True, sets x-axis upper limit at sonic point radius

        Returns:
            None
        """
        energy_plot(self.windsoln,ax,alpha,all_terms,
                    CII_line_cool,CIII_line_cool,OII_line_cool,OIII_line_cool,
                    legend,sub_sonic)
        return
    
    def six_panel_plot(self,Mdot_legend=True,c='k',ls='-',label='',label_dim=[0,1.3,2],
                        ion_label=True,first_plotted=True,ax=0):
        '''
        Plots density (g/cm3), temperature (K), velocity (10 km/s), ionization fraction, column density (g/cm2), and number density (1/cm2), as a function of r (Rp).
            
        Args:
            soln - windsoln object (sim.windsoln)
            Mdot_legend - Bool; if True, put Mdot in legend of plot. Else, just prints.
            c - str; line color 
            ls - str; line style
            label - str; line label 
            label_dim - list; default=[0,1.3,2]. Location of label and ncols [x,y,ncols]. 
            first_plotted - Bool; True if this the first of many OR the ONLY SixPlot 
                            to be plotted on the same axes. 
            ax - matplotlib axis obj; if first_plotted=False, provide axis object so this 
                will be be plotted on desired figure with other simulations for comparison
        Returns:
            ax - axes object (if first_plotted=True)
            
        Example:
            ax1 = SixPlot(sim1.windsoln, first_plotted=True)
            SixPlot(sim2.windsoln, ax=ax1)
            SixPlot(sim3.windsoln, ax=ax1)
        '''
        six_panel_plot(self.windsoln,Mdot_legend,c,ls,label,label_dim,
                        ion_label,first_plotted,ax)
        return
    
    def quick_plot(self, Mdot_legend=True, c='k', ls='-', label='',label_dim=[0,1.3,2],
                ion_label=True,first_plotted=True, ax=0): 
        """
        Produces a velocity, density, temperature, and neutral fraction plot from intermediate solutions while ramping. For a more aesthetic four-panel plot use `quickplot()`.

        Args:
            windsoln: windsoln object to be plotted.
            ax (np.ndarray): The ax array on which plots are placed (must have shape = (2,2)).
            ax_Ys (matplotlib.axes.Axes, optional): Axis for plotting ionization fraction. Defaults to None.
            label (str, optional): Label to be added to ax[0][0]. Defaults to None.
            alpha (float, optional): Alpha of lines plotted. Defaults to 0.8.
            pvar (list of lists of str, optional): Variables to plot. Defaults to [['v','rho'],['T','Ys_HI']].
            norm (list of lists of floats, optional): Values to norm (divide) pvars with. Defaults to [[1e5,1e0],[1e0,1e0]].
            radius_prefix (str, optional): SI prefix for x-axis (None defaults to units of Rp). Defaults to None.
            sub_sonic (bool, optional): Only plot sub-sonic region. Defaults to False.
            past_rmin (bool, optional): Only plot past Rmin region. Defaults to False.
            sonic_vert (bool, optional): Add vertical line at sonic point. Defaults to True.

        Returns:
            None
        """
        quick_plot(self.windsoln,Mdot_legend,c,ls,label,label_dim,
                        ion_label,first_plotted,ax)
        return

    def integrate_out(self, quiet=False):
        '''Forces solution to integrate out past sonic point to the Rmax. May not preserve self-consistency of per-species column density at sonic point (Ncol_sp) as it may be necessary to raise Ncol_sp in order to integrate outwards.
        By default Rmax is set to currently-computed R_cori, but can be manually set Rmax to desired value via:
        Example:
            sim.inputs.write_planet_params(sim.windsoln.Rmin,Rmax,*sim.windsoln.bcs_tuple[2:])
            sim.run_wind()

        Args:
            quiet (bool, optional): If True, suppress output messages.
        Returns:
            int: 0 if successful, 1 if relaxation error, 4 if Ncol_sp had to be raised
        '''
        self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=True)
        result = self.run_wind()
        if result==4:
            factor = 2
            self.raise_Ncol_sp(by_factor=2)
            self._raster_print(f"  Raising Ncol_sp x{factor:d} to allow outward integration to the currently computed R_cori.")
            while self.run_wind()==4:
                factor+=2
                self.raise_Ncol_sp(by_factor=factor)
                self._raster_print(f"  Raising Ncol_sp x{factor:d} to allow outward integration to the currently computed R_cori.")
            else:
                if not quiet:
                    self._normal_print(f"Success! Ncol_sp raised x{factor:d} from self consistent value to allow outward integration")
                return 0
        elif result==1:
            if not quiet:
                self._normal_print("Relaxation Error")
            return 1
        elif result == 0:
            if not quiet:
                self._normal_print("Success! Integrated outwards.")
            return 0
        
        
    def _make_string_from_scinum(self, x):
        exp = np.log10(x)
        base = x/(10**math.floor(exp))
        if base != int(base):
            return int(x), 0
        base = int(base)
        exp = int(np.log10(x/base))
        return base, exp     
    