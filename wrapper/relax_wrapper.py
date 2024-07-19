#!/usr/bin/env python3
"""
relax_wrapper.py:
    Contains the class for wrapping the relaxation code. Handles loading,
    saving, running, and converging boundary conditions for the relaxation
    simulations. The goal is find the solution of interest to the end-user
    and automatic the choices one must make to be intelligently picked and
    consistent, i.e., the boundary conditions and ramping thru parameter
    space. Additionally the solution should be readily accessible in python.
    Handling of the solution is done by the other major class wind_solution
    found in wrapper_utils/windsoln.py.
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
import pylab as pl
from IPython import display

from .wrapper_utils import constants as const
from .wrapper_utils.atmosphere import atmosphere
from .wrapper_utils.system import system
from .wrapper_utils.physics import physics
from .wrapper_utils.inputs import input_handler
from .wrapper_utils.windsoln import wind_solution
from .wrapper_utils.plots import four_panel_plot
from .wrapper_utils.plots import energy_plot
from .wrapper_utils.metals import metal_class
# from .wrapper_utils import metals as metal
import McAstro.atoms.atomic_species as McAtom
import pandas as pd
import time

McAstro_loaded = False
if os.path.exists('McAstro/planets/insolation/glq_rates.py'):
    McAstro_loaded = True
    from .wrapper_utils.spectrum import spectrum


class wind_simulation:
    def __init__(self, csv_file='inputs/guess.inp', name='Init. Planet',
                 expedite=True):
        self.strat_kappa = 1e-2
        self.inputs = input_handler()
        self.first_print = True
        self.skip = False

    def load_nonexpedited(self, csv_file=None):
        if csv_file is None:
            csv_file = self.last_load
        self.windsoln = wind_solution(file=csv_file,expedite=False)
        return


    def load_uservars(self, csv_file=None):
        '''
        Loads all user variables. Use this for plotting values yourself.
        It does not rewrite input parameter files, so can be run at the
        same time as simulations without interrupting the self.
        
        Suggested syntax:
        sim = wind_sim()
        self.load_uservars(FILENAME)
        soln = self.windsoln
        soln.soln['Ncol_HI'], etc.
        '''
        if csv_file is None:
            csv_file = self.last_load
        self.windsoln = wind_solution(file=csv_file, add_uservars=True)
        return


    
    def load_planet(self, csv_file, expedite=True, name='Loaded Planet',print_atmo_composition=True):
        '''
        Loading planet solution file into inputs/guess.inp as new guess. 
        Also loads planet parameters into input parameter folders. Ramping 
        to new planet with ramp_to() will rewrite param folders for the 
        desired new planet, but the guess will remain the same.
        '''
        # load header and data
#         self.failed_deeper_bcs_ramp = False
        planet = wind_solution(file=csv_file, expedite=expedite)
        for j in range(planet.nspecies): #had a problem with adding spaces to e.g., 'He I'
            planet.species_list[j] = (planet.species_list[j]).replace(' ','')
            
        #Changing Nspecies, M (# pts in relaxation region), RHOSCALE (convergence condition) and remaking C code    
        f = open('src/defs.h', 'r') 
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
        
        f = open('src/defs.h','w')
        h = (open('src/defs-master.h','r')).readlines()
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
            print(f'\rNspecies has changed from %d to %d. Remaking C code...' %(nspecies_def, nspecies_new),
                 end='                                                                                      ')
            remake = True
        #For regridded solutions, the number of points in relaxtion region may change
        if og_length != new_length:
            print("Number of points in relaxation region has changed from %d to %d. Remaking C code..." 
                  %(og_length, new_length))
            remake=True
        #For solutions with very high rho at the lower boundary, the rho convergence condition should be raised
        if new_RHOSCALE != og_RHOSCALE:
            print("RHOSCALE (convergence condition) has changed from %d to %d. Remaking C code..." 
                  %(og_RHOSCALE, new_RHOSCALE)) 
            remake=True
        #Future users can add additional parameters
        if new_N_ADD_PARAMS != og_N_ADD_PARAMS:
            print("N_ADD_PARAMS (num. of additional params added by user) has changed from %d to %d. Remaking C code..." 
                  %(og_N_ADD_PARAMS, new_N_ADD_PARAMS)) 
            remake=True
        if remake == True:
            sub = Popen('make', stdout=PIPE, stderr=PIPE) 
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
        # load planet as guess.inp if not already
        if csv_file != 'inputs/guess.inp':
            sub = Popen(["cp", csv_file, 'inputs/guess.inp'],
                        stdout=PIPE, stderr=PIPE)
            sub.wait()
            output, error_output = sub.communicate()
            if error_output:
                print(f'ERROR: Failed loading {csv_file:s} as guess.inp\n')
                print(error_output)
                return
        # If all is sucessful update wind_simulation object
        sub = Popen(["cp", 'inputs/guess.inp', 'saves/windsoln.csv'], stdout=PIPE, stderr=PIPE)
        sub.wait()
        output, error_output = sub.communicate()
        if error_output:
            print(f'ERROR: Failed copying {csv_file:s} into saves/windsoln.csv\n')
            print(error_output)
            return
        
        self.guess = planet
        self.windsoln = planet
        self.system = system(*planet.planet_tuple, name=name)
        self.mu = planet.calc_mu_base()
        self.atmosphere = atmosphere(self.system,
                                     self.guess.T_rmin*self.guess.scales[3],
                                     planet.calc_mu_base(),
                                     kappa_opt=self.strat_kappa)
        self.atmosphere.windbase_radius(self.guess.rho_rmin
                                        *self.guess.scales[1])
        self.physics = physics(*planet.physics_tuple)
        if print_atmo_composition==True:
            print("Atmosphere Composition")
            print('  Species:   '+(',        '.join('%s' %sp.replace(' ','') for sp in planet.species_list)))
            print('  Mass frac: '+(', '.join('%.2e' %hx for hx in planet.HX)))
        
        self.ramp_class = "system"
        if McAstro_loaded:
            self.load_spectrum()
        self.last_load = csv_file

        

    def load_spectrum(self, generate=False, wl_norm=1e-7):
        if self.windsoln.spec_src_file != 'scaled-solar':
            self.spectrum = spectrum(lisird=False,spectrum_file=self.windsoln.spec_src_file,
                                     wl_norm=wl_norm)
        else:
            self.spectrum = spectrum(date=self.windsoln.spec_date, wl_norm=wl_norm)
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
            self.spectrum.generate(kind=self.windsoln.spec_kind, savefile='inputs/spectrum.inp')
        return

    

    def make_string_from_scinum(self, x):
        exp = np.log10(x)
        base = x/(10**math.floor(exp))
        if base != int(base):
            return int(x), 0
        base = int(base)
        exp = int(np.log10(x/base))
        return base, exp

    
    
    
    def save_planet(self, name=None, folder='', overwrite=False,polish=False,
                    repo='saves/'):
        '''Note: It is not necessary to define a folder. All files will be saved in the saves/ repo.
        '''
        if polish == True:
            self.polish_bcs()
        if folder != '':
            if len(folder) > 1 and folder[-1] != '/':
                folder += '/'
            # Force all saves into saves directory
            if folder[0] != '/' and folder[0:len(repo)] != repo and folder[0]!='.':
                folder = repo+folder
            os.makedirs(folder, exist_ok=True)
        else:
            if name[0] != '/' and name[0]!='.' and name[0:len(repo)] != repo:
                folder = repo+folder
        if name is not None:
            if name[-4:] != '.csv':
                name = folder+name+'.csv'
            else:
                name = folder+name
        else:
            Mp, Rp, Mstar, semimajor, Ftot, Lstar = self.windsoln.planet_tuple
            name = (folder + ("Mp:{:g}_Rp:{:g}_Mstar:{:g}_a:{:g}_Ftot:{:g}_Nsp:{:d}_Spect:{:s}"
                             .format(Mp, Rp, Mstar, semimajor, Ftot, 
                                     self.windsoln.nspecies, self.windsoln.spectrum_tuple[3])) +
                    '.csv')
        if not overwrite and os.path.isfile(name):
            print('File already exists.\n'
                  '  To overwrite use save_planet(overwrite=True).')
            return
        print("Saving %s" % name)
        sub = Popen(["cp", 'saves/windsoln.csv', name],
                    stdout=PIPE, stderr=PIPE)
        sub.wait()
        output, error_output = sub.communicate()
        if error_output:
            print(error_output)
        return
    
    
    
    def easy_output_file(self,outputs=['v','T','rho'],
                         output_file='saves/simplified_outputs/output.dat',
                         comments='', overwrite=False):
        '''Description: Writes desired solution variables to csv that can be 
                        easily shared.

           Arguments: 
               outputs - list of str; desired solution columns. 'r' not required.
               output_file - str; default='output.dat'. Can include path.
               comments - str; default=''. Planet info always included in header
        '''
        if not overwrite and os.path.isfile(output_file):
            print('File already exists.\n'
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


    def run_wind(self, retry=False, expedite=False, verbose=False):
        # Run wind relaxation code
        if expedite == True: #do not integrate out, else use default of loaded solution
            self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)
#         print(self.windsoln.flags_tuple[3])
        sub = Popen('./bin/relaxed_ae', stdout=PIPE, stderr=PIPE) 
        output, error_output = sub.communicate()
        if error_output:
            if error_output[0:55] == b'Numerical Recipes run-time error...\nstep size underflow':
                return 4
            else:
                return 1
            if verbose:
                print(error_output)
        else:
            # If successful then update our guess to new solution
            sub = Popen(["cp", 'saves/windsoln.csv', 'inputs/guess.inp'],
                        stdout=PIPE, stderr=PIPE)
            output, error_output = sub.communicate()
            if error_output:
                print(error_output)
                return 2
            self.windsoln = wind_solution(expedite=expedite)
            if self.windsoln.error:
                # If we ran into analysis errors and not runtime errors, then
                # we probaby need to converge the bcs.
                if not retry:
                    # Check retry to avoid infinite recursion
                    print("Retry triggered")
                    self.ramp_base_bcs(expedite=expedite) #if expedite = False, converges bcs, too
                    self.run_isotherm(called_in_ramp_bcs=True)
                    self.converge_Ncol_sp(expedite=expedite)
                    return self.run_wind(retry=True)
                else:
                    print("Failed on retry to fix analysis errors")
                    return 3
            # saves/windsoln.csv and inputs/guess.inp are the same
            self.guess = self.windsoln
            return 0
                     

#Ramping Functions
    def ramp_to(self, system=None, converge=True,
                make_plot=False,integrate_out=True): #physics=None, bcs=None,
        if (system is None):# and physics is None):
            print("Please provide a system to ramp towards. system=system(Mp,Rp,Mstar,a,Ftot,Lstar).")
            return 0
        fail = 0
        result = 0
        if system is not None:
            self.ramp_class = "system"
            result = self.ramp_var("Ftot", system.value("Ftot"),
                                  converge_bcs=converge, make_plot=make_plot,
                                  expedite=True,
                                  integrate_out=False)
            fail = result
            
            if (fail != 0) and (fail != 5):
                return fail
            
            # Ramps Mp and Rp simultaneously, ~constant surface gravity
            result = self.ramp_grav(system, converge=converge,
                                   make_plot=make_plot, expedite=True,
                                   integrate_out=False)
            fail = result 
            
            if (fail != 0) and (fail != 5):
                return fail
            # Ramps Mstar and semimajor simultaneously, ~constant Hill radius
            result = self.ramp_star(system, converge=converge,
                                   make_plot=make_plot, expedite=True,
                                   integrate_out=False)
            fail = result
            
            
        if integrate_out == True:
            self.converge_Rmax()

        return fail


    def ramp_var(self, var, var_end, var_class=None, delta=0.02,
                 delta_additive=False, converge_bcs=True, make_plot=True,
                 expedite=False,integrate_out=True):
        if var_class is None:
            var_class = self.ramp_class
        if var_class == "system":
            var_val = self.system.value(var)
        elif var_class == "physics":
            var_val = self.physics.value(var)
        else:
            print("Unrecognized variable class for ramping: %s"
                  .format(var_class))
            return -100
        print(f'\rRamping {var:s} from {var_val:.3e} to {var_end:.3e}.',
              end='                                                     \n')
        if var_val == var_end or abs(var_val-var_end)/var_end < 1e-10:
            print(f'  {var:s} already done.',
                  end='                                                 \n')
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
            # See if can relax to partially ramped system
#             if converge_bcs:
#                 self.run_isotherm() #TEST - can afford to run each step
            result = self.run_wind(expedite=True)
            if result != 0:
                failed += 1
                delta /= 2.
                # Try converging. Can be expensive, so try only when smaller steps isn't working.
                if failed == 4: 
                    print("\r...Intermediate ramping BCs activated.",
                         end="                                           ")
#                     self.isotherm_start(run_wind=False) #updates where bolometric heat/cool dominate
                    self.ramp_base_bcs(expedite=True)
                    self.run_isotherm()
                    self.converge_Ncol_sp(expedite=True)
            elif result == 2:
                print("\nFailed to copy. Not good.")
                return 2
            else:
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge_bcs and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    
                    # converge boundary conditions on partial solution
                    self.ramp_base_bcs(intermediate=True,tolerance=0.1) #only converge if >10% diff in goal BC               
                    self.run_isotherm()
#                     self.converge_Ncol_sp() #MAYBE
                    conv_cntr = 0
                print("\r  Success %s:%.6e, delta:%.4g" %
                      (var, val_temp, flip*(val_temp-var_val)/var_val),
                      end="                                               ")
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
        print("\r  Final: ",
              end="\n")
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
        if expedite:
            expedite_flag_tuple = self.windsoln.flags_tuple
            expedite_flag_tuple[3] = 1
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out)
        result = self.run_wind(expedite=True, verbose=True)
#         if type(result) is not int: #a correct solution returns 0, failure when verbose=True returns [1,error]
#         while result != 0:
        fail = 0
        while result == 4:
            fail+=1
            print(f'\r Temporarily increasing Ncol_sp for numerical integration reasons. Failure {fail:d}',
                  end='                                                                                                                                                                              ')
            self.windsoln.bcs_tuple[5][:] = self.windsoln.Ncol_sp*1.2 
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            if fail > 10:
                print("\nERROR: Failed to integrate outwards too many times. ",
                      "\n  Recommendation: Manually raise total Ncol_sp and/or",
                      "\n                  Change numerical integrators (see docs)",
                      "\n                  Don't integrate outwards (integrate_out=False). ")
                return 1
            result = self.run_wind(expedite=True, verbose=True)
        else:
            if result == 0:
                if converge_bcs == True:
                    return self.polish_bcs(integrate_out) 
            else:
                print(f'\nERROR: Failure at the end game... result={result}')
                return 1

        return 0
    


    def ramp_grav(self, system, delta=0.02, converge=True, make_plot=True,
                  expedite=False,integrate_out=True):
#         self.failed_deeper_bcs_ramp = False
        var_Mp = self.system.value('Mp')
        var_Rp = self.system.value('Rp')
        var_sg = var_Mp/var_Rp**2
        srt_Mp = var_Mp
        srt_Rp = var_Rp
        srt_sg = var_sg
        end_Mp = system.value('Mp')
        end_Rp = system.value('Rp')
        end_sg = end_Mp/end_Rp**2
        print("\rRamping {:s} from {:.3e} to {:.3e} AND "
              "{:s} from {:.3e} to {:.3e}."
              .format('Mp', srt_Mp, end_Mp, 'Rp', srt_Rp, end_Rp), end='\n')
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
#             energy_plot(self.windsoln)
            self.ramp_base_bcs(intermediate=True,tolerance=0.1) #only converge if >10% diff in goal BC  
            self.run_isotherm()
#             print(self.run_isotherm(expedite=expedite)) #TEST
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
            elif var_Mp != end_Mp:
                # If we never needed to ramp Rp, ramp Mp linearly
                temp_Mp = var_Mp*(1.+M_delta)
                if (M_flip*temp_Mp >= M_flip*end_Mp):
                    temp_Mp = end_Mp
                temp.assign('Mp', temp_Mp)
                print("\r  Trying: {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mp', temp_Mp, M_flip*(temp_Mp-var_Mp)/var_Mp),
                      end="                                               ")
            elif var_Rp != end_Rp:
                # If we never needed to ramp Mp, ramp Rp linearly
                temp_Rp = var_Rp*(1.+R_delta)
                if (R_flip*temp_Rp >= R_flip*end_Rp):
                    temp_Rp = end_Rp
                temp.assign('Rp', temp_Rp)
                print("\r  Trying: {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Rp', temp_Rp, R_flip*(temp_Rp-var_Rp)/var_Rp),
                      end="                                               ")
            else:
                print("\nERROR: should be impossible surf_grav ramp condition")
                return
            # Write ramped variables to input file and try updating relaxation
            self.inputs.write_planet_params(*temp.system_tuple())
            result = self.run_wind(expedite=True)
            if result != 0:
                failed += 1
                if var_Mp != end_Mp:
                    M_delta /= 2.
                else:
                    R_delta /= 2.
                # Try converging. Can be expensive, so only try every once in a while
                if failed == 4:
                    print("...Intermediate Ramping BCs activated.")
#                     self.isotherm_start(run_wind=False) #updates where bolometric heat/cool dominate
                    self.ramp_base_bcs() 
                    self.run_isotherm()
#                     self.converge_Ncol_sp()
#                 self.converge_all_bcs(expedite=True) #not Ncol, because wasn't converged to before
            elif result == 2:
                print("\nERROR: Failed to copy windsoln to guess, not great.")
                return 2
            else: # Successfully relaxed with updated parameters
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    # converge boundary conditions on partial solution
#                     self.isotherm_start(run_wind=False) #updates where bolometric heat/cool dominate
                    self.ramp_base_bcs(intermediate=True,tolerance=0.1) #only converge if >10% diff in goal BC
                    self.run_isotherm()
#                     self.converge_Ncol_sp(expedite=expedite) #may be uneccesary
                    conv_cntr = 0
                # update our system to partially ramped system
                if (var_Mp != end_Mp) and (var_Rp != end_Rp):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mp', temp_Mp, 'Rp', temp_Rp,
                            M_flip*(temp_Mp-var_Mp)/var_Mp),
                          end="                                           ")
                    var_Mp = temp_Mp
                    var_Rp = temp_Rp
                    self.system.assign('Mp', var_Mp)
                    self.system.assign('Rp', var_Rp)
                elif var_Mp != end_Mp:
                    print("\r  Success %s:%.6e, M_delta:%.4g"
                          %('Mp', temp_Mp, M_flip*(temp_Mp-var_Mp)/var_Mp),
                          end="                                           ")
                    var_Mp = temp_Mp
                    self.system.assign('Mp', var_Mp)
                elif var_Rp != end_Rp:
                    print("\r  Success %s:%.6e, R_delta:%.4g"
                          %('Rp', temp_Rp, R_flip*(temp_Rp-var_Rp)/var_Rp),
                          end="                                           ")
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
                # Reset planet_params.inp to last sucessful state
                self.inputs.write_planet_params(*self.system.system_tuple())
                if expedite:
                    # return original flags
                    self.inputs.write_flags(*original_flag_tuple,integrate_out)
                print("\nERROR: Failing to converge gravity ramp.")
                return 101
        print("\r  Final: ")
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
        if expedite == True: #b/c was shut off before when expediting
            expedite_flag_tuple = self.windsoln.flags_tuple
            expedite_flag_tuple[2] = self.windsoln.bolo_heat_cool
            expedite_flag_tuple[3] = 1 #integrate_out will overwrite whatever user desires
            self.inputs.write_flags(*expedite_flag_tuple,integrate_out=integrate_out)
        result = self.run_wind(expedite=True, verbose=True)
        fail = 0
        while result == 4:
            fail+=1
            print(f'\r Temporarily increasing Ncol_sp for numerical integration reasons. Failure {fail:d}',
                  end='                                                                                                                                                                              ')
            self.windsoln.bcs_tuple[5][:] = self.windsoln.Ncol_sp*1.2 
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            if fail > 10:
                print("\nERROR: Failed to integrate outwards too many times. ",
                      "\n  Recommendation: Manually raise total Ncol_sp and/or",
                      "\n                  Change numerical integrators (see docs)",
                      "\n                  Don't integrate outwards (integrate_out=False). ")
                return 1
            result = self.run_wind(expedite=True, verbose=True)
        else:
            if result == 0:
                if converge == True:
                    return self.polish_bcs(integrate_out)
#                     if self.polish_bcs(integrate_out) == 0:
#                         return 0
#                     else:
#                         return 5
            else:
                print(f'\nERROR: Failure at the end game... result={result}')
                return 1

        return 0


    def ramp_star(self, system, delta=0.02, converge=True, make_plot=True,
                  expedite=False,integrate_out=True):
        '''Description: Ramps stellar mass (Mstar), semimajor axis (semimajor), and stellar bolometric 
                        luminosity (Lstar). If Mstar and semimajor axis change, ramps linearly along 
                        Hill radius rate of change.
                        If semimajor axis and Lstar change, ramps linearly along optical flux (essentially 
                        along skin temperature since T_skin propto Fopt) rate of change.
                        '''
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
                temp_Lstar = temp_Fopt*(temp_adist**2) #will this be going in the right direction? Yes, should be
                temp.assign('Lstar', temp_Lstar)
                print("\r  Trying: {:s}:{:.6e} & {:s}:{:.6e} &  {:s}:{:.6e}, M_delta:{:.4g}"
                      .format('Mstar', temp_Mstar, 'semimajor', temp_adist,'Lstar', temp_Lstar,
                              M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                      end="                                                                          ")
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
                self.ramp_base_bcs(intermediate=True,tolerance=0.1)
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
            else:
                print("\n ERROR: should be impossible surf_grav ramp condition")
                return 1
            # Write ramped variables to input file and try updating relaxation
            self.inputs.write_planet_params(*temp.system_tuple())
            result = self.run_wind(expedite=True,verbose=True)
            fail=0
            while result == 4:
                fail+=1
                print(f'\r Temporarily increasing Ncol_sp for numerical integration reasons. Failure {fail:d}',
                  end='                                                                          ')
                self.windsoln.bcs_tuple[5][:] = self.windsoln.Ncol_sp*1.2 
                self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                if fail > 10:
                    print("\nERROR: Failed to integrate outwards too many times. ",
                          "\n  Recommendation: Manually raise total Ncol_sp and/or",
                          "\n                   Change numerical integrators (see docs)",
                          "\n                   Don't integrate outwards (integrate_out=False). ")
                    return 1
                result = self.run_wind(expedite=True, verbose=True)

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
                    print(f"\r...Intermediate Ramping BCs activated (Instance {count:.0f}).",
                         end="                                                              ")
                    self.converge_Ncol_sp(expedite=expedite)
                    self.ramp_base_bcs() 
                    self.run_isotherm()
            elif result == 2:
                print("\nERROR: Failed to copy windsoln to guess. Not great.")
                return 2
            elif result == 0: # Successfully relaxed with updated parameters
                if make_plot:
                    plot_cntr += 1
                conv_cntr += 1
                failed -= 1
                if converge and (conv_cntr >= conv_every_n or
                                prcnt_chng >= conv_every_pc):
                    # converge boundary conditions on partial solution
                    # First update atmosphere to update Rmin location
                    self.run_isotherm() #TRY HERE OR RUNNING EVERY TIME
                    self.ramp_base_bcs(intermediate=True,tolerance=0.1) #only converge if >10% diff in goal BC
                    conv_cntr = 0
                # update our system to partially ramped system
                if (var_Mstar != end_Mstar) and (var_adist != end_adist) and (var_Lstar != end_Lstar):
                    print("\r  Success %s:%.6e & %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar, 'semimajor', temp_adist, 'Lstar', temp_Lstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                          end="                                                       ")
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
                    var_Mstar = temp_Mstar
                    var_adist = temp_adist
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('semimajor', var_adist)
                elif (var_Lstar != end_Lstar) and (var_adist != end_adist):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar, 'semimajor', temp_adist,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                          end="                                           ")
                    var_Mstar = temp_Mstar
                    var_adist = temp_adist
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('semimajor', var_adist)
                elif (var_Lstar != end_Lstar) and (var_Mstar != end_Mstar):
                    print("\r  Success %s:%.6e & %s:%.6e, M_delta:%.4g"
                          %('Lstar', temp_Lstar, 'Mstar', temp_Mstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                          end="                                                ")
                    var_Mstar = temp_Mstar
                    var_Lstar = temp_Lstar
                    self.system.assign('Mstar', var_Mstar)
                    self.system.assign('Lstar', var_Lstar)
                elif var_Mstar != end_Mstar:
                    print("\r  Success %s:%.6e, M_delta:%.4g"
                          %('Mstar', temp_Mstar,
                            M_flip*(temp_Mstar-var_Mstar)/var_Mstar),
                          end="                                           ")
                    var_Mstar = temp_Mstar
                    self.system.assign('Mstar', var_Mstar)
                elif var_adist != end_adist:
                    print("\r  Success %s:%.6e, R_delta:%.4g"
                          %('semimajor', temp_adist,
                            R_flip*(temp_adist-var_adist)/var_adist),
                          end="                                           ")
                    var_adist = temp_adist
                    self.system.assign('semimajor', var_adist)
                elif var_Lstar != end_Lstar:
                    print("\r  Success %s:%.6e, L_delta:%.4g"
                          %('Lstar', temp_Lstar,
                            L_flip*(temp_Lstar-var_Lstar)/var_Lstar),
                          end="                                           ")
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
                print("\nERROR: Failing to converge on star ramp.")
                return 101
        print("\r  Final: ")
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
        result = self.run_wind(expedite=True, verbose=True)
#         print(result)
        #a correct solution returns 0
        while result == 4:
            print(f'\r Temporarily increasing Ncol_sp for numerical integration reasons. {self.windsoln.Ncol_sp}',
                  end='                                                                                               ')
            self.windsoln.bcs_tuple[5][:] = self.windsoln.Ncol_sp*1.2 
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            if fail > 10:
                print("\nERROR: Failed to integrate outwards too many times. ",
                      "\n  Recommendation: Manually raise total Ncol_sp and/or",
                      "\n                   Change numerical integrators (see docs)",
                      "\n                   Don't integrate outwards (integrate_out=False). ")
                return 1
            result = self.run_wind(expedite=True, verbose=True)
        else:
            if result == 0:
                if converge:
                    return self.polish_bcs(integrate_out)
#                     if self.polish_bcs(integrate_out) != 0:
#     #                     print('Solution generated, but failed to polish solution BCs.')
#                         return 5
#                     else:
#                         return 0
                else:
                    return 0
            else:
                print(f'\nERROR: Failure at the end game... result={result}')
                return 1
    
    
#Metals functions    
    def metallicity(self,metals_list,Z):
        self.metals = metal_class(self.windsoln)
        return self.metals.metallicity(metals_list,Z)
        
        
        
    
    def add_metals(self, desired_species_list, Z=1, custom_mfs=[], Ncol_sp_tot=0, integrate_out=True):
        '''Description:
                Adds NEW metals to existing solutions. To ramp mass fractions or metallicity
                of species already in the solution use _______. Skips adding
                if species is already present. This function will converge to a self-
                consistent Ncol_sp for each species at the final step.
           Inputs:
               desired_species_list - list of str; list of ALL species in the self
                                      Give element + ionization state, e.g., 'HI',
                                      'h1', 'h I' for neutral hydrogen.
               Z            - float; metallicity. Default = 1.
               custom_mfs   - list or array; user specified mass fractions
               Ncol_sp_tot  - float; total column density at sonic point (1/cm2). 
                                     Default = 0 sets the Ncol_sp_tot to the current 
                                     value of the loaded solution. 
               integrate_out  - bool; can turn off integration for speed. 
                                      RECCOMMENDED True. Can occasionally be difficult 
                                      for Numerical Recipes ode integrator to integrate 
                                      out without a solution that is nearby.        
           Returns:
               None
        '''
        self.metals = metal_class(self.windsoln)
        
        desired_species = desired_species_list
        if Ncol_sp_tot == 0:
             Ncol_sp_tot = sum(self.windsoln.Ncol_sp)
        McAtom.formatting_species_list(desired_species)
        unspaced_list = [sp.replace(' ','') for sp in desired_species]

        #If user defined custom mass fractions, use those values
        if len(custom_mfs) > 0: 
            if len(custom_mfs) != len(desired_species):
                sys.exit("ERROR: Mass fraction and species list must be the same length.")
            if np.round(np.sum(custom_mfs),5) != 1:
                print('WARNING: Total Mass Fraction must sum to 1. ',
                      'sum(ZX) = %.3f' %np.sum(goal_mass_fracs))
            species_mf_dict = dict(zip(unspaced_list,custom_mfs))

        new_species = np.setdiff1d(unspaced_list,
                                   self.windsoln.species_list,
                                  assume_unique=True)
        for ns in new_species:
            self.metals = metal_class(self.windsoln)
            self.metals.add_species_to_guess(ns)
            if self.run_wind() != 0:
                print("Failed to add new species. That shouldn't have happened.",
                      "Raising Ncol_sp.")
#                 sys.exit(2)
                self.converge_Rmax()


            goal_mfs = self.metals.metallicity(self.windsoln.species_list,Z)
            current_mfs = self.windsoln.HX

            percent = 0.1
            if goal_mfs[-1]*10 >= 1:
                percent = 0.01
            delta = np.copy((goal_mfs - current_mfs)*percent)
            start_delta = np.copy((goal_mfs - current_mfs)*percent)
            

            #starting
            Ncols = np.copy(self.windsoln.Ncol_sp)
            Ncols[-1]*=10

            converge_Ncol = True
            first_run = True
            fail_counter = 0 
            while abs(sum((goal_mfs - self.windsoln.HX)/goal_mfs)) >1e-9:
                self.run_isotherm() #TEST
                diff = abs(sum((goal_mfs - self.windsoln.HX)/goal_mfs))
                print(f'\rAvg. mass fraction diff from goal: {diff:.2e}',
                      end='                                       ')

                #permanently take smaller steps if repeatedly need to lower stepsize
                if fail_counter >= 2:
                    delta = start_delta / 5*fail 
                #otherwise
                step_mfs = self.windsoln.HX + delta
                step_mfs[0] += 1-sum(step_mfs)

                #writing to input files
                self.inputs.write_physics_params(step_mfs,
                                                *self.windsoln.physics_tuple[1:])
                self.inputs.write_bcs(*self.windsoln.bcs_tuple[0:5],
                                       Ncols,
                                       self.windsoln.bcs_tuple[-1])

                converge_Ncol = True
                fail = 0
                start_mfs = np.copy(self.windsoln.HX)
                while self.run_wind(expedite=True) != 0:
                    fail+=1
                    if fail >= 2: #if repeatedly triggered, permanently take smaller steps
                        fail_counter += 1 
                    if fail > 10:
                        print(f'Failure to ramp {self.windsoln.species_list[-1]:s} to desired metallicity.')
                        print(f'Goal: {goal_mfs}, Current: {self.windsoln.HX}')
                        return 1
                    #if it will not solve, try converging Ncol one time
                    if (converge_Ncol == True) and (first_run == False):
                        self.converge_Ncol_sp(expedite=True,warning=False)
                        converge_Ncol = False
                    #otherwise, try taking smaller steps in mass fraction space
                    else:
                        print(f"\rFail {fail}: Trying delta / {(5*fail):d}.",
                              end="                                                  ")
                        step_mfs = start_mfs + delta/(5*fail)
                        step_mfs[0] += 1-sum(step_mfs)            
                        self.inputs.write_physics_params(step_mfs,
                                                *self.windsoln.physics_tuple[1:])

                #improves odds of ramping successfully if Ncol_sp is self consistently converged at some point
                if first_run == True:
                    self.converge_Ncol_sp(expedite=True,quiet=True)
                    first_run = False

                Ncols = self.windsoln.Ncol_sp 

            print(f"{ns:s} successfully added and ramped to mass fraction {goal_mfs[-1]:.2e}. Now converging Ncol_sp.")
            self.converge_Ncol_sp(quiet=True)
        return 0
    
    
    
    def remove_metals(self, remove_species_list, run_wind=True):
        self.metals = metal_class(self.windsoln)
        try:
            self.metals.remove_species_from_guess(remove_species_list)
#             self.windsoln.nspecies -= len(remove_species_list)
            self.load_planet('inputs/guess.inp',print_atmo_composition=False)
        except ValueError:
            print('One or more of the species you are attempting to remove is not present in simulation.')
            return 1
        if run_wind == True:
            if self.run_wind() != 0:
                warning = '''WARNING: Running sim with reduced number of species was unsuccessful.
                Consider reloading original solution and manually ramping down to 0 the mass fractions of
                the species you want to remove. Recall that total mass fraction must sum to 1.
                '''
                print(warning)
            else:
                print(','.join('%s' %sp for sp in remove_species_list)+" removed and new windsoln generated.")         
        return

    
    
    def ramp_metallicity(self,goal_Z=1,custom_mfs=[]):
        '''Description: Ramps up the metallicity of the species present in the simulation.
                        Can do in multiples of solar Z or set custom mass fractions.
           
           Arguments: 
               goal_Z - float; default=1. Multiples of Lodders (2008) solar metallicity.
               custom_mfs - list of floats; default=[]. If [], will default to goal_Z. 
                                           If not empty, will use the custom_mfs.
        '''
        if len(custom_mfs) != 0:
            if len(custom_mfs) != len(self.windsoln.species_list):
                sys.exit("ERROR: Mass fraction and species list must be the same length.")
            if np.round(np.sum(custom_mfs),5) != 1:
                print('WARNING: Total Mass Fraction must sum to 1. ',
                      'sum(ZX) = %.3f \n' %np.sum(custom_mfs))
            print("NOTE: Goal mass fractions will override any goal Z metallicity provided.")
            print(self.windsoln.species_list, goal_mass_fracs)
            current_HX = np.copy(self.windsoln.HX)
            self.inputs.write_physics_params(goal_mass_fracs,*self.windsoln.physics_tuple[1:])
            if self.run_wind(expedite=True) == 0:
                print(f'Mass fractions successfully ramped to: {self.windsoln.HX}')
                return self.converge_Ncol_sp(expedite=True)
            else:
                print('\r Taking smaller steps...')
                percent = 0.2
                delta = self.windsoln.HX - goal_mass_fracs
                while abs(sum((goal_mass_fracs - self.windsoln.HX)/self.windsoln.HX)) > 1e-6:
                    self.inputs.write_physics_params(self.windsoln.HX+percent*delta,
                                                         *self.windsoln.physics_tuple[1:])
                    if sum(abs(percent*delta)) > sum(abs(goal_mass_fracs - self.windsoln.HX)):
                        self.inputs.write_physics_params(self.windsoln.HX+(goal_mass_fracs - self.windsoln.HX),
                                                         *self.windsoln.physics_tuple[1:])
                    self.inputs.write_flags(integrate_out=False)
                    fail = 0
                    while self.run_wind() != 0:
                        fail+=1
                        step = self.windsoln.HX+percent*delta/(fail+1)
                        print(f'\r Fail {fail:d}: Attempting {step}',
                              end="                                           ")
                        self.inputs.write_physics_params(step,
                                                         *self.windsoln.physics_tuple[1:]) 
                        if fail>10:
                            print(f"Failed to ramp mass fractions. Current: {self.windsoln.HX}."
                                  f" Goal: {custom_mfs}")
                    self.converge_Ncol_sp(expedite=True)
                print(f'Mass fractions successfully ramped to: {self.windsoln.HX}')
                
        else: 
            grid = np.zeros(200)
            for i in range(200):
                grid[i] = abs(self.metallicity(self.windsoln.species_list,Z=i+1)[0]-
                              self.windsoln.HX[0])
            start_Z = np.where(grid==min(grid))[0][0]+1
            print("Starting metallicity: %d xSolar"%start_Z)
            current_Z = start_Z
            while (1-current_Z/goal_Z) > 1e-5:
                if goal_Z - current_Z > 5:
                    current_Z += 5
                elif goal_Z - current_Z < 5:
                    current_Z = goal_Z
                mfs = self.metals.metallicity(self.windsoln.species_list,Z=current_Z)

                self.inputs.write_physics_params(mfs,self.windsoln.species_list,
                                                 self.windsoln.molec_adjust)
                fail = 1
                while self.run_wind(expedite=True) != 0:
                    delta = goal_Z - current_Z
                    step_Z = current_Z + delta/(2**fail)
                    print(f'\r Failed: Attempting to ramp Z from {current_Z:.1f} to {step_Z:.1f}',
                         end='                                                            ')
                    if fail>10:
                        print('Failed at Z = ',current_Z)
                        sys.exit(1)
                    fail+=1
                print(f'\r Success! Attemping to ramp Z from {current_Z-2:.1f} to {current_Z:.1f}',
                 end='                                                                               ')
                stepsize=2
            print('Success! Ramped to goal metallicity, Z = %.0f x Solar'%current_Z)
            self.converge_Ncol_sp(expedite=True)
            return



#Polishing Boundary Conditions functions (many of these these enforce self-consistency, and are not neccessary for
#precision, but are for maximal accuracy (within the inherent uncertainty in the model))
    def polish_bcs(self,converge_Rmax=True):
        '''Description: 
        '''
        #Checking that base bcs (Rmin, rho, T) have been converged
        print('Polishing up boundary conditions...')
        self.ramp_base_bcs() 
        goal_bcs = self.base_bcs()
        curr_bcs = np.array(self.windsoln.bcs_tuple[:4])   
        avg_diff = sum(np.array((abs(curr_bcs-goal_bcs)/goal_bcs)[[0,2,3]]))

        #Checking that the molecular to atomic transition is occuring at the correct radius
#         isotherm = self.run_isotherm(polish=True)
        isotherm = self.run_isotherm() #temporary fix
        width = 1
        self.windsoln.add_user_vars()
        idx = np.searchsorted(self.windsoln.soln_norm['v'],
                              self.windsoln.erf_drop[0])
        while any(-self.windsoln.soln['heat_advect'][idx-10:idx+10]>
                  self.windsoln.soln['heat_ion'][idx-10:idx+10]):
#         while any(-self.windsoln.soln['heat_advect'] > self.windsoln.soln['heat_ion']):
            width+=1
            if width > 10:
                print("Warning: 20 scaleheights is an unlikely width for the error function",
                      " transitioning between molecular and atomic regions.",
                      "\nCheck energy_plot(). Stopping here.")
                break
            print(f"\r...Smoothing transition from molecular to atomic region. Erf width = {width:d} Hsc",
                 end="                                                                ")
            isotherm = self.run_isotherm(polish=True,width_factor=width/2)
            if isotherm != 0:
                print("Failed to smooth transition from molecular to atomic region.",
                      "Unphysical kinks in wind profile may be present at base of wind. ",
                      "\n       Mass loss rate relatively unaffected.")
                break
            self.windsoln.add_user_vars()
        
#         #if not capturing all of photoionization heating, print this warning
#         heat = self.windsoln.soln['heat_ion']
#         cool = self.windsoln.soln['cool_PdV']
#         try:
#             #finds where photoion heating first starts to dominate over PdV cool
#             idx = np.where(heat>-cool)[0][0]
#         except IndexError:
#             idx = 0
#         if idx <= 10:
#             P = self.windsoln.soln['P'][0]
#             note = f'''\nNOTE: Simulation currently does not capture all photoionization heating. Error in mass loss rate may be up to 10% in this case. \nTo set sim base deeper in wind, use ramp_base_bcs(user_override_press=True, base_press=__) where the base pressure is in units of microbars. This is an expensive calculation and may take a while - consider doing it on a H-He atmosphere. \n         Current sim base pressure: {10**np.floor(np.log10(P))} '''
#             print(note)
        #Converging Rmax (sets Rmax=r_cori) and/or Ncol_sp
        rcori_result = 'Success'
        if converge_Rmax==True:
            Ncol = self.converge_Ncol_sp()
            Rmax = self.converge_Rmax()
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
            bc_result = 'Failed'
        if isotherm != 0 and self.windsoln.bolo_heat_cool != 0:
            fails += 1
            warn += ' bolometric heating/cooling (run_isotherm),'
            iso_result = 'Failed'
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
            return 0
        else:
            if converge_Rmax == True:
                print(warn2 %(bc_result,iso_result,ncol_result,rcori_result))
            else:
                print(warn2 %(bc_result,iso_result,ncol_result))
            return 5
    def ramp_base_bcs(self,base_press=1,
                      user_override_press=False,
                      Kappa_opt=4e-3,Kappa_IR=1e-2,molec_adjust=2.3, 
                      adiabat=False, expedite=False, 
                      intermediate=False, tolerance=0.1):
        '''
        Description: Ramps Rmin, mass density at Rmin (rho_rmin), and temperature
                    at Rmin (T_rmin) to the physical values computed in base_bcs(). 
                    Default is the values at base_press = 1 microbar. 
                    If the loaded solution base BC is at higher pressure, though,
                    the base_press will stay that value in subsequent iterations.
        Arguments:
            base_press - int (or float); default=1. Accepted = multiples of 10.
                                        Desired pressure at base of sim in units
                                        of microbars. Default will be overridden
                                        if loaded solution has different base press,
                                        unless user_override_press=True.
            user_override_press - bool; default=False. To lower or raise pressure
                                        of base, set True.
            Kappa_opt - float; default = 4e-3. Optical opacity. Sets bolometric
                                               heating/cooling in molecular region
                                               below wind.
            Kappa_IR  - float; default = 1e-2. IR opacity. See above.
            molec_adjust - float; default=2.3. Dimensionless mean MOLECULAR weight.
                                  Mu of H2 is 2.3*mH.
            adiabat - bool; default=False. A mostly defunct feature that computes 
                            base bcs assuming atmo is adiabatic below wind.
            expedite - bool; default=False.
            
            
        '''        
        rho_scale = self.windsoln.scales_dict['rho']
        T_scale = self.windsoln.scales_dict['T']
        #Unless user wants to override the base pressure, take the pressure to be
        #the base pressure saved to solution
        if user_override_press == False:
            P = self.windsoln.soln['rho'][0]*const.kB*self.windsoln.soln['T'][0]/(self.windsoln.molec_adjust*const.mH)
            rounded_base = np.round(P/10)*10
            if rounded_base < 1:
                rounded_base = 1
            if rounded_base != base_press:
                base_press = rounded_base #=P would prob be more accurate in case someone sets it to 2
        
        goal_bcs = self.base_bcs(self.windsoln.molec_adjust,
                                 Kappa_opt,Kappa_IR,adiabat,base_press,
                                 user_override_press)

        curr_bcs = np.array(self.windsoln.bcs_tuple[:4])        
        avg_diff = np.array((abs(curr_bcs-goal_bcs)/goal_bcs)[[0,2,3]])
        if (intermediate==True) and all(avg_diff < tolerance):
            return
        
        if sum(abs((goal_bcs - curr_bcs)/goal_bcs)[[0,2,3]]) > 1e-2:
            print(f"\rAttempting to ramp Rmin:{curr_bcs[0]:.2f}->{goal_bcs[0]:.2f}Rp, rho:{curr_bcs[2]*rho_scale:.3e}->{goal_bcs[2]*rho_scale:.3e}g/cm3, T:{curr_bcs[3]*T_scale:.0f}->{goal_bcs[3]*T_scale:.0f}K",
                  end="                                                                                                                                                ")
            self.inputs.write_bcs(*goal_bcs,*self.windsoln.bcs_tuple[4:])

            if self.run_wind() == 1: 
                print(f"..Initial jump failed. Ramping variables individually.")
                if abs((self.windsoln.Rmin - goal_bcs[0])/goal_bcs[0]) > 4e-3:
                    self.ramp_Rmin(goal_bcs[0])
                if abs((self.windsoln.T_rmin  - goal_bcs[3])/goal_bcs[3]) > 1e-4:
                    self.ramp_T_rmin(goal_bcs[3])
                if abs((self.windsoln.rho_rmin -  goal_bcs[2])/goal_bcs[2]) > 1e-4:
                    self.ramp_rho_rmin(goal_bcs[2])

        curr_bcs = np.array(self.windsoln.bcs_tuple[0:4])
        if sum(abs((goal_bcs - curr_bcs)/goal_bcs)[[0,2,3]]) > 1e-2: 
            print("Failed to ramp bcs. This can occur when attempting to set base at too high or too low of pressure (e.g., for very high or very low gravity planets).")
            print(f"  Rmin: Current {curr_bcs[0]:.2f}, Goal {goal_bcs[0]:.2f}Rp; Rho: Curr {curr_bcs[2]*rho_scale:.3e}, Goal {goal_bcs[2]*rho_scale:.3e}g/cm3; T: Curr {curr_bcs[3]*T_scale:.0f}, Goal {goal_bcs[3]*T_scale:.0f}K")
            return 1
        else:
            print("\rSuccessfully ramped base boundary conditions.",
                  end="                                                           " )            
            return 0
       
    
    
    
    def base_bcs(self,molec_adjust=2.3,Kappa_opt=4e-3,
                 Kappa_IR=1e-2,adiabat=False,
                 base_press=1,user_override_press=False): 
        ''' Sets base of simulation below XUV tau=1 surface and assumes bolometric 
        heating & cooling of molecules dominates there.

        Computes the density and temperature at either R_ubar (microbar 
        pressure radius) or R_IR (*vertical* tau=1 surface to IR radiation 
        - between this radius and the nanobar radius where wind is launched, 
        the temperature profile is an isotherm at the skin temperature). 
        If R_IR is below Rp, it takes R_mbar to be the base of the simulation.

        Assumes that Rp is the *slant path* optical tau = 1 surface.      
        Computes skin temperature by balancing bolometric heating and cooling.
        Computes effective temperature from stellar bolometric luminosity. 
        For highly irradiated planets (most planets for which atmo escape is 
        significant), it is an isotherm at T_eff=L_star/(16*pi*sigma_sb*a^2)
        between R_IR and Rp or an isotherm at T_skin between R_mbar and Rp. 
        If adiabat = True, it computes the R_IR assuming an adiabat between 
        T_eff and T_skin.

        For high metallicity atmospheres, users should increase the molecular
        adjustment factor (molec_adjust). Likewise, they should change Kappa_opt
        and Kappa_IR.

        Parameters:
            molec_adjust - float; default=2.3. Accounts for increased mean weight
                                               of molecules
            Kappa_opt    - float; default=4e-3. Optical opacity
            Kappa_IR     - float; default=1e-2. IR opacity
            adiabat      - Bool; default=False. If true, computes IR BCs assuming
                                                an adiabat

        Returns:
             R - float; radius at base of simulation in units of Rp
             Rmax - float; set by original boundary conditions, units of Rp
             rho - float; density in units of RHO0 from defs.h 
                         (1e-15 ~nBAR*MH/K/T0)
             T - float; effective temperature in units of T0 (1e4 K)


        If R_IR > Rp:
            Returns: R_IR, Rmax, rho_IR, T_eff
        Else:
            Returns: Rp, Rmax, rho_Rp, T_skin
        '''
#         #computing current pressure at base. If current base pressure > 1, set current as base pressure
#         #this is because some high flux planets may need lower lower BC to capture all heating

        if user_override_press == False:
            P = self.windsoln.soln['rho'][0]*const.kB*self.windsoln.soln['T'][0]/(self.windsoln.molec_adjust*const.mH)
            rounded_base = np.round(P/10)*10
            if rounded_base < 1:
                rounded_base = 1
            if rounded_base != base_press:
                base_press = rounded_base #=P would prob be more accurate in case someone sets it to 2
        
        Rp = self.windsoln.Rp
        if molec_adjust <= 1:
            print("WARNING: Molecular adjustment factor should be >1 to account for increased",
                  "mean weight due to molecules, instead of atoms, lower in atmosphere (below the wind).")
        F_opt = self.windsoln.Lstar/(4*np.pi*self.windsoln.semimajor**2)
        T_skin = (F_opt*(Kappa_opt+Kappa_IR/4)/(2*const.sig_SB*Kappa_IR))**0.25
#         print(T_skin)
        T_eff = (F_opt/(4*const.sig_SB))**0.25
        mu = molec_adjust*const.mH/np.sum(self.windsoln.atomic_masses[0]*self.windsoln.HX/self.windsoln.atomic_masses) #assuming no ionized
        #rho_Rp derived from optical slant path tau=1 geometry
        rho_Rp = np.sqrt(mu*const.G*self.windsoln.Mp / (8 * Rp**3 *const.kB*T_skin)) / Kappa_opt
        P_Rp = rho_Rp*const.kB*T_skin / mu #in barye = 1e-6 bar

        #take solution base as millibar pressure (=1000 barye), to ease computational burden
        cs2 = (const.kB*T_skin/(2.3*const.mH))
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
            R_IR = Rp**2 /( Hsc * np.log(rho_IR/rho_Rp*np.exp(Rp/Hsc)) )
        if R_IR > Rp:
            print('Using vertical tau=1 for IR photons, R_IR.')
            return R_IR/Rp,self.windsoln.Rmax,rho_IR/1e-15,T_eff/1e4
        else:
            return R_mbar/Rp,self.windsoln.Rmax,rho_mbar/1e-15,T_skin/1e4

        
        
    def converge_Rmax(self,final_converge_Ncol=True):
        """
        Description:
            Self-consistently sets Rmax to be the Coriolis length. 
            Does not worry about converging other boundaries. Cannot 
            be expedited as we specifically need to have the Coriolis 
            length calculated.
            
        Arguments:
            final_converge_Ncol - bool; default=True. Converges Ncol_sp self-
                                        consistently again after Rmax has been 
                                        set to Rcori.
        """
#         print(self.windsoln.R_cori)
#         self.run_wind()#need to populate saves/windsoln.csv with windsoln
        #otherwise could cause trouble if new solution is and old windsoln.csv
        #is loaded 
        try: # Check if r_cori has be calculated
            self.windsoln.R_cori
        except AttributeError:
            if self.windsoln.integrate_outward == 0:
                self.inputs.write_flags(*self.windsoln.flags_tuple,
                                        integrate_out=True)
                out = self.run_wind(expedite=False)
                attempt = 0
                while out == 4:
                    attempt += 1
                    self.windsoln.bcs_tuple[-2] = self.windsoln.Ncol_sp*10*attempt
                    print(f"\rRaising Ncol to allow outward integration. Attempt {attempt:d}.",
                          end="                                                            ")
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                    if attempt > 6:
                        print("Failed to integrate outwards even after raising Ncol. Suspicious.")
                        return 1
                    out = self.run_wind()
                else:
                    if out == 1:
                        print("Relaxation error.")
                        return 1
                    elif out == 0:
    #                     self.load_nonexpedited('saves/windsoln.csv')
                        print("\rSuccessfully integrated outwards. Now converging Rmax...",
                              end="                                                                                                                                 ")
                        self.load_nonexpedited('inputs/guess.inp') #shouldn't it be windsoln?
            else:
                self.windsoln.add_user_vars()
    
        # First while statement extends domain far past true Coriolis length
        while (self.windsoln.Rmax < self.windsoln.R_cori):
            print("\r..Rmax {:.4e}, r_Cori {:.4e}"
                  .format(self.windsoln.Rmax, self.windsoln.R_cori),
                  end="                                                   ")
            bcs_tuple = self.windsoln.bcs_tuple
            # Just go far enough to ensure true R_cori is within domain
            bcs_tuple[1] = min(self.windsoln.R_cori+10,
                               1.5*self.windsoln.R_cori)
#             print(bcs_tuples)
            self.inputs.write_bcs(*bcs_tuple)
            failed = 0
            if self.run_wind(expedite=False) != 0:
                print("\nERROR: Failed to integrate outwards to r_Cori: {:.3e} "
                      "reverting back to {:.3e}."
                      .format(self.windsoln.R_cori, self.windsoln.Rmax))
                bcs_tuple[1] = self.windsoln.Rmax
                self.inputs.write_bcs(*bcs_tuple)
                return 1
        # Second while statement whittles down to true Coriolis length
        while (abs(1.-self.windsoln.Rmax/self.windsoln.R_cori) > 1e-2):
            print("\r..Rmax {:.4e}, r_Cori {:.4e}"
                  .format(self.windsoln.Rmax, self.windsoln.R_cori),
                  end="                                                   ")
            delta_Rmax = self.windsoln.Rmax-self.windsoln.R_cori
            bcs_tuple = self.windsoln.bcs_tuple
            bcs_tuple[1] = self.windsoln.Rmax-0.99*delta_Rmax
            self.inputs.write_bcs(*bcs_tuple)
            failed = self.run_wind(expedite=False)
            if failed:
                print('\nERROR: Failed to even get relax code to run')
                return 1
        if final_converge_Ncol == True:
            current,goal = self.self_consistent_Ncol(warning=False)
            self.windsoln.bcs_tuple[5] = goal
            self.inputs.write_bcs(*self.windsoln.bcs_tuple) 
            result = self.run_wind()
            fail=0
            while result != 0:
                if result == 4:
                    print("Rmax done, but Ncol cannot be lowered to self-consistent value without introducing integration errors. Maximum error in Mdot in this case is 10 percent.")
                    return 1
                else:
                    fail += 1
                    self.windsoln.bcs_tuple[5] = current+(goal-current)/(fail+1)
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple) 
                    result = self.run_wind()
            print(f"\rSuccessfully converged Rmax to {self.windsoln.Rmax:.6f} (Ncol also converged).",
                  end="                                                       ")
        else:
            print(f"\rSuccessfully converged Rmax to {self.windsoln.Rmax:.6f}. It is recommended to run converge_Ncol_sp() as Ncol will have updated.",
                  end="                                                       ")
        return 0
    

    def converge_Ncol_sp(self,expedite=False,quiet=False):
        """
        Description:
            Self-consistently converges Ncol at the sonic point, such that the 
            column density boundary condition for each species, Ncol_sp,
            matches the number density at the sonic point.

        Inputs:
            expedite - bool; default = False. True does not integrate past Sonic Point.
            quiet - bool; default = False. Prints warnings about Ncol_sp being only estimates if
                                            integrate_out = False.
        """
        scale = self.windsoln.scales_dict['Ncol_HI'] #scales the same across species
        #self-consistent Ncol finder
        Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False)
        Ncol_goal[Ncol_goal==0] = 1e-10*Ncol_goal[0]
#         print(Ncol_goal)
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
                print(f'Too many iterations (current average diff {avg:.2e}).')
                return 1
            avg = np.abs(np.mean((Ncol_goal-Ncol_current)/Ncol_goal))
            print(f'\r ...Ramping Ncol Iter {niter:d}: Current average diff {avg:.2e}', 
                  end='                                                          ')

            self.windsoln.bcs_tuple[5][:] = Ncol_goal  
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            fail = 1
            result = self.run_wind(expedite=expedite)
            while result != 0:
                attempt = 0
                if result == 4:
                    print("\rIntegration error. Turning off outward integration.",
                         end="                                                          ")
                    self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)
                    result = self.run_wind(expedite=expedite)
                else:
                    delta = (Ncol_goal - Ncol_current)/(2*fail)
                    Ncol_step = Ncol_current + delta
                    self.windsoln.bcs_tuple[5][:] = Ncol_step 
                    self.inputs.write_bcs(*self.windsoln.bcs_tuple)
                    print(f'\r...Fail {fail:d}: Attempting Ncol_sp = { Ncol_step*scale}.',
                          end='                                                                        ')

                    if fail > 10:
                        print(f"Failed at Ncol_sp = {Ncol_current} ",
                              f"& Summed neutral number density = {Ncol_goal}.")
                        return 2
                    result = self.run_wind(expedite=expedite)
                    fail += 1

            Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False) 
            Ncol_current = np.copy(self.windsoln.Ncol_sp)

        Ncol_current,Ncol_goal = self.self_consistent_Ncol(warning=False) 
        avg = np.abs(np.mean((Ncol_goal-Ncol_current)/Ncol_goal))
        if quiet == False:
            print(f'\r Success! Average Ncol difference: {avg:.2e}', 
                  end='                                                   ')

        #Final converence of all species at once (MAY NOT WANT THIS BECAUSE WILL MOVE IT again)
        if niter != 0:
            if expedite == True:
                if quiet == False:
                    print('\r...Attempting Final Ncol Convergence. (Note: This is an estimate. Cannot converge precisely without outward integration).',
                          end='                                                              ')  
            if expedite == False:
                if quiet == False:
                    print('\r...Attempting Final Ncol Convergence.',
                          end='                                                              ')  

            self.windsoln.bcs_tuple[5][:] = Ncol_goal
            self.inputs.write_bcs(*self.windsoln.bcs_tuple)
            if self.run_wind(expedite=expedite) == 4:
                print('Failure at last Ncol convergence.',
                      f' Last working solution should be fine for most purposes. Average diff: {avg:.2e}')
                return 6
        return 0

    
    
    def self_consistent_Ncol(self,method=1,warning=True):
        '''Description: Computes the self-consistent column density sonic point boundary condition
                        from the neutral number density of a given species.
           
           Returns:
               current - array of current Ncol_sp values
               goals - array of self-consistent Ncol_sp values
        '''
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
#             n_neutral *= self.windsoln.soln['Ys_'+self.windsoln.species_list[j]]
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
                print("This is an estimate of Ncol. Cannot converge precisely without outward integration.") 
            return currents, goals#*8 #to account for the fact that we can't cumsum n out to Rmax
        if self.windsoln.flags_tuple[-1] == 1:
            return currents, goals 
        
        
        
    def run_isotherm(self,width_factor=1.,return_idx=False,called_in_ramp_bcs=False,
                     polish=False):
        '''Description: Ramps the complementary error function that governs the 
                        transition from the molecular to atomic regions to the atmosphere.
                        Criterion for transition: when photoionization heating begins to
                        dominate over the PdV cooling and a wind will launch.
                        
                        Below the wind, for an average planet, molecules have not 
                        photodissociated so mu should be mean molecular weight instead of
                        mean atomic weight and the erf enforces this. Molec. opacities mean 
                        that bolometric heating and cooling dominate the energy budget and 
                        create an isotherm in molecular region below the wind. In the optically 
                        thin atomic wind, bolometric heating and cooling are negligible, so the 
                        erf also enforces the drop off bolometric heating and cooling.
                        
                        If a user needs to call this function for some reason, width_factor is the
                        only argument they should need to change. For numerical or physical reasons
                        it may be neccessary for the transition to occur over more than 1 scaleheight.
                        In this case, increase width_factor.
                        
           Arguments: 
               width_factor - float; default = 1. Sets the width in scaleheights over which the 
                                                  erf drops off. 
               return_idx - bool; default=False.
               called_in_ramp_bcs - bool; default=False.
               polish - bool; default=False
               
        '''
        if (self.windsoln.bolo_heat_cool == 0) and (polish==False):
            return
        if (polish==False) and (self.skip==True):
            return
        
        #TEST
        #loading last working solution to do this. Risky.
        self.load_planet('saves/windsoln.csv',print_atmo_composition=False)
        v_drop,rate = self.erf_velocity(called_in_ramp_bcs=called_in_ramp_bcs,
                                        polish=polish,
                                        width_factor=width_factor)    
       
        diffs = np.array([abs(v_drop-self.windsoln.erf_drop[0])/v_drop,
                         abs(rate-self.windsoln.erf_drop[1])/rate])
        if any(diffs) > 1e-5:
            current_erfs =  np.array([self.windsoln.erf_drop[0],self.windsoln.erf_drop[1]])
            goal_erfs    = np.array([v_drop,rate])
            if (self.windsoln.bolo_heat_cool == 1):
                if (self.skip == False) or (polish == True):
                    while np.mean(abs(goal_erfs - current_erfs)/goal_erfs) > 1e-5:
                        delta = goal_erfs - current_erfs
                        fail = 0
                        self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                              current_erfs+delta)
                        while self.run_wind(expedite=True) != 0:
                            fail += 1
                            if fail >= 10:
                                current_erfs = self.windsoln.erf_drop
                                print(f"Too many attempts to ramp error function. Skipping subsequent attempts until polishing. Goal: {goal_erfs[0]:.3e}, {goal_erfs[1]:.3e}.",
                                      f" Current: {current_erfs[0]:.3e}, {current_erfs[1]:.3e}.")
                                self.skip = True
                                return #self.turn_off_bolo()
                            elif fail <= 2: #try smoothing out erf
                                step_erfs = np.copy(current_erfs)
                                step_erfs[1] *= 5*fail
                                self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                      step_erfs)
                                print(f"\r Erf Fail {fail:d}: Smoothing transition: {step_erfs[1]}",
                                      end="                                      ")
                            elif (fail>2) & (fail<10) : #try taking a smaller step
                                step_erfs = current_erfs+delta/(fail+1)
                                print(f"\r Erf Fail {fail:d}: goal {goal_erfs[0]:.3e} {goal_erfs[1]:.3e}, step {step_erfs[0]:.3e} {step_erfs[1]:.3e}",
                                      end="                                      ")
                                self.inputs.write_bcs(*self.windsoln.bcs_tuple[:6],
                                                      step_erfs)
                        current_erfs = self.windsoln.erf_drop
                    return 0
                else:
                    return 1
            #on final polishing run of isotherm, check if bolometric heating/
            #cooling should be added back in
            if (self.windsoln.bolo_heat_cool == 0) and (polish==True):
#                 self.load_planet('inputs/guess.inp',expedite=False)
                self.windsoln.add_user_vars()
                heat = self.windsoln.soln['heat_ion']
                cool = self.windsoln.soln['cool_PdV']
#                 heat,cool = self.quick_calc_heat_cool()
                if len(np.where(-cool[:20]<heat[:20])[0]) > 3:
                    print("NOTE: Photoionization heating dominates down to base",
                          f"of sim ({self.windsoln.Rmin:.3f} Rp).",
                          "\n      Bolometric heat/cooling still turned off. Max error in dM/dt ~ 10%. \n       See documentation for workaround.")
                    return 1
                else: #if appropriate, ramp back in bolometric heating/cooling
                    while self.windsoln.bolo_heat_cool < 1:
                        flags = self.windsoln.flags_tuple
                        delta = 1-self.windsoln.bolo_heat_cool
                        flags[2] += delta
                        self.inputs.write_flags(*flags)
                        fail = 0
                        while self.run_wind() != 0:
                            fail+=1
                            flags[2] = self.windsoln.bolo_heat_cool + 0.1/fail
                            self.inputs.write_flags(*flags)
                            if fail>10:
                                print("Warning: Bolometric heating/cooling failed to ramp back in.")
                                energy_plot(self.windsoln)
                                return 1
                    print("  Bolometric heating/cooling successfully ramped back in.")
                    return 0
                            
            else:
                return 1
        else:
            return 0
        
        
    def quick_calc_heat_cool(self):
        ''' Description: Inexpensively calculates the photoionization heating and PdV cooling as 
                         a function of radius. Used to compute the arguments in the complementary
                         error function that enforces the molecular to atomic transition.
                
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
    #     print(frac_in_heat)
        for s, species in enumerate(self.windsoln.species): #48s
            E_matrix = np.tile(self.windsoln.E_wl,(len(background_ioniz_frac),1))
            Ncol_matrix = np.tile(Ncol_arr.iloc[:,s],(self.windsoln.npts,1)).T
            sig_matrix = np.tile(sigmas.iloc[:,s],(len(background_ioniz_frac),1))

            f = np.nan_to_num(sig_matrix*Ncol_matrix / taus) #frac of incoming photon energy that will interact with species s
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
        self.load_planet('saves/windsoln.csv',print_atmo_composition=False)
#         if failed_bolo_turn_off == True:
#             print("Previously failed to turn off bolometric heating and cooling, so not trying here.")
#             return
        if self.windsoln.bolo_heat_cool == 1:
            print('..Turning off bolometric heating/cooling.')
        #turning off bolo_heat_cool
        while self.windsoln.bolo_heat_cool > 0:
            flags = self.windsoln.flags_tuple
            bolo_flag = np.copy(flags[2])
            delta = -flags[2]
            flags[2] = 0 #turning off bolo heating and cooling
            self.inputs.write_flags(*flags)
            fail = 0
            while self.run_wind(expedite=True) != 0:
                fail+=1
                if fail>5:
                    print(f"WARNING: Failed to turn off bolometric heating/cooling. This is unusual. Current multiplicative factor: {self.windsoln.bolo_heat_cool}. Goal: 0.")
                    self.windsoln.flags_tuple[2] = 1
                    self.inputs.write_flags(*self.windsoln.flags_tuple)
                    return 1
                delta /= 2
                flags[2] = bolo_flag + delta
                print(f'\rTurning off bolometric heating & cooling. Trying factor of {flags[2]:.2f}',
                      end='                                                  ')
                self.inputs.write_flags(*flags) 
                
        return 0
    
    
        
    def erf_velocity(self,return_idx=False,called_in_ramp_bcs=False, 
                     polish=False,width_factor=1):
        '''Description: Defines the drop-off radius at which the complementary 
                        error function that governs the drop off of the bolometric
                        heating and cooling and the drop off of the mean molecular weight
                        in the isothermal part of the wind as photoionization
                        heating begins to dominate and the atmosphere becomes
                        atomic and non-isothermal.
                        
           Arguments: return_idx - bool; default=False. If true, returns radial
                                         index where drop-off occurs.
                      called_in_ramp_bcs - bool; avoids recursion of ramping bases
                                                  if called in ramp_base_bcs()
           Returns: v[drop_index] - velocity of drop-off in units of cm/s.
        '''
                    
        
        #now turns off boloheat/cool if sim base not deep enough
        heat, cool = self.quick_calc_heat_cool()
        if polish == True:
            self.windsoln.add_user_vars()
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
        def rate_calc(idx):
            #Approximate pressure scaleheight at drop radius in units of Rp
            Hsc =  const.kB*self.windsoln.soln['T'][idx]*self.windsoln.soln['r'][idx]
            Hsc /= (const.mH*self.windsoln.molec_adjust*const.G*self.windsoln.Mp) 
            #approximates mu as molecular value
            #Compute the desired rate of erfc drop-off as 
            if idx >= 10:
                slope = (v[idx+10] - v[idx-10])/(r[idx+10] - r[idx-10])
            else:
                slope = (v[idx+10] - v[idx])/(r[idx+10] - r[idx])
            width = width_factor*Hsc 
            rate = slope*width
            return rate
            
        #If the base of the simulation is still within the photoionization heating
        #region, try to lower the base of the simulation.
#         if (drop_index<=10) and (called_in_ramp_bcs==False) and (self.failed_deeper_bcs_ramp==False) and (expedite==False):
#             rmin = self.windsoln.Rmin
#             print(f"Current sim base ({rmin:.3f}Rp) may not capture all",
#                   " photoion heating.\n",
#                  "   Setting lower BC radius deeper in atmosphere.")
#             #while the photoionization doesn't drop below PdV cooling, 
#             #lower base BCs x2 attempts
#             fail=0
#             while len(np.where(-cool[:20]>heat[:20])[0]) < 3:
#                 P =  self.windsoln.soln['rho'][0]*const.kB
#                 P *= self.windsoln.soln['T'][0]/(self.windsoln.molec_adjust*const.mH)
#                 OOM_base_press_old = 10**np.floor(np.log10(P))
#                 if OOM_base_press_old >= 100:
#                     print("Lowering base past 100 microbars unlikely to lower",
#                           " Rmin sufficiently to capture all photoionization",
#                           " heating for this high gravity planet.")
#                     print("  Bolometric heating/cooling turned off.")
#                     self.failed_deeper_bcs_ramp = True #don't waste time trying
#                     self.turn_off_bolo()
#                     return v[0],rate_clac(0)
                
#                 #don't try to lower more than 2 times b/c very expensive
#                 if fail>=2:
#                     print(f"Current base P: {OOM_base_press_old:.0f} microbars.",
#                           "\nStopping here, but user can ramp manually",
#                           " to higher pressure.\n  Still not capturing all ",
#                           "photoion heating, so bolometric heating/cooling",
#                           " turned off.")
#                     self.failed_deeper_bcs_ramp = True #don't waste time trying
#                     self.turn_off_bolo()
#                     if return_idx==True:
#                         return v[0],rate_calc(0),0
#                     else:
#                         return v[0],rate_calc(0)
                
#                 #ramping to 10xlower base pressure
#                 result = self.ramp_base_bcs(base_press=OOM_base_press_old*10,
#                                           user_override_press=True)
#                 if result == 0:
#                     heat, cool = self.quick_calc_heat_cool()
#                 else:
#                     print("Failed to ramp BCs lower in atmo.")
#                     self.failed_deeper_bcs_ramp = True #don't waste time trying
#                     self.turn_off_bolo()
                    
#                     if return_idx==False:
#                         return v[0],rate_calc(0)
#                     else:
#                         return v[0],rate_calc(0),0
#                 fail+=1
        
        #If this is being called inside of the ramping function, do this 
        #to avoid recursion
        if (drop_index<=10) and (called_in_ramp_bcs==True):
            if return_idx==False:
                return v[0],rate_calc(0)
            else:
                return v[0],rate_calc(0),0
        
        #If it has previously failed to ramp deeper 
        elif (drop_index<=10):# and (self.failed_deeper_bcs_ramp==True):
            self.turn_off_bolo()
            if return_idx==False:
                return v[0],rate_calc(0)
            else:
                return v[0],rate_calc(0),0
            
#         #If expediting, should skip ramping base deeper and just turn off bolo.    
#         elif (drop_index<=10) and (expedite==True):
# #             self.failed_deeper_bcs_ramp = True #don't waste time trying
#             self.turn_off_bolo()
#             if return_idx==False:
#                 return v[0],rate_calc(0)
#             else:
#                 return v[0],rate_calc(0),0
            
        if drop_index > 500: #low grav superearths need no drop off
            rate = 10*rate_calc(drop_index)
        else:
            rate = rate_calc(drop_index)
        #If none of the above conditions are triggered, simply return v
        if return_idx == True:
            return v[drop_index],rate,drop_index
        else:
            return v[drop_index],rate   
        

        

        
    def ramp_T_rmin(self, goal_T):
        """
        Description:
            Ramps normalized temperature at Rmin. 
        Inputs:
            goal_T - float; goal temperature at Rmin in units of 1e4 K
        """
        T_scale = self.windsoln.scales_dict['T']
        if goal_T > 1:
            print("WARNING: T should be in units of 1e4 K. Returning...")
            return 
        while (abs(1.-self.windsoln.T_rmin/goal_T) > 1e-10):
            bcs_tuple = self.windsoln.bcs_tuple
            bcs_tuple[3] = goal_T
            self.inputs.write_bcs(*bcs_tuple)
            failed = 0
            print(f'\r..T_rmin {self.windsoln.T_rmin*T_scale:.0f}, goal {goal_T*T_scale:.0f}',
                  end="                                                   ")
            while self.run_wind(expedite=True) != 0:
                failed += 1
                delta_T = (bcs_tuple[3]-self.windsoln.T_rmin)/2
                bcs_tuple[3] = self.windsoln.T_rmin+delta_T
                self.inputs.write_bcs(*bcs_tuple)
                print(f'\r..T_rmin {self.windsoln.T_rmin*T_scale:.0f}, '
                      f'try {bcs_tuple[3]*T_scale:.0f}',
                      end="                                               ")
                if failed > 15:
                    print("\nStruggling to substep towards new T_rmin")
                    return 1
        print(f"\r   Successfully converged T_rmin to {self.windsoln.T_rmin*T_scale:.0f}",
              end="                                                       ")
        self.run_isotherm(called_in_ramp_bcs=True) #TEST
        return 0
    
    
    
    def ramp_Rmin(self,goal_Rmin):
        """ Description:
                Ramps normalized Rmin, where Rmin is the base of the simulation and 
                should ideally be a radius "below the wind", a.k.a., in the molecular region
                below ~1 microbar. 
            Arguments:
                The desired Rmin value in normalized units of Rp.
        """
        if goal_Rmin > 30:
            print("WARNING: Rmin should be in units of Rp. Returning...")
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
            print(f'\r..Goal Rmin: {goal_Rmin:.4f}Rp. Current: {self.windsoln.Rmin:.4f}Rp. Trying {step_Rmin:.4f}',
                  end="                                                                                           ")
            fail = 0
            result = self.run_wind(expedite=True)
            while result != 0:
                fail+=1
                if fail==1:
                    self.run_isotherm(called_in_ramp_bcs=True) #TEST
                if result == 4:
                    self.inputs.write_flags(*self.windsoln.flags_tuple,
                                            integrate_out=False)
                    print("Turning off outward integration temporarily.")
                else:
                    bcs_tuple = self.windsoln.bcs_tuple
                    step_Rmin = current_Rmin + delta/(fail+1)
                    bcs_tuple[0] = step_Rmin
                    self.inputs.write_planet_params(*self.windsoln.planet_tuple)
                    print(f'\r..Fail {fail:d}: Rmin {current_Rmin:.5g}Rp, trying {step_Rmin:.5g}Rp',
                          end="                                                                                                                                    ")
                    self.inputs.write_bcs(*bcs_tuple)
                if fail > 10:
                    print(f"   Failed at {current_Rmin} Rp. Goal {goal_Rmin} Rp.")
                    return 1
                result = self.run_wind(expedite=True)
            current_Rmin = self.windsoln.Rmin
        print(f"\r   Successfully converged Rmin to {self.windsoln.Rmin:.6f} Rp.\n",
              end ="                                                                                       ")
        self.run_isotherm(called_in_ramp_bcs=True) #TEST
        return 0
    
    
    
    def ramp_rho_rmin(self,goal_rho):
        """ Description:
                Ramps to normalized mass density at Rmin.
            Arguments:
                The desired rho value in normalized units of RHO0 (found in src/defs.h) which is default 1e-15.
        """
        rho_scale = self.windsoln.scales_dict['rho']
        if goal_rho < 10:
            print("WARNING: Rho should be in units of RHO0 =  %.0e ~nBAR*MH/K/T0"%rho_scale)
            return
        failed = 0
        OOM_old = np.copy(np.floor(np.log10(self.windsoln.rho_rmin)))
        while (abs(1.-self.windsoln.rho_rmin/goal_rho) > 1e-10):
            OOM_current = np.floor(np.log10(self.windsoln.rho_rmin))
            if OOM_old != OOM_current:
                print(f'\r..Order of mag of rho_rmin has changed, updating convergence scale.',
                      end="                                                                      ")            
                new_RHOSCALE = 10**np.floor(np.log10(self.windsoln.rho_rmin*0.001))
                h = (open('src/defs.h','r')).readlines()
                f = open('src/defs.h','w')
                for idx,hline in enumerate(h):
                    splitline = hline.split()
                    if len(splitline) >= 2:
                        line_var = splitline[0]+' '+splitline[1]
                    if line_var == '#define RHOSCALE':
                        index3 = idx
                h[index3] = '#define RHOSCALE %.1f\n' %new_RHOSCALE
                f.writelines(h)
                f.close() 
                sub = Popen('make', stdout=PIPE, stderr=PIPE) 
                output, error_output = sub.communicate() #FIX (put output check)
#                 print(error_output)
                OOM_old = np.copy(np.floor(np.log10(self.windsoln.rho_rmin)))

            if failed > 2:
                bcs_tuple = self.windsoln.bcs_tuple
                bcs_tuple[2] = self.windsoln.rho_rmin+(bcs_tuple[2]-self.windsoln.rho_rmin)/(10**failed) #here
                self.inputs.write_bcs(*bcs_tuple)
                print(f'\r..Proceeding with smaller stepsize: rho_rmin {self.windsoln.rho_rmin:.5g}, goal {goal_rho:.5g}',
                      end="                                                   ")
                fail2 = 0
                while self.run_wind(expedite=True) == 1: 
                    fail2 += 1
                    delta_rho = (bcs_tuple[2]-self.windsoln.rho_rmin)/10**(failed+fail2)
                    bcs_tuple[2] = self.windsoln.rho_rmin+delta_rho
                    self.inputs.write_bcs(*bcs_tuple)
                    print(f'\r..Attempt {failed:d}: rho_rmin {self.windsoln.rho_rmin*rho_scale:.3e} g/cm3, '
                          f'try {bcs_tuple[2]*rho_scale:.3e}',
                          end="                                               ")
                    if fail2 > 15:
                        print("\n  Struggling to substep towards new rho_rmin")
                        return 1
            else:
                bcs_tuple = self.windsoln.bcs_tuple
                bcs_tuple[2] = goal_rho
                self.inputs.write_bcs(*bcs_tuple)
                failed = 0
                print(f'\r..Goal rho_rmin: {goal_rho*rho_scale:.3e} g/cm3. Current: {self.windsoln.rho_rmin*rho_scale:.3e}.',
                      end="                                                   ")

                while self.run_wind(expedite=True) != 0:
                    failed += 1
                    if failed == 1:
                        self.run_isotherm(called_in_ramp_bcs=True) #TEST
                    delta_rho = (bcs_tuple[2]-self.windsoln.rho_rmin)/10
                    bcs_tuple[2] = self.windsoln.rho_rmin+delta_rho
                    self.inputs.write_bcs(*bcs_tuple)
                    print(f'\r..Attempt {failed:d}: rho_rmin {self.windsoln.rho_rmin*rho_scale:.3e} g/cm3, '
                          f'try {bcs_tuple[2]*rho_scale:.3e}',
                          end="                                               ")
                    if failed > 15:
                        print("\nStruggling to substep towards new rho_rmin")
                        return 1
        print(f"\r   Successfully converged rho_rmin to {self.windsoln.rho_rmin*rho_scale:.5e} g/cm3",
              end="                                                       ")
        self.run_isotherm(called_in_ramp_bcs=True) #TEST
        return 0
    
    
    
#Spectrum tools    
    def flux_norm(self,goal_flux_in_range,eV_range=[13.6,100],ramp=False,plot=True,
                 integrate_out=True,converge_bcs=True):
        '''Description: 
                Computes the total flux (across the loaded spectrum range) that is
                neccessary to achieve the goal_flux_in_range in the eV_range identified.
                If ramp = True, also ramps the solution to that Ftot.
                
           Arguments:
               goal_flux_in_range - float; desired flux in ergs/s/cm2 in range given by eV_range
               eV_range - arr or list; range to normalize to in electron volts (13.6-100 is EUV).
               
           Returns:
               0 - if ramp successful or already done.
               1,4,etc. - error codes if ramp unsuccesful.
               goal_total_flux - if ramp = False (does not ramp).
        '''
        euv = np.array(eV_range)*const.eV
        E = self.windsoln.E_wl
        flux_per_bin = self.windsoln.Ftot*self.windsoln.wPhi_wl*E
        current_euv_flux = sum(flux_per_bin[(E>euv[0])&(E<euv[1])])
        ratio = goal_flux_in_range/current_euv_flux
        goal_total_flux = ratio*self.windsoln.Ftot
        if ramp == True:
            if abs(goal_total_flux-self.windsoln.Ftot)/self.windsoln.Ftot > 0.01:
#                 print(f"Ramping flux from {self.windsoln.Ftot:.0f} to {goal_total_flux:.0f} ergs/s/cm2")
                return self.ramp_var("Ftot",goal_total_flux,make_plot=plot,
                                    converge_bcs=converge_bcs,integrate_out=integrate_out)
            else:
                print("Flux ramping already done.")
                return 0
        else:
            return goal_total_flux

        
        
    def ramp_spectrum(self,Fnorm=0.0,norm_spec_range=[],
                      goal_spec_range=[],units='eV',normalize=True,
                      kind='full',plot=False):
        #norm range and flux in norm range
        '''Description: Ramping stellar spectrum wavelength/energy range to 
                        new wavelength/energy range.
        Arguments:
            Fnorm - float; default=0.0 ergs/s/cm2 AT SEMIMAJOR AXIS OF PLANET.
                              When Fnorm=0.0, flux is normalized to the current value 
                              in norm_spec_range. Else, ramps to given Fnorm.
                              If fails, can ramp independently using self.ramp_var('Ftot')
            norm_spec_range  - list or array; default=[] is desired range over
                                which to normalize in units of 'units'.
            goal_spec_range - list or array; default=[] is current range.
                                        Else, custom upper and lower limits
                                        of spectrum in units of 'units'.
            units - str; default='eV'. Options = 'cm', 'nm'. Units of range values
            normalize - bool; default=True. If true, flux in new range
                              will be normalized to Fnorm in norm_spec_range.
            kind - str; 'full' or 'mono' - spectrum frequency type
            plot - bool; default = False.

        Returns:
            0  - ramping successful
            1  - ramping with smaller stepsize unsuccesful
            2  - need to provide Ftot or set norm_flux = True
        '''
        wl_norm = 1e-7
        if (kind == 'mono') and (self.windsoln.spec_kind == 'multi'):
            print("WARNING: Trouble ramping from multifrequency to monofrequency solutions. Start from a monofrequency solution and ramp flux using ramp_var('Ftot',...).")
            return
        
        if self.windsoln.nspecies > 4:
            print("WARNING: Ramping spectrum with large number of metals can be expensive and fail.\n  Suggestion: Ramp spectrum for H,He version of planet, then add metals.")
               
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
                print('WARNING: Spectrum range limits should be in ascending order')
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
                    print('When Fnorm=0.0, flux is normalized to the current value in norm_spec_range.')  
                    print('However, bounds of normalization range exceed the current',
                          f'bounds of spectrum. (Norm: [{norm_span[0]:.2f}{norm_span[1]:.2f}]nm',
                          f'Current: [{curr_span[0]:.2f}{curr_span[1]:.2f}]nm.)')
                    print('So, flux will be normalized to total flux in current spectrum range.')
                    norm_span[0] = max(norm_span[0],curr_span[0])
                    norm_span[1] = min(norm_span[1],curr_span[1])
                #set Fnorm to the current value in the normalization range
                E = self.windsoln.E_wl
                flux_per_bin = E*self.windsoln.Ftot*self.windsoln.wPhi_wl
                Fnorm = sum(flux_per_bin[(E<E_top_norm) & (E>E_bot_norm)])
                
                print(f'Spectrum will be normalized such that',
                      f'sum(Flux[{norm_span[0]:.2f},',
                      f'{norm_span[1]:.2f}]nm) = {Fnorm:.0f} ergs/s/cm2.')
        

        if kind != self.windsoln.spec_kind:
            print("Warning: Ramper sometimes has difficulty changing from 'full' to 'mono'.")
            print("Consider ramping from existing monofrequency solution.")
        
        if self.windsoln.spec_src_file == 'scaled-solar':
            spec = spectrum(date=self.windsoln.spec_date)
        else:
            spec = spectrum(lisird=False,spectrum_file=self.windsoln.spec_src_file)
        for sps in self.windsoln.species_list:
            spec.add_species(sps)

        print(f'Goal: {goal_span} nm')
        spec.set_resolved(*goal_span) #
        spec.set_normalized(*goal_span) 
        spec.set_window(*goal_span,kind=kind)
        spec.generate(kind=kind, savefile='inputs/spectrum.inp')

        if self.run_wind() == 0:
            if normalize == True:
                print(f'\rRamped spectrum wavelength range, now normalizing spectrum. \n ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm',
                     end='                                                     ') 
                ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                          const.hc/(norm_span[0]*wl_norm)/const.eV]
                return self.flux_norm(Fnorm,ranges,
                                   ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
            else:
                return 0

        #If it cannot make it in one leap, ramp spectrum range
        delta = np.copy(goal_span - self.windsoln.spec_resolved/self.spectrum.wl_norm)
        gen_fail = 0 
        last_spec = self.windsoln.spec_resolved
        while abs(sum((goal_span-self.windsoln.spec_resolved/self.spectrum.wl_norm)/goal_span)) > 1e-3:
            avg_diff = abs(sum((goal_span-self.windsoln.spec_resolved/self.spectrum.wl_norm)/goal_span))
            curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
            step_span = self.windsoln.spec_resolved/self.spectrum.wl_norm + delta / 2
            print(f'\rFail 1: Current [{curr[0]:.1f},{curr[1]:.1f}]. Attempting [{step_span[0]:.1f},{step_span[1]:.1f}]',
                  end='                                          ')
#             spec = spectrum()
#             for sps in self.windsoln.species_list:
#                 spec.add_species(sps)
            spec.set_resolved(*step_span) #
            spec.set_normalized(*step_span) 
            spec.set_window(*step_span,kind='full')
            spec.generate(kind='full', savefile='inputs/spectrum.inp')
            
            if np.array_equal(self.windsoln.spec_resolved,last_spec):
                if gen_fail > 5:
                    print("Spectrum ramping cannot proceed. Try ramping spectrum for an H,He atmosphere then add metals.")
                    return 1
                gen_fail += 1
                delta /= 3 
            last_spec = np.copy(self.windsoln.spec_resolved)

            fail = 1
            while self.run_wind(expedite=True) != 0:
                fail += 1
                factor = 1/(2**fail)
                if fail > 6:
                    curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                    eV = const.hc/self.windsoln.spec_resolved / const.eV
                    print(f'Failure: Too many attempts. Final range: [{curr[0]:.2f},{curr[1]:.2f}]nm ([{eV[0]:.2f},{eV[0]:.2f}]eV).')
                    print('          To ramp manually, see tutorial.')
                    if normalize == True:
                        print(f'\rNow normalizing spectrum. \n ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm',
                             end='                                                     ') 
                        ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                                  const.hc/(norm_span[0]*wl_norm)/const.eV]
                        result =  self.flux_norm(Fnorm,ranges,
                                               ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
                        curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                        print(f'Final Ftot at planet across [{curr[0]:.2e},{curr[1]:.2e}]nm = {self.windsoln.Ftot:.0f} ergs/s/cm2.')
                    return 1
                curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                step_span = curr + delta*factor
                print(f'\rFail {fail}: Current [{curr[0]:.1f},{curr[1]:.1f}]. Attempting {step_span}nm',
                      end='                                              ')
                spec.set_resolved(*step_span) #
                spec.set_normalized(*step_span) 
                spec.set_window(*step_span,kind='full')
                spec.generate(kind='full', savefile='inputs/spectrum.inp') 
#             print(f'Success! Current: [{self.windsoln.spec_resolved/ self.spectrum.wl_norm}]nm')
        else:
            if normalize == True:
                print(f'\rRamped spectrum wavelength range, now normalizing spectrum. \n ..Fnorm = {Fnorm:.0f} ergs/s/cm2. Norm range = [{norm_span[0]:.2f},{norm_span[1]:.2f}]nm',
                     end='                                                     ') 
                ranges = [const.hc/(norm_span[1]*wl_norm)/const.eV,
                          const.hc/(norm_span[0]*wl_norm)/const.eV]
                result =  self.flux_norm(Fnorm,ranges,
                                       ramp=True,plot=plot,integrate_out=True,converge_bcs=False)
                curr = self.windsoln.spec_resolved/self.spectrum.wl_norm
                print(f'Final Ftot at planet across [{curr[0]:.2e},{curr[1]:.2e}]nm = {self.windsoln.Ftot:.0f} ergs/s/cm2.')
                return result
            else:
                return 0
        
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
            lines = [v1, s0, v2, s1, v3]
            for i, span in enumerate(compare):
                old = span[0]
                new = span[1]
                for o_s, n_s in zip(old, new):
                    if abs(o_s-n_s)/o_s > 1e-2:
                        ax2.arrow(o_s, 0.75+0.1*i, n_s-o_s, 0, alpha=0.3,
                                  color=sc[i], width=0.015, head_width=0.05,
                                  head_length=3, length_includes_head=True)
                        if i == 2 and not oldplot:
                            oldplot = True
                            old_mk = ((self.spectrum.data_norm['wl']>=ons[0])
                                      &(self.spectrum.data_norm['wl']<=ons[1]))
                            s3, = ax.plot(self.spectrum.data_norm['wl'][old_mk],
                                          old_phi_smth[old_mk][old_mk], lw=1,
                                          label='Old smoothed')
                            lines += [s3]
            ax.set_yscale('log')
            ax.set_title('Spectrum ramp')
            ax.set_xlabel('Wavelength ($\lambda$) [nm]')
            ax.set_ylabel('Normalized spectral photon irradiance\n'
                          r'$\left(\phi_{\lambda}\right)$ [s$^{-1}$ cm$^{-2}$]')
            fig.tight_layout(pad=0.3)
            fig.subplots_adjust(bottom=0.3, top=0.9)
            labels = [l.get_label() for l in lines]
            fig.legend(lines, labels, bbox_to_anchor=(0.5, 0.0), loc='lower center',
                       ncol=3)
            plt.show()
            

        return 0
        
        
        
    def ramp_to_user_spectrum(self,spectrum_filename,species_list=[],
                              updated_F=0.0,norm_spec_range=[],goal_spec_range=[],
                              units='eV',normalize=True,plot=True,ramp_range=True):
        
        '''Description: Ramping stellar spectrum wavelength/energy range to 
                        new wavelength/energy range.
        Arguments:
            spectrum_filename - str; Name of the formatted spectrum file saved 
                                    in McAstro/stars/additional_spectra/
            species_list - list of str; default = [] sets to species list of 
                                    the loaded windsoln.
            updated_F - float; default=0.0 ergs/s/cm2 AT SEMIMAJOR AXIS OF PLANET.
                              When =0.0, flux is normalized to the user spectrum's og value
                              in norm_spec_range. Else, ramps to given updated_F.
                              If fails, can ramp independently using self.ramp_var('Ftot')
            norm_spec_range  - list or array; default=[] is desired range over
                                which to normalize in units of 'units'.
            goal_spec_range - list or array; default=[] is current range.
                                        Else, custom upper and lower limits
                                        of spectrum in units of 'units'.
            units - str; default='eV'. Options = 'cm', 'nm'. Units of range values
            normalize - bool; default=True. If true, flux in new range
                              will be normalized to Fnorm in norm_spec_range.
            kind - str; 'full' or 'mono' - spectrum frequency type
            plot - bool; default = False.
            '''
        #generating a smoothed and binned version of the user-input code
        wl_norm=1e-7
        spec = spectrum(lisird=False,spectrum_file=spectrum_filename,
                        wl_norm=wl_norm,just_loading=False)

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
        spec.generate(kind='full', savefile='inputs/goal_spectrum.inp')
        spec.generate(kind='full', savefile='inputs/spectrum.inp')
        
        old_E = self.windsoln.sim_spectrum['E']
        old_wPhi = self.windsoln.sim_spectrum['wPhi']
        
        self.inputs.write_flags(*self.windsoln.flags_tuple,
                                integrate_out=False)
        if self.run_wind() != 0:
            new = np.genfromtxt('inputs/goal_spectrum.inp',
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
                    spec.generate(kind='full', savefile='inputs/spectrum.inp')
                    if self.run_wind() != 0:
                        fail = 0
                        percent = 0.1

                        current_wPhi = np.copy(self.windsoln.sim_spectrum['wPhi'])
                        step_wPhi = np.array(current_wPhi+delta_wPhi*(percent))
                        #if next step is larger than the remainder, set step to goal
                        if avg < abs(np.average(delta_wPhi*percent)):
                            step_wPhi = wPhi
                        if plot == False:
                            print(f"\r..Trying a {percent*100:.0f}% step",end="                      ")
                        step_wPhi[0] = step_wPhi[1]
                        step_wPhi[-1] = step_wPhi[-2]
                        self.inputs.write_spectrum(*self.windsoln.spectrum_tuple[:3],
                                                  spectrum_filename,
                                                  *self.windsoln.spectrum_tuple[4:10],
                                                  step_wPhi,
                                                  *self.windsoln.spectrum_tuple[11:]) 
                        self.inputs.write_flags(*self.windsoln.flags_tuple,
                                                integrate_out=False)
                        while self.run_wind() != 0:
                            fail+=1
                            print(f"\rFail {fail:d}: delta = {percent/(fail+1):.3f}",
                                 end ="                                                   ")
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
                                print("Failed to ramp to new stellar spectrum.")
                                return 1

                        self.converge_Ncol_sp(expedite=True,quiet=True)
                        self.run_isotherm()
                        self.inputs.write_flags(*self.windsoln.flags_tuple,integrate_out=False)

                        avg = abs(np.average((wPhi - self.windsoln.sim_spectrum['wPhi'])))
                        if plot == True:
                            pl.clf()
                            pl.plot(const.hc/old_E/1e-7,start_wPhi(old_E),label='Original')
                            pl.plot(const.hc/old_E/1e-7,wPhi,label='Goal')
                            pl.plot(const.hc/old_E/1e-7,self.windsoln.sim_spectrum['wPhi'],
                                         label='Current',ls='--')
                            pl.xlabel('Wavelength (nm)')
                            pl.ylabel('Photon Density')
                            pl.legend()
                            pl.yscale('log')
                            pl.title(f"\r Success! Average difference now {avg:.2e}")
                            display.display(pl.gcf())
                            display.clear_output(wait=True)
                        else:
                            print(f"\r Success! Average difference now {avg:.2e}",
                                 end='                                                  ')
            except ValueError:
                print(f"Succeeded in full jump to {spectrum_filename}")

        #Ramping to desired spectral range and flux
        self.converge_Ncol_sp(expedite=True,quiet=True)
        self.run_isotherm()
        print(f'\nSuccess ramping to {spectrum_filename} spectrum shape!')
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
            print('     Now ramping flux and spectral range.')
            if len(goal_spec_range) == 0:
                goal_spec_range = self.windsoln.spec_resolved/self.spectrum.wl_norm
                units = 'nm'
            if updated_F == 0:
                f = open("McAstro/stars/spectrum/additional_spectra/"+spectrum_filename,"r")
                updated_F = float(f.readlines()[1])/(self.windsoln.semimajor/const.au)**2
                norm_spec_range = [12.4,91.1]
                f.close()
            if self.ramp_spectrum(updated_F,norm_spec_range,
                                  goal_spec_range,units,normalize,
                                  kind='full',plot=plot) == 0:
                print("Success! Ramped to user-input stellar spectrum "+spectrum_filename)
                return 0
            else:
                print("Ramped successfully to user-input spectrum, but ramping spectral range and/or flux was not successful.")
                return 2
        else:
            return 0

        
        
    def format_user_spectrum(self,wl,flux_1au,wl_units='cm',spectrum_name='',
                            comment='',overwrite=False):
        '''
        Description: Users may import a custom spectrum. This function
                    formats the spectrum to be readable by the code.

        Arguments:
        wl - list/array; wavelengths.
        wl_units - str; default='cm'. Options = 'nm','m','A' (angstroms)
        flux_1au - list/array; flux per wavelength bin AT 1 AU FROM STAR
                    in ergs/s/cm2. 
                    Stellar spectra are often given from Earth in pc. 
                    To scale, multiply by (dist to star in au / 1 au)^2.
                    NOT FLUX DENSITY (to convert from flux density, multiply
                    by delta(wl) - be careful of units).
        spec_name - str; name of file. Will be saved in solution spectrum.
        comment  - str; any comments (name of spectrum and other info). 
                    No need to add newline character.
        overwrite - bool; if True overwrite existing file of same name.
        '''
        spectrum_name = spectrum_name.replace(' ','_')
        wl = np.array(wl)
        flux_1au = np.array(flux_1au)

        if len(wl) != len(flux_1au):
            print("Wavelength and flux at 1 au arrays must be same length.")
            return    

        file = 'McAstro/stars/spectrum/additional_spectra/'+spectrum_name

        if (exists(file)==True) and (overwrite == False):
            print("To overwrite existing file, set overwrite=True.")
            return

        g = open(file,'w')
        if len(comment) > 0:
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
            print("Invalid units: Options = 'A','nm','cm',or 'm'.")
            return
        if np.median(flux_1au) < 1:
            print("WARNING: Flux at 1 au from star should be in units of ergs/s/cm2 (check that your input is not the flux at Earth).")
        
        total_EUV_flux = sum(flux_1au[(wl>min(np.max(wl),1.2e-6))&(wl<max(np.min(wl),9.1e-6))])
        g.write('%.5f\n'%total_EUV_flux)
        
        #Converting to ergs/s/cm2/cm flux density
        delta_lam = np.diff(wl_cm,prepend=wl_cm[0]-(wl_cm[1]-wl_cm[0]))
        flux_dens_1au = flux_1au*delta_lam
        
        df = pd.DataFrame(np.column_stack((
            wl_cm,
            flux_dens_1au,
            np.zeros_like(wl_cm),
            const.c/wl_cm,
            flux_dens_1au*wl_cm**2/const.c,
        )), columns=['wl','F_wl','unc','nu','F_nu'])

        df.to_csv(g, header=True, index=False)
        g.close()

        print('Spectrum saved in Wind-AE readable format at',file)

        return