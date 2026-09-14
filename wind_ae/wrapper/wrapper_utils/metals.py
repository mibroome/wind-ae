#!/usr/bin/env python3
"""
metals.py:
Contains the class for ramping in new species to existing solutions,
as well as computing mass fractions of species for a given metallicity. 
"""

import sys 
import os
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import importlib.resources as pkg_resources

from wind_ae.wrapper.wrapper_utils import constants as const
from wind_ae.wrapper.wrapper_utils.spectrum import spectrum
from wind_ae.wrapper.wrapper_utils.inputs import input_handler
import wind_ae.McAstro.atoms.atomic_species as McAtom
from scipy.interpolate import CubicSpline
from scipy.special import exp1

class metal_class:
    def __init__(self, windsoln_object, workdir=None):
        self.path = str(pkg_resources.files('wind_ae'))+'/'
        # wpath: writable path for inputs/ (spectrum.inp, guess.inp, etc.)
        # Inherits the wind_simulation's workdir so parallel instances stay isolated.
        if workdir is None:
            self.workdir = None
            self.wpath = self.path
        else:
            self.workdir = str(workdir).rstrip('/')
            self.wpath = self.workdir + '/'
        self.windsoln = windsoln_object
        self.df = pd.read_csv(pkg_resources.files('wind_ae.wrapper.wrapper_utils').joinpath('dere_table29.dat'),
                              sep=r'\s+', names=list(range(45)))  
        self.df = self.df.rename(columns={0:'Z', 1:'Ion', 2:'NS', 3:'I', 4:'Tmin'})
        # filepath = pkg_resources.files('wind_ae.McAstro.atoms').joinpath('Lodders.dat')
        # self._lodders = pd.read_csv(filepath,comment='#',delimiter=' ')
        
        # filepath = pkg_resources.files('wind_ae.McAstro.atoms').joinpath('atomic_masses.csv')
        # self._amu = pd.read_csv(filepath, comment='#')
        self._Zgrid = pd.read_csv(pkg_resources.files("wind_ae").joinpath("wrapper/wrapper_utils/metallicity_grid.csv"))                                    

    def load_spectrum_add_metals(self, generate=False,wl_norm=1e-7):
        """Identical behavior to load_spectrum() function in relax_wrapper.py.
        Loads the high resolution spectrum referenced in the guess solution so that it can be re-smoothed when new metals are added to the solution. 
        
        Args: 
            generate (bool): if True, generates a new smoothed spectrum stored in inputs/spectrum.inp
            wl_norm (float): converts nm to cm. Should not need to be changed.
        """
        if self.windsoln.spec_src_file != 'scaled-solar':
            self.spectrum = spectrum(lisird=False,
                                     spectrum_file=self.windsoln.spec_src_file,
                                     wl_norm=wl_norm)
        else:
            self.spectrum = spectrum(date=self.windsoln.spec_date, wl_norm=wl_norm)
        for j in range(self.windsoln.nspecies):
            self.windsoln.species_list[j] = (self.windsoln.species_list[j]).replace(' ','')
        # A monofrequency windsoln has a degenerate (single-wavelength)
        # spec_resolved/spec_normalized/spec_window, which the glq_spectrum
        # smoothing/binning pipeline underneath add_species()/set_resolved()/
        # set_normalized()/set_window()/generate() cannot handle (see
        # spectrum.write_mono_spectrum() for details). Skip all of that and
        # write the (possibly updated) species' cross sections directly at the
        # existing wavelength instead.
        if self.windsoln.spec_kind.lower() == 'mono':
            if generate:
                mono_nm = self.windsoln.spec_resolved[0]/wl_norm
                self.spectrum.write_mono_spectrum(mono_nm, self.windsoln.species_list,
                                                  savefile=self.wpath+'inputs/spectrum.inp')
            return
        for name in self.windsoln.species_list:
            self.spectrum.add_species(name)
        soln_resolved = self.windsoln.spec_resolved/wl_norm
        self.spectrum.set_resolved(*soln_resolved)
        soln_normalized = self.windsoln.spec_normalized/wl_norm
        self.spectrum.set_normalized(*soln_normalized)
        soln_window = self.windsoln.spec_window/wl_norm
        self.spectrum.set_window(*soln_window, kind=self.windsoln.spec_kind)
        if generate:
            self.spectrum.generate(kind=self.windsoln.spec_kind, savefile=self.wpath+'inputs/spectrum.inp')
        return
    
    def _generate_rate_coeffs(self):
        """ Identical to the _generate_rate_coeffs() function in relax_wrapper.py.
        
        Generates rate coefficients for secondary ionization and populates 
        array in inputs/rate_coeffs.inp. Rate coefficients are generated from interpolation 
        over the Dere (2007) table.
        """
        new_spec = np.genfromtxt(self.wpath+'inputs/spectrum.inp',skip_header=9,delimiter=',')
        E_wl = new_spec[:,0]
        f = open(self.wpath+'inputs/spectrum.inp','r')
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

        # Write rate coefficients to inputs/rate_coeffs.inp (read at runtime by
        # load_rate_coeffs() in glq_rates.c) instead of src/rate_coeffs.h.
        # Format: comment header then nspecies^2 rows of npts comma-separated
        # floats; row index j*nspecies+m is the rate for a photoelectron
        # released from species j to secondarily ionize species m.
        # (The old rate_coeffs.h declared nspecies^2+nspecies rows but only
        # filled nspecies^2; the extra were zero-padding, never accessed.)
        nspecies = len(species_new)
        f = open(self.wpath+'inputs/rate_coeffs.inp', 'w')
        f.write('# Rate coefficients (cm^3/s) from Dere (2007)\n')
        f.write('# Row j*nspecies+m: photoelectron from species j ionizes species m\n')
        f.write('# NSPECIES: %d, NPTS: %d, NROWS: %d\n' % (nspecies, len(E_wl), nspecies**2))
        for j in range(nspecies):
            for m in range(nspecies):
                Rs = Interpolater(j, m)
                f.write(','.join('%.4e' % rs for rs in Rs) + '\n')
        f.close()
    

    def add_species_to_guess(self, new_species_names):
        """
        Adds new species to the list of species in the wind solution object and updates
        all species-dependent properties.

        Args:
            new_species_names (list of str): The names of the new species to add. (e.g., ['O3','c1','Ne II'])
        """
        ns_unspaced = []
        ns_spaced = []
        for new_species_name in new_species_names:
            try:
                new_species_spaced = McAtom.formatting_species_list([new_species_name])[0] 
                #adds spaces and converts, e.g., fe6 to Fe IX
                new_species_unspaced = new_species_spaced.replace(' ','')
                ns_unspaced = np.append(ns_unspaced,new_species_unspaced)
                ns_spaced = np.append(ns_spaced,new_species_spaced)
            except IndexError:
                print("ERROR: Species must include an ionization number either in Roman or Arabic numerals (e.g., H1, heIV, Fe 6)", 
                      file=sys.stderr) # Python 3.x
                sys.exit(1)

            if self.windsoln.spec_src_file != 'scaled-solar':
                spec = spectrum(lisird=False,spectrum_file=self.windsoln.spec_src_file)#,
            else:
                spec = spectrum(date=self.windsoln.spec_date)
            if new_species_unspaced in self.windsoln.species_list:
                print(new_species_unspaced+' already in simulation. Current species list:')
                print('   ',self.windsoln.species_list)
                return 5
            self.windsoln.species_list.append(new_species_unspaced)
        species_list_unspaced = np.copy(self.windsoln.species_list)
    #     McAtom.formatting_species_list(self.windsoln.species_list) #reintroduces space
        #generate new spectrum
        self.load_spectrum_add_metals(generate=True)
        new_spec = np.genfromtxt(self.wpath+'inputs/spectrum.inp',skip_header=8,delimiter=',')
        #rewriting the newly generated spectrum columns to add new species column
        try:
            self.windsoln.E_wl     = new_spec[:,0]
            self.windsoln.wPhi_wl  = new_spec[:,1]
            self.windsoln.npts     = len(self.windsoln.wPhi_wl)
            self.windsoln.sigma_wl = new_spec[:,2:]
        except IndexError:
            self.windsoln.E_wl     = new_spec[0]
            self.windsoln.wPhi_wl  = new_spec[1]
            self.windsoln.npts     = 1
            self.windsoln.sigma_wl = new_spec[2:]


        #rewriting Ys and Ncol columns to add new species column
        #takes last Ys column and duplicates it for the new species
        last_Ncol = np.copy(self.windsoln.soln_norm.iloc[:,3+2*len(self.windsoln.Ys_rmin)])
        for new_species_unspaced,new_species_spaced in zip(ns_unspaced,ns_spaced):
            prior_nspecies = len(self.windsoln.Ys_rmin)
            self.windsoln.soln_norm.insert(4+prior_nspecies, 'Ys_'+new_species_unspaced, 
                                          self.windsoln.soln_norm.iloc[:,3+prior_nspecies], 
                                          allow_duplicates=True) 
            #takes last Ncol column (snapshotted once before this loop, so every
            #new species in a multi-species add is seeded from the same
            #pre-existing column rather than chaining off each other's
            #already-tiny 1e-10 seed) and *1e-10 it for the new species
            self.windsoln.soln_norm.insert(5+2*prior_nspecies, 'Ncol_'+new_species_unspaced, 
                                        last_Ncol*1e-10,
                                        allow_duplicates=True) 
            #Adding other species details
            ma = McAtom.atomic_species(new_species_spaced)

            # Detect whether this new species is a chain child of an existing species.
            # A chain child shares the same element Z and has one more electron stripped
            # (Ne_new == Ne_existing - 1).  If so, its HX must equal the parent's HX
            # (both refer to the full element mass fraction); otherwise start near 0.
            new_Z, new_Ne = McAtom.spectroscopy_to_atomic_notation(new_species_spaced)
            chain_hx = 1e-10  # default: tiny seed for a truly new element
            for _p_idx, _p_sp in enumerate(self.windsoln.species_list[:-1]):  # excludes just-appended entry
                _p_sp_spaced = McAtom.formatting_species_list([_p_sp])[0]
                _p_Z, _p_Ne = McAtom.spectroscopy_to_atomic_notation(_p_sp_spaced)
                if new_Z == _p_Z and new_Ne == _p_Ne - 1:
                    chain_hx = self.windsoln.HX[_p_idx]
                    break

            self.windsoln.HX = np.append(self.windsoln.HX, chain_hx)
            self.windsoln.ion_pot = np.append(self.windsoln.ion_pot,ma.verner_data['E_th']*const.eV)
            self.windsoln.atomic_masses = np.append(self.windsoln.atomic_masses,ma.mass) 

            self.windsoln.Ys_rmin = np.append(self.windsoln.Ys_rmin,1.) #Neutral fraction
            self.windsoln.Ncol_sp = np.append(self.windsoln.Ncol_sp,1e-10) #0 mass fraction in the new species for now
            self.windsoln.scales.insert(4+prior_nspecies,1.)
            self.windsoln.scales.insert(5+2*prior_nspecies,1e17)

            self.windsoln.nspecies = len(self.windsoln.Ys_rmin) #reassigning new nspecies

            self.windsoln.physics_tuple = [self.windsoln.HX, 
                                           self.windsoln.species_list,
                                           self.windsoln.molec_adjust,
                                           self.windsoln.atomic_masses,
                                           self.windsoln.kappa_opt,
                                           self.windsoln.kappa_IR,
                                           self.windsoln.gamma] 
            self.windsoln.bcs_tuple = [self.windsoln.Rmin, 
                                       self.windsoln.Rmax, 
                                       self.windsoln.rho_rmin, 
                                       self.windsoln.T_rmin,
                                       self.windsoln.Ys_rmin, 
                                       self.windsoln.Ncol_sp,
                                       self.windsoln.erf_drop]
            self.windsoln.spectrum_tuple = [self.windsoln.npts, 
                                            self.windsoln.nspecies,
                                            self.windsoln.spec_date,
                                            self.windsoln.spec_src_file,
                                            self.windsoln.spec_kind,
                                            self.windsoln.spec_window,
                                            self.windsoln.spec_resolved,
                                            self.windsoln.spec_normalized,
                                            self.windsoln.ion_pot, 
                                            self.windsoln.E_wl, 
                                            self.windsoln.wPhi_wl,
                                            self.windsoln.sigma_wl,
                                            self.windsoln.species_list]

        #Rewriting guess with new species, next step is ramping up the mass fraction
        inputs = input_handler(workdir=self.workdir)
        inputs.write_physics_params(*self.windsoln.physics_tuple)
        inputs.write_bcs(*self.windsoln.bcs_tuple)
        inputs.write_spectrum(*self.windsoln.spectrum_tuple)
        self.rewrite_guess()

        new_length = len(self.windsoln.soln['q'][self.windsoln.soln['q']<=1])

        #writing rate coefficient file
        self._generate_rate_coeffs()

        # M is #define M g_m — a runtime global read from tech_params.inp at startup.
        # No defs.h rewrite or recompile needed; write_solver_params is sufficient.
        new_RHOSCALE = 10**np.floor(np.log10(self.windsoln.rho_rmin * 0.01))
        _it = int(1e3) if self.windsoln.conduction == 1 else 100
        input_handler(workdir=self.workdir).write_solver_params(new_length, _it, new_RHOSCALE)
        return
    
    
    
    def remove_species_from_guess(self,remove_species):
        """
        Removes specified species from the wind solution object. Likely to break if removing Helium,
        so it is recommended to ramp directly from pure-H solutions provided.

        Args:
            remove_species (list of str): The names of the species to remove. (e.g., ['O3','c1','Ne II'])
        """
        remove_species_spaced = McAtom.formatting_species_list(remove_species)
        remove_species_unspaced = [rs.replace(' ','') for rs in remove_species_spaced]

        self.windsoln.species_list = [sp.replace(' ','') for sp in self.windsoln.species_list]  
        # if remove_species_unspaced in self.windsoln.species_list:
        for j in range(len(remove_species)):
            rsu = remove_species_unspaced[j]
            idx = self.windsoln.species_list.index(rsu)
            #physics
            self.windsoln.scales  = np.delete(self.windsoln.scales,
                                             [4+idx,(4+self.windsoln.nspecies)+idx])
            self.windsoln.ion_pot = np.delete(self.windsoln.ion_pot,idx)
            self.windsoln.atomic_masses = np.delete(self.windsoln.atomic_masses,idx)
            #bcs
            self.windsoln.Ys_rmin = np.delete(self.windsoln.Ys_rmin, idx)
            self.windsoln.Ncol_sp = np.delete(self.windsoln.Ncol_sp, idx)
            #solution
            self.windsoln.soln_norm = self.windsoln.soln_norm.drop(columns = ['Ys_'+rsu, 'Ncol_'+rsu])
            #species_list
            self.windsoln.species_list.remove(rsu)
            self.windsoln.nspecies = len(self.windsoln.species_list)
            
        self.load_spectrum_add_metals(generate=True)
        new_spec = np.genfromtxt(self.wpath+'inputs/spectrum.inp',skip_header=8,delimiter=',')
        #rewriting the newly generated spectrum columns to add new species column
        self.windsoln.E_wl     = new_spec[:,0]
        self.windsoln.wPhi_wl  = new_spec[:,1]
        self.windsoln.npts     = len(self.windsoln.wPhi_wl)
        self.windsoln.sigma_wl = new_spec[:,2:]
        
        self.windsoln.HX = self.metallicity(self.windsoln.species_list)
#         print(self.metallicity(self.windsoln.species_list))

        self.windsoln.physics_tuple = [self.windsoln.HX, 
                                       self.windsoln.species_list, 
                                       self.windsoln.molec_adjust,
                                       self.windsoln.atomic_masses,
                                       self.windsoln.kappa_opt,
                                       self.windsoln.kappa_IR,
                                       self.windsoln.gamma] 
        self.windsoln.bcs_tuple = [self.windsoln.Rmin, 
                                   self.windsoln.Rmax, 
                                   self.windsoln.rho_rmin, 
                                   self.windsoln.T_rmin,
                                   self.windsoln.Ys_rmin, 
                                   self.windsoln.Ncol_sp,
                                   self.windsoln.erf_drop]
        self.windsoln.spectrum_tuple = [self.windsoln.npts, 
                                        self.windsoln.nspecies, 
                                        self.windsoln.spec_date,
                                        self.windsoln.spec_src_file,
                                        self.windsoln.spec_kind,
                                        self.windsoln.spec_window,
                                        self.windsoln.spec_resolved,
                                        self.windsoln.spec_normalized,
                                        self.windsoln.ion_pot, 
                                        self.windsoln.E_wl, 
                                        self.windsoln.wPhi_wl,
                                        self.windsoln.sigma_wl, 
                                        self.windsoln.species_list]

        #Rewriting guess with removed species
        inputs = input_handler(workdir=self.workdir)
        inputs.write_physics_params(*self.windsoln.physics_tuple)
        inputs.write_bcs(*self.windsoln.bcs_tuple)
        inputs.write_spectrum(*self.windsoln.spectrum_tuple)
        self.rewrite_guess()

        #writing rate coefficient file (species list/spectrum changed above)
        self._generate_rate_coeffs()

        # NSPECIES is now read at runtime from phys_params.inp — no defs.h rewrite or make needed here.
        return

        
    def rewrite_guess(self,regrid=False,degree=1,ramping_metals=False):
        '''Reformats wind_ae/inputs/guess.inp (the guess needed for relaxation routine). 
            Only called when new species is added to solution. 

            Args:
                regrid (bool): If True, regrids the solution with a dynamic grid size.
                degree (int): The degree of the polynomial for regridding.
                ramping_metals (bool): If True, applies ramping to metal species.

        '''
        for i,sp in enumerate(self.windsoln.species_list):
            self.windsoln.species_list[i] = (sp).replace(' ','')
        g = open(self.wpath+"inputs/guess.inp","w")
        #Header
        self.windsoln.nspecies = len(self.windsoln.Ys_rmin)
        g.write("#nspecies: %d\n" %self.windsoln.nspecies)
        Ncol_string = ''
        Ys_string = ''
        string = ''
        floats = ''
        for n in range(self.windsoln.nspecies):
            floats += '{:.17e},'
            string += '{:s},'
            Ncol_string += 'Ncol_'+self.windsoln.species_list[n]+','
            Ys_string   += 'Ys_'+self.windsoln.species_list[n]+','
        g.write("#vars: r,rho,v,T,"+Ys_string+Ncol_string+"q,z\n")
        g.write('#scales: '+(','.join('%.17e' %i for i in self.windsoln.scales))+ '\n')
        g.write('#plnt_prms: '+(','.join('%.17e' %i for i in (self.windsoln.planet_tuple)))+ '\n')
        g.write('#phys_prms: '+(floats+string+floats).format(*self.windsoln.HX,*self.windsoln.species_list,
                                                             *self.windsoln.atomic_masses)[:-1]
                +',{:.5f},{:.2e},{:.2e},{:.5f}\n'.format(self.windsoln.molec_adjust,self.windsoln.kappa_opt,
                                                          self.windsoln.kappa_IR,self.windsoln.gamma))
        g.write('#bcs: '
                +(','.join('{:.17e}' for i in range(6+2*self.windsoln.nspecies))).format(*self.windsoln.bcs_tuple[:4], *self.windsoln.Ys_rmin, *self.windsoln.Ncol_sp,self.windsoln.erf_drop[0], self.windsoln.erf_drop[1])+'\n')
        g.write('#tech: '+(','.join('%.17e' %i for i in self.windsoln.tech_tuple))+ '\n')
        g.write('#flags: {:d},{:.5f},{:d},{:.5f},{:d},{:d},{:d},{:.5f}\n'.format(*self.windsoln.flags_tuple))
        g.write('#add_prms: %d' %self.windsoln.add_tuple[0])
        if self.windsoln.n_add_prms > 0:
            for i in range(self.windsoln.n_add_prms):
                g.write(',%s,%.17e'%(self.windsoln.add_param_name[i],self.windsoln.add_param_val[i]))
        g.write('\n')
        g.write('#spec_prms: {:d},{:d},'.format(self.windsoln.npts,self.windsoln.nspecies)+self.windsoln.spec_date
                +','+self.windsoln.spec_src_file+','+self.windsoln.spec_kind+
                ',{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'.format(*self.windsoln.spec_window,
                                                                           *self.windsoln.spec_resolved, 
                                                                           *self.windsoln.spec_normalized))
        g.write(floats[:-1].format(*self.windsoln.ion_pot)+'\n')   
        #Spectrum
        #just adds new species at this stage. In other functions, can expand or change range,norm,etc.
        g.write(r"# $hc/\lambda_i$, $w_i\Phi_{\lambda_i}/F_{tot}$,"+
                (','.join(r"$\sigma_{\lambda_i,"+sp+r"}$" for sp in self.windsoln.species_list))+"\n")
        if self.windsoln.npts != 1:
            for n in range(self.windsoln.npts):
                g.write(('## {:.17e},{:.17e},'+floats[:-1]+'\n').format(self.windsoln.E_wl[n],self.windsoln.wPhi_wl[n],
                                                                        *self.windsoln.sigma_wl[n]))
        else:
            g.write(('## {:.17e},{:.17e},'+floats[:-1]+'\n').format(self.windsoln.E_wl,self.windsoln.wPhi_wl,
                                                                        *self.windsoln.sigma_wl))
            
        if ramping_metals == True:
            factor = self.windsoln.Ncol_sp[-1]/self.windsoln.Ncol_sp[-2]
            self.windsoln.soln_norm.iloc[:,-5] = factor*self.windsoln.soln_norm.iloc[:,-6]
            self.windsoln.soln_norm.to_csv(g, header=False, index=False, float_format='%.17e')

            # for k in range(len(self.windsoln.soln_norm)):
            #     #multiply the guess Ncol column for the last added species by the mass frac to speed solving (hopefully)
            #     self.windsoln.soln_norm.iloc[k,-3] = factor*self.windsoln.soln_norm.iloc[k,-4]
            #     g.write(','.join('%.17e' %i for i in self.windsoln.soln_norm.iloc[k,:]))
            #     g.write('\n')
        #Guess Solution
        if regrid == False:
            self.windsoln.soln_norm.to_csv(g, header=False, index=False, float_format='%.17e')
            # for k in range(len(self.windsoln.soln_norm)):
            #     g.write(','.join('%.17e' %i for i in self.windsoln.soln_norm.iloc[k,:]))
            #     g.write('\n')
            # g.close()
        g.close()

        return

    def metallicity(self,species_list,Z=1):
        '''Given a metallicity (Z, where Z=1 is solar), returns an array of mass 
            fractions for the species in species_list. NOTE: Mass fraction must 
            sum to 1, so delta=1-sum(mass fraction Z) is added to the mass fraction
            of hydrogen.
            Relative abundances from Lodders (2009).

           Args:
               species_list (list of str): The list of species for which to compute metallicity
                                            (e.g., ['O3','c1','Ne II'])
               Z (int or float): The metallicity (default is 1 for solar).
                                  Non-integer values are linearly interpolated
                                  between adjacent grid rows.

            Returns:
                numpy.ndarray: An array of mass fractions for the species in species_list.
        '''

        try:
            species_list = McAtom.formatting_species_list(species_list) 
            #adds spaces and converts, e.g., fe6 to Fe IX
        except IndexError:
            print("ERROR: Species must include an ionization number either in Roman or Arabic numerals (e.g., H1, heIV, Fe 6)", 
                  file=sys.stderr) # Python 3.x
            sys.exit(1)

        el_list = [species.split()[0] for species in species_list]

        # When chain ionization states are present (e.g. CI and CII both in the
        # species list), el_list contains duplicate element names.  Look up the
        # metallicity grid using only unique elements, then map back to the full
        # list so that chain-linked slots receive the same mass fraction.
        unique_els = []
        seen = set()
        for el in el_list:
            if el not in seen:
                unique_els.append(el)
                seen.add(el)

        solar_mass_fracs = self._Zgrid[unique_els].to_numpy()
        # Support fractional Z by linearly interpolating between adjacent grid
        # rows.  Integer Z follows the original path exactly.
        Z_lo = int(np.floor(Z))
        Z_hi = int(np.ceil(Z))
        if Z_lo == Z_hi:  # integer Z — original behaviour
            Z_unique = solar_mass_fracs[Z_lo-1, :] / np.sum(solar_mass_fracs[Z_lo-1, :])
        else:
            frac = Z - Z_lo
            mf_lo = solar_mass_fracs[Z_lo-1, :] / np.sum(solar_mass_fracs[Z_lo-1, :])
            mf_hi = solar_mass_fracs[Z_hi-1, :] / np.sum(solar_mass_fracs[Z_hi-1, :])
            Z_unique = (1.0 - frac) * mf_lo + frac * mf_hi
        el_to_mf = dict(zip(unique_els, Z_unique))
        Z_array = np.array([el_to_mf[el] for el in el_list])

        # N_per_N_H = 10**(lodders['A']-12)
        # mn = np.zeros_like(N_per_N_H)
        # for i,el in enumerate(lodders['name']):
        #     mass = (amu['mass'][amu['Symbol'] == el]).iloc[0]
        #     if (el != 'H') and (el != 'He'):
        #         mn[i] = mass*Z*N_per_N_H[i]
        #     else:
        #         mn[i] = mass*N_per_N_H[i]
        # solar_mass_frac = mn/np.sum(mn)
        # Z_array = []
        # for element in el_list:
        #     Z_array = np.append(Z_array,solar_mass_frac[lodders['name']==element])
        # Z_array = Z_array/np.sum(Z_array)
        # temp_delta = (1-np.sum(Z_array))
        # Z_array[0] += temp_delta

        return Z_array

#     def rewrite_mass_frac(self,goal_mass_fracs,Ncol_sp_tot,integrate_out=False):
#         '''Description: 
#                 Rewrites the physics parameters input file (phys_params.inp)
#                 and boundary conditions input file (bcs.inp) for a new mass fraction.
#                 Column density in each species is taken to be Ncol_sp_tot*mass_fraction
#                 and converged to self consistent values at the last step.
#            Inputs:
#                mass_frac      - array_like; mass fraction in each species
#                Ncol_sp_tot    - float; total column density at sonic point (1/cm2)
#                integrate_out  - bool; can turn off integration for speed
#            Returns:
#                None
#                 '''
#         if len(goal_mass_fracs) != len(self.windsoln.species_list):
#             sys.exit("ERROR: Mass fraction and species list must be the same length.")
#         if np.round(np.sum(goal_mass_fracs),5) != 1:
#             print('WARNING: Total Mass Fraction must sum to 1. sum(ZX) = %.3f' %np.sum(goal_mass_fracs))
#         for i,sp in enumerate(self.windsoln.species_list):
#             self.windsoln.species_list[i] = (sp).replace(' ','')


#         self.windsoln.HX = np.array(goal_mass_fracs)
#         self.windsoln.Ncol_sp = Ncol_sp_tot*self.windsoln.HX

#         self.windsoln.physics_tuple = [self.windsoln.HX,
#                                        self.windsoln.species_list,
#                                        self.windsoln.molec_adjust, 
#                                       self.windsoln.atomic_masses] 
#         self.windsoln.bcs_tuple = [self.windsoln.Rmin, 
#                                    self.windsoln.Rmax, 
#                                    self.windsoln.rho_rmin,
#                                    self.windsoln.T_rmin,
#                                    self.windsoln.Ys_rmin, 
#                                    self.windsoln.Ncol_sp,
#                                    self.windsoln.erf_drop]

#         inputs = input_handler(workdir=self.workdir)
#         inputs.write_physics_params(*self.windsoln.physics_tuple)
#         inputs.write_bcs(*self.windsoln.bcs_tuple)
#         self.rewrite_guess() #does not change columns of guess, because just has to be close enough
#         gmf = np.array(goal_mass_fracs)

#         physics_tuple = [gmf, self.windsoln.species_list,
#                                        self.windsoln.molec_adjust, 
#                                       self.windsoln.atomic_masses] 
#         bcs_tuple = [self.windsoln.Rmin, self.windsoln.Rmax, 
#                                    self.windsoln.rho_rmin, 
#                                    self.windsoln.T_rmin,
#                                    self.windsoln.Ys_rmin,
#                                    Ncol_sp_tot*gmf,
#                                    self.windsoln.erf_drop]

#         inputs = input_handler(workdir=self.workdir)
#         inputs.write_physics_params(*physics_tuple)
#         inputs.write_bcs(*bcs_tuple)
#         return