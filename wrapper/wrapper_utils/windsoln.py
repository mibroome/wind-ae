#!/usr/bin/env python3

import re
import numpy as np
import pandas as pd
from scipy import integrate, interpolate, optimize
import matplotlib.pyplot as plt
import McAstro.atoms.atomic_species as McAtom
import time
import math

from . import constants as const

class wind_solution:
    """
    Class containing relaxation result, object typically referred to as the
    planet. Also contains functions for analyzing relaxation result
    """
    def __init__(self, file='saves/windsoln.csv', expedite=False,
                 add_uservars=False):
        self.error = 0
        with open(file, 'r') as f:
            i_spec = 0
            for line in f:
                line = re.findall(r'[^,\s]+', line)
                if line[0] == '#nspecies:':
                    self.nspecies = int(line[1:][0])
                elif line[0] == '#vars:':
                    self.varz = line[1:]
                elif line[0] == '#scales:':
                    line[1:] = [float(x) for x in line[1:]]
                    self.scales = line[1:]
                    self.scales_dict = {v:s for v, s in
                                        zip(self.varz, self.scales)}
                elif line[0] == "#plnt_prms:":
                    plnt_prms = [float(x) for x in line[1:]]
                    # planet parameters
                    self.Mp = plnt_prms[0]
                    self.Rp = plnt_prms[1]
                    self.Mstar = plnt_prms[2]
                    self.semimajor = plnt_prms[3]
                    self.Ftot = plnt_prms[4]
                    self.Lstar = plnt_prms[5]
                elif line[0] == "#phys_prms:":
                    phys_prms = [x for x in line[1:]]
                    # physics paramters
                    mass_fractions = np.zeros(self.nspecies)
                    atom_masses = np.zeros(self.nspecies)
                    species_names = []
                    for s in range(self.nspecies): #CHANGED
                        mass_fractions[s] = float(phys_prms[s])
                        species_names.append(phys_prms[self.nspecies+s]) 
                        atom_masses[s] = float(phys_prms[2*self.nspecies + s])
                    self.HX = mass_fractions
                    self.species_list = species_names
                    self.atomic_masses = atom_masses
                    self.molec_adjust = float(phys_prms[-1])
                elif line[0] == "#bcs:":
                    bcs = [float(x) for x in line[1:]]
                    # boundary conditions
                    self.Rmin = bcs[0]
                    self.Rmax = bcs[1]
                    self.rho_rmin = bcs[2]
                    self.T_rmin = bcs[3]
                    Ys_array = np.zeros(self.nspecies)
                    Ncol_array = np.zeros(self.nspecies)
                    for s in range(self.nspecies):
                        Ys_array[s] = float(bcs[4+s])
                        Ncol_array[s] = float(bcs[4+self.nspecies + s])
                    self.Ys_rmin = Ys_array
                    self.Ncol_sp = Ncol_array
                    self.erf_drop = np.array(bcs[4+2*self.nspecies:4+2*self.nspecies+2])
                elif line[0] == "#tech:":
                    tech = [float(x) for x in line[1:]]
                    # tech parameters
                    self.breezeparam = tech[0]
                    self.rapidity = tech[1]
                    self.erfn = tech[2]
                    self.mach_limit = tech[3]
                elif line[0] == "#flags:":
                    flags = line[1:]
                    # flags
                    self.lyacool = int(flags[0])
                    self.tidalforce = int(flags[1])
                    self.bolo_heat_cool = float(flags[2])
                    self.integrate_outward = int(flags[3])
#                     print(self.integrate_outward)
                elif line[0] == "#add_prms:":
                    add_prms = line[1:]
                    self.n_add_prms = int(add_prms[0])
                    self.add_param_name = [0]
                    self.add_param_val  = [0]
                    if self.n_add_prms>0:
                        add_param_name_array = []
                        add_param_val_array  = np.zeros(self.n_add_prms)
                        self.add_param_dict = {}
                        for i in range(self.n_add_prms):
                            add_param_name_array = np.append(add_param_name_array,add_prms[2*i])
                            add_param_val_array[i]  = float(add_prms[2*i+1])
                            self.add_param_dict[add_prms[2*i]] = float(add_prms[2*i+1])
                        self.add_param_name = add_param_name_array
                        self.add_param_val  = add_param_val_array
                elif line[0] == "#spec_prms:": #moved up to allow access to nspecies, may cause read problems
                    spec = [x for x in line[1:]]
                    self.npts = int(spec[0])
#                     self.nspecies = int(spec[1])
                    self.spec_date = spec[2]
                    self.spec_src_file = spec[3]
                    self.spec_kind = spec[4]
                    self.spec_window = np.asarray([float(spec[5]),
                                                   float(spec[6])])
                    self.spec_resolved = np.asarray([float(spec[7]),
                                                     float(spec[8])])
                    self.spec_normalized = np.asarray([float(spec[9]),
                                                       float(spec[10])])
                    if (self.spec_kind != 'full' and
                        self.spec_window[0] != self.spec_window[1]):
                        print('WARNING: Source is monochromatic, '
                              'but window bounds not equal.\n'
                              f'         window: [{self.spec_window[0]}, '
                              f'{self.spec_window[1]}].')
                    self.ion_pot = np.zeros(self.nspecies)
                    for s in range(len(self.ion_pot)):
                        self.ion_pot[s] = float(spec[11+s])
                    sft = s
                    self.E_wl = np.zeros(self.npts)
                    self.wPhi_wl = np.zeros(self.npts)
                    self.sigma_wl = np.zeros((self.npts, self.nspecies))
#                 elif line[1] == "$hc/\\lambda_i$":
#                     self.species = ['']*self.nspecies
#                     for s in range(self.nspecies):
#                         self.species[s] = f'{line[4*(s+1)].split("}")[0]}'
                    self.species = (self.species_list).copy()
                elif line[0] == "##":
                    spec_data = [float(x) for x in line[1:]]
                    self.E_wl[i_spec] = spec_data[0]
                    self.wPhi_wl[i_spec] = spec_data[1]
                    for s in range(self.nspecies):
                        self.sigma_wl[i_spec][s] = spec_data[2+s]
                    i_spec += 1
                elif line[0][0] == "#": #Skipping other comments
                    continue
                else: # Assuming all comment strings are in header
                    if i_spec != self.npts:
                        print("ERROR: Failed reading simulations' spectrum.\n"
                              f'       Only read in {i_spec} points, '
                              f'but should have read in {self.npts} points.')
                    else:
                        word = ''
                        for s in range(self.nspecies):
                            word += f'sigma_{self.species_list[s]},'
                        self.sim_spectrum = pd.DataFrame(
                            np.column_stack((self.E_wl, self.wPhi_wl,
                                             self.sigma_wl)),
                            columns=['E', 'wPhi', *word.split(',')[:-1]])
                    break
        # Read in solution data as pandas data frame
        self.soln_norm = pd.read_csv(file, comment='#', header=None,
                                     names=self.varz)
        self.soln = self.soln_norm*self.scales
        
        # Mean molecular weight
        if self.species_list[0] != 'HI':
            print("WARNING: Calculation of dimensionless mean molecular weight, mu, assumes Hydrogen is the first species listed.")
#         mu_denominator = 0
#         for j in range(self.nspecies):
#             col_name = 'Ys_'+self.species_list[j]
#             mu_denominator += (self.atomic_masses[0] / self.atomic_masses[j]) *self.HX[j]*(2-self.soln[col_name]) 
#         self.soln['mu'] = const.mH/mu_denominator
        # Make tuples of input files for ease of use
        self.planet_tuple = [self.Mp, self.Rp, self.Mstar, self.semimajor,
                             self.Ftot,self.Lstar]
        self.physics_tuple = [self.HX, self.species_list, self.molec_adjust, self.atomic_masses] 
        self.bcs_tuple = [self.Rmin, self.Rmax, self.rho_rmin, self.T_rmin,
                          self.Ys_rmin, self.Ncol_sp, self.erf_drop]
        self.flags_tuple = [self.lyacool, self.tidalforce,
                            self.bolo_heat_cool, self.integrate_outward]
        self.tech_tuple = [self.breezeparam, self.rapidity, self.erfn,
                           self.mach_limit]
        self.add_tuple = [self.n_add_prms,self.add_param_name,self.add_param_val]
        self.spectrum_tuple = [self.npts, self.nspecies, self.spec_date,self.spec_src_file,
                               self.spec_kind, self.spec_window,
                               self.spec_resolved, self.spec_normalized,
                               self.ion_pot, self.E_wl, self.wPhi_wl,
                               self.sigma_wl, self.species]
        
        # set surface Rp values (may not be 0 if integrated inwards)
        self.rmin_index = len(self.soln[self.soln_norm['r'] < self.Rmin])
        # set critical values
        self.R_sp = 1.+self.soln_norm['z'][0] #sonic point radius, critical point BC in Parker Wind
        self.crit_index = len(self.soln[self.soln_norm['r'] < self.R_sp])-1
        self.v_crit = self.soln['v'][self.crit_index]
        # Mass loss
        self.calc_massloss()
        # Skip unnecessary calculations when users is expediting windsoln
        if not expedite:
            # Add user vars to dataframe
            self.add_user_vars()
            # collisionality calculations
#             self.Kn_hb_crit = self.soln['Kn_hb'][self.crit_index]
#             self.Kn_Co_crit = self.soln['Kn_Co'][self.crit_index]
#             self.calc_r_exo()
#             self.calc_tau_one()
#             self.calc_Jeans()
            if self.integrate_outward == 1:
#                 print('integrating outward')
                self.calc_fits() 
                self.calc_Coriolis()
                self.calc_vert_extent()
                self.calc_ballistic()
        elif add_uservars:
            self.add_user_vars()
           

    def add_user_vars(self):
        self.soln = self.soln.drop(columns=self.soln.columns[6+2*self.nspecies:])

#         start = time.time()
        fe = np.genfromtxt('wrapper/wrapper_utils/recombo_vars.dat',delimiter=',',usecols=(0,1,2,3))[0:10]
        rrec = np.genfromtxt('wrapper/wrapper_utils/recombo_vars.dat',delimiter=',',skip_header=10,usecols=(0,1,2,3))[0:349]
        rnew = np.genfromtxt('wrapper/wrapper_utils/recombo_vars.dat',delimiter=',',skip_header=359,usecols=np.arange(6))
        Verner = pd.read_csv('McAstro/atoms/Verner.csv', comment='#')
        def alpha_rec(species,T):
            '''Adapted from the CLOUDY rrfit fortran agorithm. Computes the temperature-dependent recombination coefficient
            for a species with Z protons and N_e electrons.

            Parameters:
            species: str, species name in any format, e.g., 'C II' or 'h1'
            T: float, temperature of the gas in K
            '''
            species_name_spaced = McAtom.formatting_species_list([species])[0]
            ma = McAtom.atomic_species(species_name_spaced)
            iz = ma.Z #float((Verner['Z']).item())
            iN = ma.Ne #float((Verner['Ne']).item())
            T = np.array(T)
            if iN<3 or iN==11 or (iz>5 and iz<9) or iz==10:
                idx = np.where((rnew[:,0]==iz) & (rnew[:,1]==iN) )
                tt = np.sqrt(T/rnew[idx,4])
                r = rnew[idx,2] / ( tt*np.power(tt+1.0,1.0-rnew[idx,3]) * np.power(1.0+np.sqrt(T/rnew[idx,5]),1.0+rnew[idx,3]) );
            else:
                tt=T*1.0e-04
                if iz==26 and iN<13:
                    idx = np.where(fe[:,0]==iN)
                    r=fe[idx,1]/np.power(tt,(fe[idx,2]+fe[idx,3]*np.log10(tt)));
                else:
                    idx = np.where((rrec[:,0]==iz) & (rrec[:,1]==iN))
                    r=rrec[idx,2]/np.power(tt,rrec[idx,3])

            return r[0][0] 
#         start = time.time()
        self.gamma = 5./3. #Read this in from somewhere?
        self.Omega = np.sqrt(const.G*(self.Mp+self.Mstar)/self.semimajor**3)
        self.semimajor_normed = self.semimajor/self.Rp
        self.rHill = (self.semimajor_normed*(self.Mp/(3.*self.Mstar))**(1./3.))
        self.v_esc = np.sqrt(2*const.G*self.Mp/self.Rp)
        # Calculate sim_spectrum totals/means
        self.Phi0 = self.Ftot*self.sim_spectrum['wPhi'].sum()
        self.F0 = self.Ftot*(self.sim_spectrum['E']
                            *self.sim_spectrum['wPhi']).sum()
        self.E0 = self.F0/self.Phi0
        self.wl0 = const.hc/self.E0
        self.nu0 = const.c/self.wl0

        # number densities
        species_copy = self.species_list.copy()
        spaced = McAtom.formatting_species_list(species_copy)
        n_tot = np.zeros_like(self.soln['rho'])
        n_tot_neutral = np.zeros_like(self.soln['rho'])
        self.soln['n_e'] = np.zeros_like(self.soln['rho'])
        for j in range(self.nspecies):
            element_name = ((spaced[j]).split())[0]
            lowest_state = ((spaced[j]).split())[1] #will have to adapt for elements with more than 1 ionization state
            #converting to arabic numbers to make future multiple-ionization-state version of the code easier
            highest_state = McAtom.arabic_to_roman(McAtom.roman_to_arabic(lowest_state)+1)
            total = 'n_'+element_name
            neutral = 'n_'+element_name+lowest_state
            ionized = 'n_'+element_name+highest_state
            self.soln[total] = self.HX[j]*self.soln['rho']/self.atomic_masses[j]
            self.soln[neutral] = self.soln['Ys_'+(self.species_list[j]).replace(' ','')]*self.soln[total]
            self.soln[ionized] = self.soln[total] - self.soln[neutral]
            self.soln['n_e'] += self.soln[ionized] #print and see if same
            n_tot += self.soln[total]
            n_tot_neutral += self.soln[neutral]
        self.soln['n_tot'] = n_tot
        
        # Mean molecular weight
        if (self.species_list[0]).replace(' ','') != 'HI':
            print("WARNING: Calculation of dimensionless mean molecular weight, mu, assumes Hydrogen is the first species listed.")
        mu_denominator = 0
        for j in range(self.nspecies):
            col_name = 'Ys_'+(self.species_list[j]).replace(' ','')
            mu_denominator += (self.atomic_masses[0] / self.atomic_masses[j]) *self.HX[j]*(2-self.soln[col_name]) 
        
        v = self.soln_norm['v']
        #Approximate pressure scaleheight in units of Rp at the base of the sim 
        Hsc0 = const.kB*self.T_rmin*1e4*self.Rp/(const.mH*2.3*const.G*self.Mp)
        self.smoothing_erf = np.zeros_like(v)
        self.smoothing_erf_mu = np.zeros_like(v)
        
        for i in range(len(self.smoothing_erf)):
            rate_mu = 10/Hsc0
            rate = 1/Hsc0
            self.smoothing_erf[i] = 1-math.erf((v[i]-self.erf_drop[0])/self.erf_drop[1])
            self.smoothing_erf_mu[i] = 1-math.erf((v[i]-self.erf_drop[0])/self.erf_drop[1])
            
#         print(np.max(self.smoothing_erf))
        self.smoothing_erf /= np.max(self.smoothing_erf)
        self.smoothing_erf_mu /= np.max(self.smoothing_erf_mu)
    
        
        #erf set above. Drops mu from molec_adjust (usually 2.3 mH) below wind to mu inside wind
        mu_atom = const.mH/mu_denominator*(1-self.smoothing_erf_mu) #(usually) atomic within wind
        mu_mol = self.molec_adjust*self.smoothing_erf_mu*const.mH #molecular below wind
        
        self.soln['mu'] = mu_atom+mu_mol
        
        # Gas pressure
        self.soln['P'] = (self.soln['rho']
                          *const.kB*self.soln['T']/self.soln['mu'])        
        # Ram pressure
        self.soln['ram'] = (0.5*self.soln['rho']*self.soln['v']**2+
                            self.soln['P'])
        # velocities
        self.soln['cs'] = np.sqrt(self.gamma*self.soln['P']/self.soln['rho'])
        self.soln['Mach'] = self.soln['v']/self.soln['cs']
        
        #Pressure scaleheight
        self.soln['Hsc'] = const.kB*self.soln['T']*(self.soln['r']**2) 
        self.soln['Hsc'] /= (self.soln['mu']*const.G*self.Mp)         
       
        # Multifrequency calculations
        background_ioniz_frac = self.soln['n_HII']/n_tot #only based on H fraction
        background_ioniz_frac[background_ioniz_frac<0] = 1e-10
        
        #Re-written for speed, so most is reshaping 
        ionization_rate = np.zeros((len(self.soln), self.nspecies))
        euv_ionization_rate = np.zeros((len(self.soln), self.nspecies))
        ionization_tau = np.zeros((len(self.soln), self.nspecies))
        ionization_eq = np.zeros((len(self.soln), self.nspecies))
        heating_rate = np.zeros((len(self.soln), self.nspecies))
        total_heating = np.zeros(len(self.soln))
        rec_coeff = np.zeros((len(self.soln), self.nspecies))

        Ncol_arr    = self.soln.iloc[:,4+self.nspecies:4+2*self.nspecies]
        Ncol_arr[Ncol_arr<0] = 0      #less than 0 because unconverged soln returns negative Ncols. This feels sus.
        sigmas = self.sim_spectrum.iloc[:,2:2+self.nspecies]
        taus = np.dot(Ncol_arr,sigmas.T)
        wPhi = np.multiply(np.tile(self.sim_spectrum['wPhi'],(len(taus),1)),np.exp(-taus))*self.Ftot

        # Adapted from Mocassin (Shull & Steenberg 1985)
        # Accounts for secondary ionizations due to highly energetic (>100eV) incoming photons
        frac_in_heat = 0.9971 * (1 - pow(1-pow(background_ioniz_frac,0.2663),1.3163))
        for s, species in enumerate(self.species): #48s
            E_matrix = np.tile(self.E_wl,(len(background_ioniz_frac),1))
            Ncol_matrix = np.tile(Ncol_arr.iloc[:,s],(self.npts,1)).T
            sig_matrix = np.tile(sigmas.iloc[:,s],(len(background_ioniz_frac),1))
            
            f = np.nan_to_num(sig_matrix*Ncol_matrix / taus) #frac of incoming photon energy that will interact with species s
            E0_matrix = E_matrix - self.ion_pot[s]

            heatfrac_matrix = (np.tile(frac_in_heat,(len(self.E_wl),1))).T
            heatfrac_matrix[E0_matrix<6.408707e-11] = 1 #zeroing where E_0 too low for secondary ionizations (100eV)
            #technically should go to lower eV, but Shull & Steenberg analytic eq. is only valid down to 100 eV
            for m in range(self.nspecies):
                sigma_m = sigmas.iloc[:,m]
                sig_matrix_m = np.tile(sigmas.iloc[:,m],(len(background_ioniz_frac),1))
                sig_matrix_m_Hydrogen = np.tile(sigmas.iloc[:,0],(len(background_ioniz_frac),1))
                if self.species[m] == 'HI':
                    C=0.3908; a=0.4092; b=1.7592
                    frac_in_ion = C * np.power(1-np.power(background_ioniz_frac,a),b) 
                    ionfrac_matrix = (np.tile(frac_in_ion,(len(self.E_wl),1))).T
                    ionfrac_matrix[E0_matrix<1.6022e-10] = 0 #zeroing where E_0 too low for secondary ionizations
                    ionfrac_matrix = np.nan_to_num(ionfrac_matrix) 
                    #Number of secondary ionizatioms (as a func of background ioniz. frac.)
                    eta = ionfrac_matrix*E0_matrix/self.ion_pot[m]    
                    #see Murray-Clay (2009) for energy equation. Here we add the factor eta per Shull & Van Steenberg (1985)
                    secondary_H_ion_rate = eta*sig_matrix_m*f*wPhi
                    ionization_rate[:,m] += np.sum(secondary_H_ion_rate,axis=1)
                elif self.species[m]=="HeI":
                    C=0.0554; a=0.4614; b=1.6660
                    frac_in_ion = C * np.power(1-np.power(background_ioniz_frac,a),b) 
                    ionfrac_matrix = (np.tile(frac_in_ion,(len(self.E_wl),1))).T
                    ionfrac_matrix[E0_matrix<1.6022e-10] = 0 #zeroing where E_0 too low for secondary ionizations
                    ionfrac_matrix = np.nan_to_num(ionfrac_matrix) 
                    eta = ionfrac_matrix*E0_matrix/self.ion_pot[m]    
                    ionization_rate[:,m] += np.sum(eta*sig_matrix_m*f*wPhi,axis=1)      
                else: 
                    #To lowest order, can estimate metal ioniz. rate as H ioniz. rate x sigma[metal]/sigma[H]                   
                    ionization_rate[:,m] += np.sum(secondary_H_ion_rate*(sig_matrix_m/sig_matrix_m_Hydrogen),axis=1)           
            heating_rate[:,s] = np.sum((heatfrac_matrix)*E0_matrix*sig_matrix*f*wPhi,axis=1)
            ionization_rate[:,s] += np.sum( sig_matrix*f*wPhi,axis=1) 
            euv_ionization_rate[:,s] = np.sum( sig_matrix*wPhi,axis=1) #optically thin ion rate of H. limiting to EUV is negligible, no fractionation into multiple species because optically thin

            n_abs = self.soln['n_'+species]
            with np.errstate(all='ignore'):
                ionization_tau[:,s] = (
                np.where(ionization_rate[:,s] != 0,
                1./ionization_rate[:,s], np.nan)
                )
            alpha = alpha_rec(species,self.soln['T'])
            rec_coeff[:,s] = alpha
            self.alpha=alpha
            temp = ionization_rate[:,s]/(n_abs*alpha) 
            with np.errstate(all='ignore'):
                ionization_eq[:,s] = (
                    np.where(temp >=1e-200,
                             0.5*temp*(np.sqrt(1.+4./temp)-1.), 0)
                )
            heating_rate[:,s] *= n_abs #/rho temporary
            ionization_rate[:,s] *= n_abs
            euv_ionization_rate[:,s] *= n_abs
            

        total_heating = np.sum(heating_rate,axis=1)

        ion_header = ''
        euv_ion_header = ''
        ion_eq_header = ''
        tau_ion_header = ''
        heat_header = ''
        for s in self.species_list:
            ion_header += f"ionization_{s.replace(' ','')},"
            euv_ion_header += f"euv_ionization_{s.replace(' ','')},"
            ion_eq_header += f"ion_eq_{s.replace(' ','')},"
            tau_ion_header += f"tau_ion_{s.replace(' ','')},"
            heat_header += f"heating_{s.replace(' ','')},"
        ion_header = ion_header.split(',')[:-1]
#         print(ion_header)
#         print(ionization_rate)
        euv_ion_header = euv_ion_header.split(',')[:-1]
        ion_eq_header = ion_eq_header.split(',')[:-1]
        tau_ion_header = tau_ion_header.split(',')[:-1]
        heat_header = heat_header.split(',')[:-1]
        self.soln = self.soln.join(
            pd.DataFrame(ionization_rate, columns=ion_header)
        )
        self.soln = self.soln.join(
            pd.DataFrame(euv_ionization_rate, columns=euv_ion_header)
        )
        self.soln = self.soln.join(
            pd.DataFrame(heating_rate, columns=heat_header)
        )
        self.soln = self.soln.join(
            pd.DataFrame(ionization_tau, columns=tau_ion_header)
        )
        self.soln = self.soln.join(
            pd.DataFrame(ionization_eq, columns=ion_eq_header)
        )
        self.soln['heat_ion'] = total_heating
        self.I_mean = 0
        for s,sp in enumerate(self.species_list):
            s_ion = f"ionization_{sp.replace(' ','')}"
            self.I_mean += (self.ion_pot[s]*integrate.simps(self.soln[s_ion],
                                                            self.soln['r']))
        # Efficiency 1 based on F0 guess, we conserve either only Phi0 or F0
        self.eff1 = 1.-self.I_mean/self.F0
        self.I_mean /= self.Phi0
        # Efficiency 2 represents the simulations total heating
        self.eff2 = integrate.simps(self.soln['heat_ion'],
                                    self.soln['r'])/self.F0
        # Knudsen number (collisionality) - Only computed for H
        ys_HI = 'Ys_'+(self.species_list[0]).replace(' ','')
        self.soln['DlnP'] = (1./np.gradient(np.log(self.soln['P']),
                                            -self.soln['r'])) #1/press. scaleheight
        # Hardbody collisions (3 angstroms) - area ~ 1e-15 cm2
        self.soln['mfp_hb'] = (1./(self.soln['rho']/const.mH      #FIX
                                   *np.pi*(2.*1e-8)**2))
        # coloumb cross section (see RMC 2009)
        self.soln['mfp_Co'] = (1./(self.soln['n_HII']       #this should be exobase
                                   *1e-13*(self.soln['T']/1e4)**-2))
        # Weighted mix of Coloumb and hardbody collisions
        self.soln['mfp_mx'] = (1./((1.-self.soln[ys_HI])/self.soln['mfp_Co']+
                                   self.soln[ys_HI]/self.soln['mfp_hb'])) 
        self.soln['Kn_hb'] = self.soln['mfp_hb']/self.soln['DlnP'] #divided by scalelength of system
        self.soln['Kn_Co'] = self.soln['mfp_Co']/self.soln['DlnP']
        self.soln['Kn_mx'] = self.soln['mfp_mx']/self.soln['DlnP']
#         add mfp by n1sig1 + n2sig2+ n3.... Collision rates from 
        #For neutral and for ions for each species
#         # Ionization
#         # Not traditional collisional. Instead, it's how far neutral goes before ionizing
#         for s, tau_name in enumerate(tau_ion_header):
#             L_ion_name = f'L_ion_{self.species[s]}'
#             Kn_ion_name = f'Kn_ion_{self.species[s]}'
#             self.soln[L_ion_name] = self.soln['v']*self.soln[tau_name] #dist before ionizing
#             self.soln[Kn_ion_name] = self.soln[L_ion_name]/self.soln['DlnP']
        # Rates
        for s,species in enumerate(self.species):
            alpha = rec_coeff[:,s]
            #recombination = alpha*ne*n+
            #n+ = n_tot*(1-Ys) = (n_0/Ys)*(1-Ys) <--in this form for ease of coding 'n_'+species
            self.soln['recomb_'+species] = alpha*self.soln['n_e']*self.soln['n_'+species]*(1-self.soln['Ys_'+species])/self.soln['Ys_'+species]
        # generic cooling
        ## lyman-alpha cooling
        if (self.species_list[0]).replace(' ','') != 'HI':
            print("WARNING: Calculation of Lyman alpha cooling assumes Hydrogen is the first species listed.")        
        self.soln['cool_lyman'] = -7.5e-19*self.soln['n_HI']*self.soln['n_e']*np.exp(-118348/self.soln['T']) #is this too many elextrons when we have other species? Should there be in some limit <= n_H? 
        
        C = pd.read_csv('wrapper/wrapper_utils/line_cooling_coeffs.dat',comment='#',delimiter=' ')
        for s,species in enumerate(self.species):
            ne = self.soln['n_e']
            T = self.soln['T']
            if species == 'CI':  #CII line cooling
#                 print('CII is triggered')
                line_names = ['cool_CII_1570000A','cool_CII_2326A','cool_CII_1334A']
                start = 8; stop = 11
                A, T_line,nc = C['A'][start:stop], C['T_line'][start:stop], C['n_c'][start:stop]
                nIONj = self.soln['n_CII']
                for idx in range(len(A.T)):
                    self.soln[line_names[idx]] = -nIONj*ne* A.iloc[idx]* np.exp(-T_line.iloc[idx]/T) / (ne*(1+nc.iloc[idx]/ne))   
            if species == 'CII': #CIII line cooling
                line_names = ['cool_CIII_1910A','cool_CIII_977A']
                start = 11; stop = 13
                A, T_line,nc = C['A'][start:stop], C['T_line'][start:stop], C['n_c'][start:stop]
                nIONj = self.soln['n_CIII']
                for idx in range(len(A.T)):
                    self.soln[line_names[idx]] = -nIONj*ne* A.iloc[idx]* np.exp(-T_line.iloc[idx]/T) / (ne*(1+nc.iloc[idx]/ne)) 
            if species == 'OI':  #OII line cooling
                line_names = ['cool_OII_834A','cool_OII_2741A','cool_OII_3727A','cool_OII_7320A']
                start = 0; stop = 4
                A, T_line,nc = C['A'][start:stop], C['T_line'][start:stop], C['n_c'][start:stop]
                nIONj = self.soln['n_OII']
                for idx in range(len(A.T)):
                    self.soln[line_names[idx]] = -nIONj*ne* A.iloc[idx]* np.exp(-T_line.iloc[idx]/T) / (ne*(1+nc.iloc[idx]/ne))
            if species == 'OII': #OIII line cooling
                line_names = ['cool_OIII_520000A','cool_OIII_5000A','cool_OIII_166A','cool_OIII_84A']
                start = 4; stop = 8
                A, T_line,nc = C['A'][start:stop], C['T_line'][start:stop], C['n_c'][start:stop]
                nIONj = self.soln['n_OIII']
                for idx in range(len(A.T)):
                    self.soln[line_names[idx]] = -nIONj*ne* A.iloc[idx]* np.exp(-T_line.iloc[idx]/T) / (ne*(1+nc.iloc[idx]/ne))
        
        ## PdV work
        self.soln['cool_PdV'] = self.soln['P']*self.soln['v']/self.soln['rho']
        self.soln['cool_PdV'] *= np.gradient(self.soln['rho'], self.soln['r'])
        ## conduction
        kappa = 4.45e4*(self.soln['T']/1e3)**(0.7)
        grad_T = np.gradient(self.soln['T'], self.soln['r'])
        grad_T[0] = grad_T[1]
        self.soln['cool_cond'] = kappa*np.gradient(grad_T, self.soln['r'])
        self.soln['cool_cond'] += (grad_T*(kappa/self.soln['r']
                                           +np.gradient(kappa, self.soln['r'])))
#         ## collisional cooling   
#         self.soln['cool_col'] = -1.3e-21*(self.HX[0]/const.mH)**2       #FIX - need collisional cooling for all
#         self.soln['cool_col'] *= (np.exp(-157809/self.soln['T'])
#                                   *np.sqrt(self.soln['T'])*(self.soln['rho'])**2
#                                   *(1.-self.soln['Ys_HI'])*self.soln['Ys_HI'])
        
        ## recombination cooling                       
#         self.soln['cool_rec'] = -2.85e-27*(self.HX[0]/const.mH)**2           #FIX
#         self.soln['cool_rec'] *= (np.sqrt(self.soln['T'])
#                                   *(5.914-0.5*np.log(self.soln['T'])
#                                     +0.01184*(self.soln['T'])**(1./3.))
#                                   *(self.soln['rho']*(1.-self.soln['Ys_HI']))**2)
        #CHECK (may be diff prefactor for each species)
        self.soln['cool_rec'] = -2.85e-27*(self.soln['n_e'])**2 *np.sqrt(self.soln['T'])* (5.914-0.5*np.log(self.soln['T']) +0.01184*(self.soln['T'])**(1./3.))
#         self.soln['cool_rec'] *= np.sqrt(self.soln['T'])
#                                   *(5.914-0.5*np.log(self.soln['T'])
#                                     +0.01184*(self.soln['T'])**(1./3.))
         #n_e squared because num.ionized = num.electrons
    
        ## Bolometric heating and cooling #FIX should be read in from somewhere
        self.soln['boloheat'] = np.zeros_like(self.soln['rho'])
        self.soln['bolocool'] = np.zeros_like(self.soln['rho'])
        
        if self.bolo_heat_cool != 0:
            kappa_opt = 4e-3*self.smoothing_erf; kappa_IR = 1e-2*self.smoothing_erf  #need pressure BC (smooth using erf) because kappas are not valid in wind
            F_opt = self.Lstar/(4*np.pi*self.semimajor**2)

            self.soln['boloheat'] = F_opt*(kappa_opt+0.25*kappa_IR)  * self.soln['rho']
            self.soln['bolocool'] = -2*const.sig_SB*self.soln['T']**4 * kappa_IR   * self.soln['rho']
    
        
#         ## free-free emission
#         self.soln['cool_free'] = -1.426e-27*1.3*(self.HX[0]/const.mH)**2          #FIX ,. Make sure that ne-ne is actually what doms
#         self.soln['cool_free'] *= (np.sqrt(self.soln['T'])
#                                    *(self.soln['rho']*(1.-self.soln['Ys_HI']))**2)
        ## gravitational cooling
        self.soln['cool_grav'] = -const.G*self.Mp
        self.soln['cool_grav'] *= self.soln['v']*self.soln['rho']/self.soln['r']**2
        ## advection
        self.soln['e_therm'] = self.soln['P']/((self.gamma-1.)*self.soln['rho'])
        self.soln['heat_advect'] = self.soln['v']*self.soln['rho']
        grad = -np.gradient(self.soln['e_therm'],self.soln['r'])
#         grad[0] = grad[1]
        self.soln['heat_advect'] *= grad
        # cumulative differential heating
        self.soln['cum_heat'] = integrate.cumtrapz(
            (self.soln['heat_ion']+self.soln['cool_lyman']) 
            /(self.soln['rho']*self.soln['v']), self.soln['r'], initial=0)
        # Bernoulli constant
        self.soln['sp_kinetic'] = 0.5*self.soln['v']**2
        self.soln['sp_enthalpy'] = self.soln['P']/self.soln['rho']
        self.soln['sp_enthalpy'] *= self.gamma/(self.gamma-1.0)
        self.soln['sp_grav'] = -const.G*self.Mp/self.soln['r']
        self.soln['bern_rogue'] = (self.soln['sp_kinetic']
                                   +self.soln['sp_enthalpy']
                                   +self.soln['sp_grav'])
        if self.tidalforce:
            self.soln['sp_grav'] -= (const.G*self.Mstar*(
                1./(self.semimajor-self.soln['r'])
                +0.5*(self.semimajor-self.soln['r'])**2./self.semimajor**3.))
        self.soln['bern_isen'] = ( self.soln['sp_kinetic']
                                  +self.soln['sp_enthalpy']
                                  +self.soln['sp_grav'])
        self.soln['bern'] = self.soln['bern_isen']-self.soln['cum_heat']
        # velocity at infinity for a Rogue planet
        sqrt_term = 2*(-const.G*self.Mp/self.Rp
                       +self.soln['sp_enthalpy'][self.rmin_index]
                       +self.soln['cum_heat'].iloc[-1])
        if sqrt_term >= 0:
            self.v_inf = np.sqrt(sqrt_term)
        else:
            self.v_inf = -0.0
        # Adiabatic profile
        self.soln['ad_prof'] = (1.+(self.soln['sp_grav'][self.rmin_index]-
                                    self.soln['sp_grav'])
                                /self.soln['sp_enthalpy'])
        # static, unionized adiabatic profile
        self.soln['static_prof'] = (1.
                                    + (self.soln['sp_grav'][self.rmin_index]
                                       - self.soln['sp_grav'])
                                    /self.soln['sp_enthalpy'][self.rmin_index])
        self.soln['static_rho'] = (self.soln['rho'][self.rmin_index]
                                   *self.soln['static_prof']**(1./
                                                               (self.gamma-1.)))
        self.soln['static_P'] = (self.soln['P'][self.rmin_index]
                                 *self.soln['static_prof']**(self.gamma/
                                                             (self.gamma-1.)))
        self.soln['static_T'] = (self.soln['T'][self.rmin_index]
                                 *self.soln['static_prof'])
        # Adiabatic exponents analysis including ionization
        ## (see Hansen, Kawaler, & Trumble pg. 137 and 
        ##  http://astro.phy.vanderbilt.edu/~berlinaa/teaching/
        ##         stellar/AST8030_7_thermodynamics.pdf)
        I_kT = 13.6*const.eV/(const.kB*self.soln['T'])
        ## del_ad: adiabatic del used for convection analysis
        nom = 1.0+self.ion_wgt(1.-self.soln['Ys_HI'])*(1.5+I_kT)                        #FIX?
        dom = 2.5+self.ion_wgt(1.-self.soln['Ys_HI'])*(3.75+5*I_kT+I_kT**2)
        self.soln['del_ad'] = nom/dom
        ## Gamma_2: adiabatic exponent for temperature's response to pressure
        ##          changes——re-expression of grad_ad
        self.soln['Gamma_2'] = 1.0/(1.0-self.soln['del_ad'])
        ## Gamma_3: adiabatic exponent for temperature's response to compression
        nom = 2.+2.*self.ion_wgt(1.-self.soln['Ys_HI'])*(1.5+I_kT)
        dom = 3.+2.*self.ion_wgt(1.-self.soln['Ys_HI'])*(1.5+I_kT)**2
        self.soln['Gamma_3'] = 1.+nom/dom
        
        if self.integrate_outward:
            self.calc_fits()
            self.calc_Coriolis()
            self.calc_vert_extent()
            self.calc_ballistic()
        
        return
    
    
    def calc_mu_base(self):
        if self.species_list[0] != "H I":
            print('WARNING: Computation of dimensionless mean molecular mass, mu, assumes that HI is the first species.')
        denominator = 0
        for j in range(len(self.HX)):
            denominator += (self.atomic_masses[0]/self.atomic_masses[j]) * self.HX[j] #assumes H is first element
        mu = const.mH/denominator
        return mu


    def ion_wgt(self, X):
        """
        Function of ionization fraction used for adiabatic exponents.
        """
        return X*(1.-X)/((1.+X)*(2.-X))


    def calc_massloss(self):
        indx = self.crit_index
        self.Mdot = (4*np.pi*self.soln['r'].iloc[indx]**2
                     *self.soln['rho'].iloc[indx]*self.soln['v'].iloc[indx])/3
        return


    def calc_tau_sp(self):
        int_bounds = self.soln_norm.r >= self.soln_norm['z']+1
        self.act_tau_sp = integrate.simps(self.soln['nHI'][int_bounds]*
                                          6.304e-18*(13.6/16)**3,
                                          x=self.soln['r'][int_bounds])
        self.last_tau = np.asarray(self.soln['tau'])[-1]
        return


    def calc_tau_one(self):
        try:
            index = self.soln_norm[self.soln.tau > 1].index[-1]
        except IndexError:
            self.error += 1
            return
        if (index+1 == len(self.soln['tau'])):
            self.error += 1
            return
        dtaudr = ((self.soln.iloc[index]['tau']-
                   self.soln.iloc[index+1]['tau'])/
                  (self.soln_norm.iloc[index]['r']-
                   self.soln_norm.iloc[index+1]['r']))
        self.r_tau1 = (self.soln_norm.iloc[index]['r']+
                       (1.-self.soln.iloc[index]['tau'])/dtaudr)
        self.tau_index = index
        return


    def calc_r_exo(self, Kn='Kn_hb'):
        try:
            index = self.soln[self.soln[Kn] > 1].index[0]
        except IndexError:
            self.r_exo = self.R_sp
            self.exo_index = self.crit_index
            return
        if (index+1 == len(self.soln[Kn])):
            self.r_exo = self.R_sp
            self.exo_index = self.crit_index
            return
        dKndr = ((self.soln.iloc[index][Kn]-
                   self.soln.iloc[index+1][Kn])/
                  (self.soln_norm.iloc[index]['r']-
                   self.soln_norm.iloc[index+1]['r']))
        self.r_exo = (self.soln_norm.iloc[index]['r']+
                       (1.-self.soln.iloc[index][Kn])/dKndr)
        self.exo_index = index
        if self.r_exo > self.R_sp:
            self.r_exo = self.R_sp
            self.exo_index = self.crit_index
        return


    def calc_Jeans(self):
        self.Jeans = ((const.G*self.Mp*self.soln['mu']
                       /(const.kB*self.soln['T']
                         *self.r_exo*self.scales[0]))[self.exo_index])
        v_mp = np.sqrt(2./self.gamma)*self.soln['cs'][self.exo_index]
        self.Mdot_Jeans = (2.*np.sqrt(np.pi)*self.soln['rho'][self.exo_index]
                           *v_mp*(self.Jeans + 1.0)*np.exp(-self.Jeans)
                           *(self.r_exo*self.scales[0])**2)
        return


    def calc_fits(self):
        self.v_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['r'], self.soln['v'], ext=2, k=1)
#         self.T_fit = interpolate.InterpolatedUnivariateSpline(
#             self.soln_norm['r'], self.soln['T'], ext=2, k=1)
        self.T_fit = interpolate.CubicSpline(self.soln_norm['r'],self.soln['T'])
        self.n_HII_fit = interpolate.CubicSpline(self.soln_norm['r'],self.soln['n_HII'])
        self.mu_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['r'], self.soln['mu'], ext=2, k=1)
        # Time since Rmin, for ease of doing temporal integrals
        r_int = np.linspace(self.soln_norm['r'].min(),
                            self.soln_norm['r'].max(), 1000000, endpoint=False)
        stime = (integrate.cumtrapz(1./self.v_fit(r_int), r_int, initial=0)
                 *self.scales[0])
        self.stream_time = interpolate.InterpolatedUnivariateSpline(
            r_int, stime, ext=2, k=1)
        return


    def calc_Coriolis(self):
        dvds_fit = self.v_fit.derivative()
        # Scale to same units as dvds
        Omega = self.Omega*self.scales[0]

        def cori_dvds(s, y):
            vx, vy = y
            v = self.v_fit(s)
            dvds = dvds_fit(s)
            dyds = [(vx*dvds+2.*Omega*vy)/v,
                    (vy*dvds-2.*Omega*vx)/v]
            return dyds
        # Start integration from r_sonic, where you "lose" pressure support
        self.cori = integrate.solve_ivp(cori_dvds,
                                        (self.R_sp,
                                         self.soln_norm['r'].iloc[-1]),
                                        [-self.v_fit(self.R_sp), 0],
                                        method='RK45', dense_output=True)
        if self.cori.status:
            print(f'Coriolis ivp failuare: {self.cori.status}')
        else:
            # Calculate the x and y coordinates of streamline
            s_int = np.linspace(self.R_sp, self.soln_norm['r'].iloc[-1],
                                1000000)
            cori_x = (-s_int[0]+integrate.cumtrapz(self.cori.sol(s_int)[0]
                                                   /self.v_fit(s_int),
                                                   x=s_int, initial=0))
            cori_y = (integrate.cumtrapz(self.cori.sol(s_int)[1]
                                         /self.v_fit(s_int),
                                         x=s_int, initial=0))
            fit_cori_x = interpolate.InterpolatedUnivariateSpline(s_int, cori_x,
                                                                  ext=0)
            fit_cori_y = interpolate.InterpolatedUnivariateSpline(s_int, cori_y,
                                                                  ext=0)
            def cori_pos_tuple(s):
                return np.array([fit_cori_x(s), fit_cori_y(s)])
            # Set velocity and position tuple (both return vectors [x, y])
            self.cori_vel = self.cori.sol
            self.cori_pos = cori_pos_tuple
            
            def calc_vel_deflect_one_rad(s, phi_0):
                if phi_0 < np.pi/2.:
                    return np.pi/4.-(np.pi+phi_0
                                     -np.arctan2(-self.cori_vel(s)[1],
                                                 -self.cori_vel(s)[0]))
                else:
                    return np.pi/4.-(phi_0-np.arctan2(self.cori_vel(s)[1],
                                                      self.cori_vel(s)[0]))
            phi_0 = np.arctan2(self.cori_vel(s_int[0])[1],
                               self.cori_vel(s_int[0])[0])
            self.r_cori = optimize.fsolve(calc_vel_deflect_one_rad, self.R_sp,
                                         args=(phi_0))[0]
        return


    def calc_vert_extent(self):
        def verticle_extent(h, e_int_leftover=0.0):
            therm_term = ((self.gamma-1.+2.*(1.-e_int_leftover))
                          /(2*(self.gamma-1.))
                          *self.soln['cs'][self.crit_index]**2
                          /(const.G*self.Mstar))*self.Rp
            q = self.Mp/self.Mstar
            return (therm_term + q*(1./h-1./self.R_sp)
                    +1./np.sqrt(h**2 + self.semimajor_normed**2)
                    -1./np.sqrt(self.R_sp**2 + self.semimajor_normed**2))
        
        self.vert_extent = optimize.fsolve(verticle_extent, self.R_sp)[0]
        return


    def calc_ballistic(self, n_ball=31):
        def ballistic_eoq(t, p, Mp, Mstar, smjr, Omg):
            x, y, z, u, v, w = p
            r_planet_15 = (x**2+y**2+z**2)**(1.5)
            r_star_15 = ((x+smjr)**2+y**2+z**2)**(1.5)
            dpdt = [
                u, v, w,
                -const.G*(Mp*x/r_planet_15+Mstar*(x+smjr)/r_star_15)
                    +Omg**2*(x+smjr)+2.*Omg*v,
                -const.G*y*(Mp/r_planet_15+Mstar/r_star_15)
                    +Omg**2*y-2.*Omg*u,
                -const.G*z*(Mp/r_planet_15+Mstar/r_star_15)
            ]
            return dpdt

        sol = [None,]*n_ball
        v_p = np.sqrt(self.v_inf**2+self.v_esc**2)
        for i, phi in enumerate(np.linspace(0, 2*np.pi, n_ball)):
            rx = self.Rp*self.Rmin*np.cos(phi)
            ry = self.Rp*self.Rmin*np.sin(phi)
            vx = v_p*np.cos(phi)
            vy = v_p*np.sin(phi)
            sol[i] = integrate.solve_ivp(
                ballistic_eoq, [0, 1e7], [rx, ry, 0, vx, vy, 0],
                args=(self.Mp, self.Mstar, self.semimajor, self.Omega),
                method='RK45', dense_output=True)
        self.ballistic_sols = sol
        self.n_ball = n_ball
        return


    def repair_interior(self, csv_file='saves/windsoln.csv'):
        """
        Description:
            Replace inward integration with idealized adiabatic static
            atmosphere
        """
        print("ERROR: Broke in multifrequency update, never repaired, glhf.")
        print("ERROR: Never repaired for multispecies - mib.")
        return
        temp = self.soln_norm.copy(deep=True)
        temp.loc[temp['r'] < self.Rmin, 'rho'] = (
            self.soln.loc[self.soln['r'] < self.Rp, 'static_rho']
            /self.scales[1])
        temp.loc[temp['r'] < self.Rmin, 'v'] = (0.0)
        temp.loc[temp['r'] < self.Rmin, 'T'] = (
            self.soln.loc[self.soln['r'] < self.Rp, 'static_T']/self.scales[3])
        temp.loc[temp['r'] < self.Rmin, 'tau'] = (
            integrate.cumtrapz(self.crosssec/const.mH*self.scales[1]
                               *temp.loc[temp['r'] < self.Rmin, 'rho'][::-1],
                               -self.scales[0]*temp.loc[temp['r']
                                                        < self.Rmin, 'r'][::-1],
                               initial=0)[::-1]
            + temp.loc[temp['r'] == self.Rmin, 'tau'].values[0])
        with open(csv_file, 'w') as file_:
                    #Header
#             self.nspecies = len(self.Ys_rmin)
            file_.write("#nspecies: %d\n" %self.nspecies)
            Ncol_string = ''
            Ys_string = ''
            string = ''
            floats = ''
            for n in range(self.nspecies):
                floats += '{:.17e},'
                string += '{:s},'
                Ncol_string += 'Ncol_'+self.species_list[n]+','
                Ys_string   += 'Ys_'+self.species_list[n]+','
            file_.write("#vars: r,rho,v,T,"+Ys_string+Ncol_string+"q,z\n")
            file_.write('#scales: '+(','.join('%.17e' %i for i in self.scales))+ '\n')
            file_.write('#plnt_prms: '+(','.join('%.17e' %i for i in (self.planet_tuple)))+ '\n')
            file_.write('#phys_prms: '+(floats+string+floats).format(*self.HX,*self.species_list,*self.atomic_masses)[:-1]+'\n')
            file_.write('#bcs: '+(','.join('{:.17e}' for i in range(5+2*self.nspecies))).format(*self.bcs_tuple[:4], *self.Ys_rmin, *self.Ncol_sp, *self.erf_drop)+'\n')
            file_.write('#tech: '+(','.join('%.17e' %i for i in self.tech_tuple))+ '\n')
            file_.write('#flags: '+(','.join('%d' %i for i in self.flags_tuple))+ '\n')
#             file_.write('#spec_prms: {:d},{:d},'.format(self.npts,self.nspecies)+self.spec_date+','+self.spec_kind+
#                     ',{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'.format(*self.spec_window,*self.spec_resolved, 
#                                                                                *self.spec_normalized))
#             file_.write(floats[:-1].format(*self.ion_pot)+'\n')  
#             file_.write('#vars: r,rho,v,T,Ys,Ncol,q,z\n')
#             file_.write('#scales: {:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'
#                         '{:.17e},{:.17e},{:.17e}\n'
#                         .format(self.scales[0], self.scales[1], self.scales[2],
#                                 self.scales[3], self.scales[4], self.scales[5],
#                                 self.scales[6], self.scales[7]))
#             file_.write('#plnt_prms: {:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'
#                         '{:.17e},{:.17e}\n'.format(*self.planet_tuple))
#             file_.write('#phys_prms: {:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'
#                         '{:.17e},{:.17e}\n'
#                         .format(*self.physics_tuple))
#             file_.write('#bcs: {:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},'
#                         '{:.17e}\n'.format(*self.bcs_tuple))
#             file_.write('#tech: {:.17e},{:.17e},{:.17e},{:.17e}\n'
#                         .format(*self.tech_tuple))
#             file_.write('#flags: {:d},{:d},{:d},{:d},{:d},{:d}\n'
#                         .format(*self.flags_tuple))
            temp.to_csv(file_, float_format='%.17e', index=False, header=False)
        return


    def calc_norm_fits(self, degree=1):
        # Sometimes playing with degree helps converge regrids, e.g.,
        # degree=2 has linear dvdr instead of constant of degree=1
        self.r_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['q'], self.soln_norm['r'], ext=2, k=degree)
        self.rho_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['q'], self.soln_norm['rho'], ext=2, k=degree)
        self.v_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['q'], self.soln_norm['v'], ext=2, k=degree)
        self.T_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['q'], self.soln_norm['T'], ext=2, k=degree)
        for j in range(self.nspecies):
            self.Ys_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
                self.soln_norm['q'], self.soln_norm['Ys_'+self.species_list[j]], ext=2, k=degree)
            self.Ncol_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
                self.soln_norm['q'], self.soln_norm['Ncol_'+self.species_list[j]], ext=2, k=degree)
        self.z_q_norm_fit = interpolate.InterpolatedUnivariateSpline(
            self.soln_norm['q'], self.soln_norm['z'], ext=2, k=degree)
        return


    def regrid(self, q_arr=None, simple=True):
        self.calc_norm_fits()
        if simple:
            dq0 = 1e-4
            dq1 = 8e-3
            dq2 = 1e-4

            q0 = 0.00
            q1 = 0.05
            q2 = 0.95
            q3 = 1.00

            # Constant dq between q0 and q1
            q_arr = list(np.arange(q0, q1, dq0))+[q1]
            q = q1
            # Linear rate of change between q1 and q2
            s = 0
            while q < (q2+q1)/2:
                q += (dq1-dq0)*s+dq0
                s = (q-q1)/((q2+q1)/2-q1)
                if s > 1:
                    break
                q_arr.append(q)
            s = 0
            q = q2
            q_arr.append(q)
            while q > (q2+q1)/2:
                os = s
                q -= (dq1-dq2)*s+dq2
                s = 1.-(2*q-(q2+q1))/(q2-q1)
                if s > 1:
                    q += (dq1-dq2)*os+dq2
                    break
                q_arr.append(q)
            # Patch two linear regimes together
            q_arr.sort()
            lq = q_arr[q_arr.index(q)-1]
            uq = q_arr[q_arr.index(q)]
            div = 1
            while True:
                dq = (uq-lq)/div
                if dq > dq1:
                    div += 1
                else:
                    break
            for i in range(div):
                q_arr.append(lq+(i+1)*dq)
            q_arr.sort()

            # Constant dq between q2 and q3
            q_arr = np.concatenate((np.array(q_arr),
                                    np.arange(q2, q3-dq2, dq2)[1:],
                                    np.array([1.0])))
        elif q_arr is None:
            dq0 = 5e-5
            dq1 = 1e-5
            dq2 = 1.07e-2

            q0 = 0.00
            q1 = 0.05
            q2 = 0.95
            q3 = 0.99
            q4 = 1.00

            q_arr = [0]
            # Gaussian rate of change between q0 and q1
            q = q0
            s = 0
            while q < q1/2:
                q += -(dq0-dq1)*np.exp(-(s-1)**2)+dq0
                s = (q-q0)/(q1/2-q0)
                q_arr.append(q)
            s = 1
            while q < q1:
                q += -(dq0-dq1)*np.exp(-(s-1)**2)+dq0
                s = 1-(q-q1/2)/(q1-q1/2)
                if s >= 1:
                    q = q1
                q_arr.append(q)
            # Linear rate of change between q1 and q2
            s = 0
            while q < (q2+q1)/2:
                q += (dq2-dq0)*s+dq0
                s = (q-q1)/((q2+q1)/2-q1)
                q_arr.append(q)
            s = 1
            while q < q2:
                q += (dq2-dq0)*s+dq0
                s = 1-(q-(q2+q1)/2)/(q2-(q2+q1)/2)
                if s > 1:
                    q = q2
                q_arr.append(q)
            # Gaussian rate of change between q2 and q3
            s = 0
            while q < q3:
                q += -(dq0-dq1)*s+dq0
                s = (q-q2)/(q3-q2)
                if s > 1:
                    q = q3
                    break
                q_arr.append(q)
            # Constant rate of change between q3 and q4
            q_arr = np.concatenate(
                (np.array(q_arr), np.linspace(q_arr[-1], q4,
                                          int((q4-(q3+dq1))/dq1+1.99))[1:])
            )

        # Create regrided solution, ensure r[0] = Rmin exactly, cannot allow
        # any rounding error, i.e, 9.99999999999999778e-01 != 1.0 in relaxation
        # r, rho, v, T, Ys, Ncol, q, z
        def Ys_q_norm_fit(q,degree=1):
            Ys_normed_fit = np.zeros((len(q),self.nspecies))
            for j in range(self.nspecies):
                Ys_name = 'Ys_'+self.species_list[j]
                s = interpolate.InterpolatedUnivariateSpline(
                    self.soln_norm['q'], self.soln_norm[Ys_name], ext=2, k=degree)
                Ys_normed_fit[:,j] = s(q)
            return Ys_normed_fit
        
        def Ncol_q_norm_fit(q,degree=1):
            Ncol_normed_fit = np.zeros((len(q),self.nspecies))
            for j in range(self.nspecies):
                Ncol_name = 'Ncol_'+self.species_list[j]
                s = interpolate.InterpolatedUnivariateSpline(
                    self.soln_norm['q'], self.soln_norm[Ncol_name], ext=2, k=degree)
                Ncol_normed_fit[:,j] = s(q)
            return Ncol_normed_fit
        
#         Ncol_normed_fit = np.zeros((len(q_arr),self.nspecies))
#         Ys_normed_fit = np.zeros((len(q_arr),self.nspecies))
#         for j in range(self.nspecies):
#             Ys_normed_fit[:,j] = Ys_q_norm_fit('Ys_'+self.species_list[j], q_arr)
#             Ncol_normed_fit[:,j] = Ncol_normed_fit('Ncol_'+self.species_list[j], q_arr)
        df = pd.DataFrame(np.column_stack((
            #strictly enforce r[0] == Rmin
            np.array([self.Rmin]+list(self.r_q_norm_fit(q_arr[1:]))),
            self.rho_q_norm_fit(q_arr),
            self.v_q_norm_fit(q_arr),
            self.T_q_norm_fit(q_arr),
            Ys_q_norm_fit(q_arr),
            Ncol_q_norm_fit(q_arr),
#             self.Ys_q_norm_fit(q_arr),
#             self.Ncol_q_norm_fit(q_arr),
            q_arr,
            self.z_q_norm_fit(q_arr)
        )), columns=self.soln_norm.columns[:(6+2*self.nspecies)])

        with open('regrid_soln.csv', 'w') as file_:
                                #Header
#             self.nspecies = len(self.Ys_rmin)
            file_.write("#nspecies: %d\n" %self.nspecies)
            Ncol_string = ''
            Ys_string = ''
            string = ''
            floats = ''
            for n in range(self.nspecies):
                floats += '{:.17e},'
                string += '{:s},'
                Ncol_string += 'Ncol_'+self.species_list[n]+','
                Ys_string   += 'Ys_'+self.species_list[n]+','
            file_.write("#vars: r,rho,v,T,"+Ys_string+Ncol_string+"q,z\n")
            file_.write('#scales: '+(','.join('%.17e' %i for i in self.scales))+ '\n')
            file_.write('#plnt_prms: '+(','.join('%.17e' %i for i in (self.planet_tuple)))+ '\n')
            file_.write('#phys_prms: '+(floats+string+floats).format(*self.HX,*self.species_list,*self.atomic_masses)[:-1]+'\n')
            file_.write('#bcs: '+(','.join('{:.17e}' for i in range(5+2*self.nspecies))).format(*self.bcs_tuple[:4], *self.Ys_rmin,*self.Ncol_sp,*self.erf_drop)+'\n')
            file_.write('#tech: '+(','.join('%.17e' %i for i in self.tech_tuple))+ '\n')
            file_.write('#flags: '+(','.join('%d' %i for i in self.flags_tuple))+ '\n')
#             file_.write(('#vars: '+('{:s},'*len(self.soln_norm.columns[:8]))[:-1]
#                          +'\n').format(*self.soln_norm.columns[:8]))
#             file_.write(('#scales: '+('{:.17e},'*len(self.scales))[:-1]
#                          +'\n').format(*self.scales))
#             file_.write(('#plnt_prms: '+('{:.17e},'*len(self.planet_tuple))[:-1]
#                          +'\n').format(*self.planet_tuple))
#             file_.write(('#phys_prms: '+('{:.17e},'*len(self.physics_tuple))[:-1]
#                          +'\n').format(*self.physics_tuple))
#             file_.write(('#bcs: '+('{:.17e},'*len(self.bcs_tuple))[:-1]
#                          +'\n').format(*self.bcs_tuple))
#             file_.write(('#tech: '+('{:.17e},'*len(self.tech_tuple))[:-1]
#                          +'\n').format(*self.tech_tuple))
#             file_.write(('#flags: '+('{:d},'*len(self.flags_tuple))[:-1]
#                          +'\n').format(*self.flags_tuple))
            df.to_csv(file_, header=False, index=False, float_format='%.17e')
        return df