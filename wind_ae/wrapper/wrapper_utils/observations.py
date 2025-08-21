#!/bin/env python3
"""
BETA
"""

import ChiantiPy as ch
import ChiantiPy.tools.io
import ChiantiPy.core

import numpy as np
from scipy.interpolate import CubicSpline

import wind_ae.McAstro.radiation.lines.lyman_alpha as McLyman
import wind_ae.McAstro.radiation.lines.lyman_alpha as McLine
import wind_ae.McAstro.radiation.radiative_transfer as McRT
import wind_ae.McAstro.atoms.atomic_species as McAtom
from wind_ae.wrapper.wrapper_utils import constants as const

# McRT.rt_ray
# McRT.radiative_transfer
# McLyman.Lya_alpha_nu
# McLyman._Lya_nu

class observation:
    def __init__(self, windsoln):
        self.windsoln = windsoln
        self.alpha_nu = McLyman.Lya_alpha_nu
        self.nu = McLyman._Lya_nu
        self.windsoln.calc_fits()
        self.n_tot_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_H'])
        self.n_neutral_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_HI'])
        self.n_ion_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_HII'])
        self.ne_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_e'])
       
        return
    
    def ray_trace(self, nu, alpha_nu, r_star, species='', continuum=False, vel_lo=-50, vel_hi=50, constfrac=0, nu0=0, f12=0, A21=0, 
                  level=0, useOfit = 0, onlyHill = 0, onlyrcrit = 0):      
        # print(species)
        species_name = McAtom.formatting_species_list([species])[0]
        ss = species_name.split()
        chianti_species = ss[0].lower()+'_'+str(McAtom.roman_to_arabic(ss[1]))   
        n0_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_'+species_name.replace(' ','')])
        
        self.v_arr = np.linspace(vel_lo, vel_hi, 101)*1e5
        self.Obscuration = np.zeros(len(self.v_arr))
        r_star = r_star/self.windsoln.Rp
        # Ray trace through spherical extent of outflow
        if onlyrcrit == 1:
            sphere_extent = min(r_star, self.windsoln.R_sp)
        elif onlyHill == 1:
            sphere_extent = min(r_star, self.windsoln.R_hill)
        else:
            sphere_extent = min(r_star, self.windsoln.Rmax)

        Rmin = self.windsoln.Rmin
        if onlyrcrit == 1:
            b_arr = np.linspace(Rmin, min(r_star, self.windsoln.R_sp), 50)
        elif onlyHill == 1:
            b_arr = np.linspace(Rmin, min(r_star, self.windsoln.R_hill), 50)
        else:
            b_arr = np.linspace(Rmin, min(r_star, self.windsoln.Rmax), 50)
        tau_arr = np.zeros((len(self.v_arr), len(b_arr)))

        if onlyrcrit == 1:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.R_sp, 20)
        elif onlyHill == 1:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.R_hill, 20)
        else:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.Rmax, 20)
        temp = self.windsoln.T_fit(pop_s)
        levelpops = np.zeros(len(pop_s))
        ne = self.ne_fit(pop_s)
        ne[ne<0] = 0
        ne = np.where(ne==0,1.0,ne) #so there is at least 1 electron
        nprot = ne #invalid for metals?
        # print("about to populate")
        for k in range(len(pop_s)):
            #can use real ion xsec
            atom = ch.core.ion(chianti_species, temperature=temp[k], eDensity=ne[k], pDensity=nprot[k])
            atom.populate()
            levelpops[k] = atom.Population['population'][0][level]
        # print(levelpops)
        for j, v in enumerate(self.v_arr):
            #if j%5 == 0: print("j = ", j)
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
            for i, b in enumerate(b_arr):
                #print("i = ", i)
                if onlyrcrit == 1:
                    s =  np.linspace(-self.windsoln.R_sp, self.windsoln.R_sp, 100)
                elif onlyHill == 1:
                    s =  np.linspace(-self.windsoln.R_hill, self.windsoln.R_hill, 100)
                else:
                    s = np.linspace(-self.windsoln.Rmax, self.windsoln.Rmax, 100)
                rs = np.sqrt(b**2 + s**2)
                if onlyrcrit == 1:
                    sim_mask = rs < self.windsoln.R_sp
                elif onlyHill == 1:
                    sim_mask = rs < self.windsoln.R_hill
                else:
                    sim_mask = rs < self.windsoln.Rmax
                rs = rs[sim_mask]
                s = s[sim_mask]
                if len(s) < 2:
                    continue
                # Setup ray through spherically symmetric outflow
                my_ray = McRT.rt_ray(len(rs))
                my_ray.nu0 = nu0
                my_ray.f12 = f12
                my_ray.A21 = A21
                my_ray.v_LOS = self.windsoln.v_fit(rs)*(s/rs)
                my_ray.T = self.windsoln.T_fit(rs)
                levelpop = np.zeros(len(rs))
                for ppp in range(len(rs)):
                    bestindex = np.argmin(np.fabs(np.array(pop_s - rs[ppp])))
                    levelpop[ppp] = levelpops[bestindex]
                my_ray.n_abs = n0_fit(rs)*levelpop
                my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
                my_ray.mu = self.windsoln.mu_fit(rs)
                my_ray.I[0] = 1e1
                # Do radiative transfer
                tau_arr[j][i] = McRT.radiative_transfer(nu_obs, my_ray, alpha_nu, species_name)
        # Do single ray through constant horizontal tube
        tau_hor = np.zeros(len(self.v_arr))
        s = np.linspace(-sphere_extent, sphere_extent, 100)
        # Setup ray through horizontal tube
        my_ray = McRT.rt_ray(len(s))
        my_ray.nu0 = nu0
        my_ray.f12 = f12
        my_ray.A21 = A21
        my_ray.v_LOS = np.zeros(len(s))
        my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
        my_ray.T = np.ones(len(s))*np.array(self.windsoln.soln['T'])[-1]
        my_ray.mu = np.ones(len(s))*np.array(self.windsoln.soln['mu'])[-1]
        my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_'+species_name.replace(' ','')])[-1]*levelpop[-1]
        my_ray.I[0] = 1e1
        for j, v in enumerate(self.v_arr):
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
           # Do radiative transfer
            tau_hor[j] = McRT.radiative_transfer(nu_obs, my_ray, alpha_nu, species_name)
#        for j, v in enumerate(self.v_arr):
#            for i, b in enumerate(b_arr):
#                if b < sphere_extent:
                    # Add tube contribution
#                    tau_arr[j][i] += (tau_hor[j]
#                                      *(1.-np.sqrt(1.-(b/sphere_extent)**2)))
        # Calculate disk average obscuration of spherical extent + horizontal tube
        # h = min(r_star, sphere_extent)
        # chi = h/r_star
        # zeta = np.sqrt(1.-chi**2)
        eta = min(1., sphere_extent/r_star)
        # w_hor = (2./const.pi)*(math.acos(zeta)+chi*zeta) - eta**2
        # RH = self.windsoln.semimajor*(self.windsoln.Mp/(3.0*self.windsoln.Mstar))**(1./3.)
        # print(tau_arr)
        # print(eta)
        # print(b_arr)
        for j in range(len(self.v_arr)):
            # Obscuration from spherical extent
            self.Obscuration[j] = eta**2-np.trapezoid(2.*b_arr*np.exp(-tau_arr[j]),
                                                  x=b_arr)/r_star**2
            #print(self.Obscuration[j], eta**2, np.trapz(2.*b_arr*np.exp(-tau_arr[j]), x=b_arr)/r_star**2, np.max(b_arr), np.max(tau_arr[j]))
            # Obscuration from horizontel extent
            #self.Obscuration[j] += (1-np.exp(-tau_hor[j]))*w_hor