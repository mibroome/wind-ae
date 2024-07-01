#!/bin/env python3
import ChiantiPy as ch
import ChiantiPy.tools.io
import ChiantiPy.core

import numpy as np
import math
from scipy.interpolate import CubicSpline

import McAstro.atoms.atomic_species as McAtom
import McAstro.radiation.lines.lyman_alpha as McLyman
import McAstro.radiation.radiative_transfer as McRT
from wrapper.wrapper_utils import constants as const


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
        self.n_H_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_H'])
        self.n_HI_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_HI'])
        self.n_HII_fit = CubicSpline(self.windsoln.soln_norm['r'],self.windsoln.soln['n_HII'])
       
        return
    
    def ray_trace(self, nu, alpha_nu, r_star, vel_lo=-50, vel_hi=50, constfrac=0, nu0=0, f12=0, A21=0, species='', level=0, abund=0, useOfit = 0, onlyHill = 0, onlyrcrit = 0):
        self.v_arr = np.linspace(vel_lo, vel_hi, 101)*1e5
        self.Obscuration = np.zeros(len(self.v_arr))
        r_star = r_star/self.windsoln.Rp
        # Ray trace through spherical extent of outflow
        if onlyrcrit == 1:
            sphere_extent = min(r_star, self.windsoln.R_sp)
        elif onlyHill == 1:
            sphere_extent = min(r_star, self.windsoln.R_Hill)
        else:
            sphere_extent = min(r_star, self.windsoln.Rmax)
        vert_extent = self.windsoln.vert_extent
        print("vertical:", vert_extent)
        Rmin = self.windsoln.Rmin
        if onlyrcrit == 1:
                b_arr = np.linspace(Rmin, min(r_star, self.windsoln.R_sp), 50)
        elif onlyHill == 1:
            b_arr = np.linspace(Rmin, min(r_star, self.windsoln.R_Hill), 50)
        else:
            b_arr = np.linspace(Rmin, min(r_star, self.windsoln.Rmax), 50)
        tau_arr = np.zeros((len(self.v_arr), len(b_arr)))
        #print(len(self.v_arr))
        #print(len(b_arr))
        if onlyrcrit == 1:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.R_sp, 20)
        elif onlyHill == 1:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.R_Hill, 20)
        else:
            pop_s = np.linspace(Rmin, 0.95*self.windsoln.Rmax, 20)
        temp = self.windsoln.T_fit(pop_s)
        levelpops = np.zeros(len(pop_s))
        ne = self.n_HII_fit(pop_s)
        ne[ne<0] = 0
        ne = np.where(ne==0,1.0,ne)
        nprot = ne
        print("about to populate")
        for k in range(len(pop_s)):
            atom = ch.core.ion(species, temperature=temp[k], eDensity=ne[k], pDensity=nprot[k])
            atom.populate()
            levelpops[k] = atom.Population['population'][0][level]
        print(levelpops)
        for j, v in enumerate(self.v_arr):
            #if j%5 == 0: print("j = ", j)
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
            for i, b in enumerate(b_arr):
                #print("i = ", i)
                if onlyrcrit == 1:
                    s =  np.linspace(-self.windsoln.R_sp, self.windsoln.R_sp, 100)
                elif onlyHill == 1:
                    s =  np.linspace(-self.windsoln.R_Hill, self.windsoln.R_Hill, 100)
                else:
                    s = np.linspace(-self.windsoln.Rmax, self.windsoln.Rmax, 100)
                rs = np.sqrt(b**2 + s**2)
                if onlyrcrit == 1:
                    sim_mask = rs < self.windsoln.R_sp
                elif onlyHill == 1:
                    sim_mask = rs < self.windsoln.R_Hill
                else:
                    sim_mask = rs < self.windsoln.Rmax
                rs = rs[sim_mask]
                #print(len(rs))
                s = s[sim_mask]
                if len(s) < 2:
                    continue
                # Setup ray through spherically symmetric outflow
                my_ray = McRT.rt_ray(len(rs))
                my_ray.nu0 = nu0
                my_ray.f12 = f12
                my_ray.A21 = A21
                my_ray.v_LOS = self.windsoln.v_fit(rs)*(s/rs)
#                print("before")
#                print("rs")
#                print(rs)
#                print("pop_s")
#                print(pop_s)
                my_ray.T = self.windsoln.T_fit(rs)
#                temp = self.windsoln.T_fit(rs)
#                temp2 = self.windsoln.T_fit(pop_s)
#                print("after")
                #ne = self.windsoln.n_HII_fit(rs)
                #ne = np.where(ne==0,1.0,ne)
                #nprot = ne #self.windsoln.n_H_fit(rs)
                if constfrac != 0:
                    my_ray.n_abs = self.n_H_fit(rs)*constfrac
                elif species=='h_1' or species=='o_1':
                    levelpop = np.zeros(len(rs))
                    for ppp in range(len(rs)):
                        #print(ppp)
                        #print(ne[ppp], nprot[ppp], my_ray.T[ppp], species)
                        #if ppp%10 == 0:
                        #if j == 0:
                        #    atom = ch.core.ion(species, temperature=my_ray.T[ppp], eDensity=ne[ppp], pDensity=nprot[ppp])
                        #    atom.populate()
                        #    levelpop_arr[ppp] = atom.Population['population'][0][level]
                            #levelpop[ppp] = atom.Population['population'][0][level]
                            #levelpop = 1.0
                        #else:
                            #levelpop[ppp] = levelpop[ppp-1]
                        #    levelpop = 1.0
                        #levelpop = 1.0
                        bestindex = np.argmin(np.fabs(np.array(pop_s - rs[ppp])))
                        levelpop[ppp] = levelpops[bestindex]
                    #print(levelpop)
                    #my_ray.n_abs = self.windsoln.n_HI_fit(rs)*abund*levelpop
                    my_ray.n_abs = self.n_HI_fit(rs)*levelpop
                    #print(my_ray.n_abs)
                elif species=='o_1':
                    levelpop = np.zeros(len(rs))
                    for ppp in range(len(rs)):
                        bestindex = np.argmin(np.fabs(np.array(pop_s - rs[ppp])))
                        levelpop[ppp] = levelpops[bestindex]
                        if useOfit == 1:
                            my_ray.nabs = self.windsoln.n_OI_fit(rs)*levelpop
                        else:
                            my_ray.n_abs = self.windsoln.n_HI_fit(rs)*abund*levelpop
                elif species!='':
                    levelpop = np.zeros(len(rs))
                    for ppp in range(len(rs)):
                        #if ppp%10 == 0:
                            #atom = ch.core.ion(species, temperature=my_ray.T[ppp], eDensity=ne[ppp], pDensity=nprot[ppp])
                            #atom.populate()
                            #levelpop[ppp] = atom.Population['population'][0][level]
                        #    levelpop = 0.5
                        #else:
                            #levelpop[ppp] = levelpop[ppp-1]
                        #    levelpop = 0.5
                        #levelpop = 0.5
                        bestindex = np.argmin(np.fabs(np.array(pop_s - rs[ppp])))
                        levelpop[ppp] = levelpops[bestindex]
                    #print(levelpop)
                    my_ray.n_abs = self.n_HII_fit(rs)*abund*levelpop
                else:
                    my_ray.n_abs = self.windsoln.n_HI_fit(rs)
                my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
                my_ray.mu = self.windsoln.mu_fit(rs)
                my_ray.I[0] = 1e1
                # Do radiative transfer
                tau_arr[j][i] = McRT.radiative_transfer(nu_obs, my_ray, alpha_nu)
                #print(tau_arr[j][i])
        # Do single ray through constant horizontal tube
        tau_hor = np.zeros(len(self.v_arr))
        s = np.linspace(-sphere_extent, sphere_extent, 100)
        # Setup ray through horizontal tube
        my_ray = McRT.rt_ray(len(s))
        my_ray.nu0 = nu0
        my_ray.f12 = f12
        my_ray.A21 = A21
        my_ray.v_LOS = np.zeros(len(s))
        if constfrac != 0:
            my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_HI'])[-1]*constfrac
        elif species=='h_1' or species=='o_1':
            my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_HI'])[-1]*abund*levelpop[-1]
        elif species=='c_2':
            my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_H'])[-1]*abund*levelpop[-1]
            print("abundance in ray trace:", abund, levelpop[-1])
        elif species!='':
            my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_HII'])[-1]*abund*levelpop[-1]
        else:
            my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_HI'])[-1]
        my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
        my_ray.T = np.ones(len(s))*np.array(self.windsoln.soln['T'])[-1]
        my_ray.mu = np.ones(len(s))*np.array(self.windsoln.soln['mu'])[-1]
        my_ray.I[0] = 1e1
        for j, v in enumerate(self.v_arr):
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
           # Do radiative transfer
            tau_hor[j] = McRT.radiative_transfer(nu_obs, my_ray, alpha_nu)
#        for j, v in enumerate(self.v_arr):
#            for i, b in enumerate(b_arr):
#                if b < sphere_extent:
                    # Add tube contribution
#                    tau_arr[j][i] += (tau_hor[j]
#                                      *(1.-np.sqrt(1.-(b/sphere_extent)**2)))
        # Calculate disk average obscuration of spherical extent + horizontal tube
        h = min(r_star, sphere_extent)
        chi = h/r_star
        zeta = np.sqrt(1.-chi**2)
        eta = min(1., sphere_extent/r_star)
        w_hor = (2./const.pi)*(math.acos(zeta)+chi*zeta) - eta**2
        RH = self.windsoln.semimajor*(self.windsoln.Mp/(3.0*self.windsoln.Mstar))**(1./3.)
        print(RH/self.windsoln.Rp, self.windsoln.Rmax, sphere_extent, r_star)
        print(eta)
        for j in range(len(self.v_arr)):
            # Obscuration from spherical extent
            self.Obscuration[j] = eta**2-np.trapz(2.*b_arr*np.exp(-tau_arr[j]),
                                                  x=b_arr)/r_star**2
            #print(self.Obscuration[j], eta**2, np.trapz(2.*b_arr*np.exp(-tau_arr[j]), x=b_arr)/r_star**2, np.max(b_arr), np.max(tau_arr[j]))
            # Obscuration from horizontel extent
            #self.Obscuration[j] += (1-np.exp(-tau_hor[j]))*w_hor