#!/usr/bin/env python3

import numpy as np

import McAstro.radiation.lines.lyman_alpha as McLyman
import McAstro.radiation.radiative_transfer as McRT

# McRT.rt_ray
# McRT.radiative_transfer
# McLyman.Lya_alpha_nu
# McLyman._Lya_nu

class observation:
    def __init__(self, windsoln):
        self.windsoln = windsoln
        self.alpha_nu = McLyman.Lya_alpha_nu
        self.nu = McLyman._Lya_nu
        return
    
    
    def ray_trace(self, nu, alpha_nu, vel_lo, vel_hi, r_star):
        self.Obscuration = np.zeros(len(v_arr))
        v_arr = np.linspace(-50, 50, 101)*1e5
        r_star = self.r_star/self.windsoln.Rp
        # Ray trace through spherical extent of outflow
        b_arr = np.linspace(1.0, min(r_star, self.windsoln.Rmax), 50)
        tau_arr = np.zeros((len(v_arr), len(b_arr)))
        for j, v in enumerate(v_arr):
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
            for i, b in enumerate(b_arr):
                s = np.linspace(-self.windsoln.Rmax, self.windsoln.Rmax, 100)
                rs = np.sqrt(b**2 + s**2)
                sim_mask = rs < self.windsoln.Rmax
                rs = rs[sim_mask]
                s = s[sim_mask]
                if len(s) < 2:
                    continue
                # Setup ray through spherically symmetric outflow
                my_ray = McRT.rt_ray(len(rs))
                # Get velocity, and project into LOS (s/rs)
                my_ray.v_LOS = self.windsoln.v_fit(rs)*(s/rs)
                my_ray.n_abs = self.windsoln.n_HI_fit(rs)
                my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
                my_ray.T = self.windsoln.T_fit(rs)
                my_ray.mu = self.windsoln.mu_fit(rs)
                my_ray.I[0] = 1e1
                # Do radiative transfer
                tau_arr[j][i] = McRT.radiative_transfer(nu_obs, my_ray,
                                                        alpha_nu)
        # Do single ray through constant horizontal tube
        tau_hor = np.zeros(len(v_arr))
        sphere_extent = max(r_star, self.windsoln.Rmax)
        s = np.linspace(-sphere_extent, sphere_extent, 100)
        # Setup ray through horizontal tube
        # Take LOS velocity to be zero, and constant values in tube
        my_ray = McRT.rt_ray(len(s))
        my_ray.v_LOS = np.zeros(len(s))
        my_ray.n_abs = np.ones(len(s))*np.array(self.windsoln.soln['n_HI'])[-1]
        my_ray.ds *= (s[1]-s[0])*self.windsoln.Rp
        my_ray.T = np.ones(len(s))*np.array(self.windsoln.soln['T'])[-1]
        my_ray.mu = np.ones(len(s))*np.array(self.windsoln.soln['mu'])[-1]
        my_ray.I[0] = 1e1
        for j, v in enumerate(v_arr):
            nu_obs = nu*np.sqrt((1.-v/const.c)/(1.+v/const.c))
            # Do radiative transfer
            tau_hor[j] = McRT.radiative_transfer(nu_obs, my_ray, alpha_nu)
        for j, v in enumerate(v_arr):
            for i, b in enumerate(b_arr):
                if b < sphere_extent:
                    # Add tube contribution:
                    # Fraction of tube is 1-sin(theta) = 1-sqrt(1-(adj/hyp)^2)
                    tau_arr[j][i] += (tau_hor[j]
                                      *(1.-np.sqrt(1.-(b/sphere_extent)**2)))
        # Calculate disk average obscuration of spherical extent + horizontal tube
        h = min(r_star, sphere_extent)
        chi = h/r_star
        zeta = np.sqrt(1.-chi**2)
        eta = min(1., sphere_extent/r_star)
        w_hor = (2./const.pi)*(math.acos(zeta)+chi*zeta) - eta**2
        for j in range(len(v_arr)):
            # Obscuration from spherical extent
            self.Obscuration[j] = eta**2-np.trapz(2.*b_arr*np.exp(-tau_arr[j]),
                                                  x=b_arr)/r_star**2
            # Obscuration from horizontel extent
            self.Obscuration[j] += (1-np.exp(-tau_hor[j]))*w_hor
