import numpy as np
from scipy import integrate

from wind_ae.McAstro.utils import constants as const

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

def ballistic_jacobian(t, p, Mp, Mstar, smjr, Omg):
    x, y, z, u, v, w = p
    
    x2 = x**2
    xa2 = (x+smjr)**2
    y2 = y**2
    z2 = z**2
    
    r_planet_2 = x2+y**2+z**2
    r_planet_25 = r_planet_2**(2.5)
    r_star_2 = (x+smjr)**2+y**2+z**2
    r_star_25 = r_star_2**(2.5)
    
    ja_kobi_x = np.array([0,0,0,1,0,0])
    ja_kobi_y = np.array([0,0,0,0,1,0])
    ja_kobi_z = np.array([0,0,0,0,0,1])
    
    ja_kobi_u = np.array([
        (-const.G*(Mp*(r_planet_2-3*x2)/r_planet_25
                   +Mstar*(r_star_2-3*xa2)/r_star_25)+Omg**2),
        (-const.G*(Mp*(-3*y*x)/r_planet_25
                     +Mstar*(-3*y*(x+smjr))/r_star_25)),
        (-const.G*(Mp*(-3*z*x)/r_planet_25
                     +Mstar*(-3*z*(x+smjr))/r_star_25)),
        0, 2*Omg, 0])
    ja_kobi_v = np.array([
        (-const.G*(Mp*(-3*x*y)/r_planet_25
                   +Mstar*(-3*y*(x+smjr)/r_star_25))),
        (-const.G*(Mp*(r_planet_2-3*y2)/r_planet_25
                     +Mstar*(r_star_2-3*y2)/r_star_25)+Omg**2),
        (-const.G*(Mp*(-3*z*y)/r_planet_25
                     +Mstar*(-3*z*y)/r_star_25)),
        -2*Omg, 0, 0])
    ja_kobi_w = np.array([
        (-const.G*(Mp*(-3*x*z)/r_planet_25
                   +Mstar*(-3*z*(x+smjr))/r_star_25)),
        (-const.G*(Mp*(-3*y*z)/r_planet_25
                     +Mstar*(-3*y*z)/r_star_25)),
        (-const.G*(Mp*(r_planet_2-3*z2)/r_planet_25
                     +Mstar*(r_star_2-3*z2)/r_star_25)),
        0, 0, 0])
        
    ja_kobi = np.array([
        ja_kobi_x, ja_kobi_y, ja_kobi_z,
        ja_kobi_u, ja_kobi_v, ja_kobi_w
    ])
    
    return ja_kobi
