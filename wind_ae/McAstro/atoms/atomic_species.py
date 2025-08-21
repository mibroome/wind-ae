import sys, os
import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator

from wind_ae.McAstro.utils import constants as const

_Zelem = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7,
              'O':8, 'F':9,'Ne':10, 'Na':11, 'Mg':12, 'Al':13,
              'Si':14, 'S':16, 'Ar':18, 'Ca':20, 'Fe':26}

filed = os.path.dirname(os.path.abspath(__file__))
_Verner = pd.read_csv(filed+'/Verner.csv', comment='#')
_atomic_masses = pd.read_csv(filed+'/atomic_masses.csv', comment='#')

def roman_to_arabic(roman):
    numerals = {'M':1000, 'D':500, 'C':100, 'L':50, 'X':10, 'V':5, 'I':1}
    arabic = 0
    last_value = 0
    for value in (numerals[c] for c in reversed(roman.upper())):
        arabic += (value, -value)[value < last_value]
        last_value = value
    return arabic


def arabic_to_roman(arabic):
    numerals = [(1000, 'M'), (900, 'CM'), (500, 'D'), (400, 'CD'),
                (100, 'C'), (90, 'XC'), (50, 'L'), (40, 'XL'),
                (10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'),
                (1, 'I')]
    roman = ''
    for (n, s) in numerals:
        (q, arabic) = divmod(arabic, n)
        roman += s*q
    return roman


def spectroscopy_to_atomic_notation(name):
    spect = name.split()
    try:
        Z = _Zelem[spect[0]]
    except KeyError:
        print("Error: {:s} is not in Verner's atomic list.".format(spect[0]))
        return -1, -1
    Ne = Z - (roman_to_arabic(spect[1])-1)
    if Ne < 0:
        print('Spectrscopic species does not exist, '
              'max ionized state for {:s} is {:s}.',
              spect[0], arabic_to_roman(Z+1))
        return -1, -1
    return Z, Ne


def formatting_species_list(species_list): #Added 4/1/2021 - MIB
    '''Converts list of species from any format to "Z Roman" format (e.g. "H I")'''
    for s in range(len(species_list)):
        name = species_list[s]
        if str.isupper(name[0]) == False:
            name = str.upper(name[0])+name[1:]
        while " " not in name:
            if (name[:2])[-1] in ['M', 'D', 'C', 'L', 'X', 'V', 'I','0','1','2','3','4','5','6','7','8','9']:
                name = name[0]+' '+name[1:]
            else:
                name = name[:2]+' '+name[2:]
        else:
            element = str((name.split())[0])
            if len(element)>1:
                element = element[0]+element[1].lower() 
            ion_num = (name.split())[1]
            if str(ion_num[-1]) not in ['M', 'D', 'C', 'L', 'X', 'V', 'I']:
                roman = arabic_to_roman(int(ion_num))
                species_list[s] = element+' '+str(roman) 
            else:
                species_list[s] = name
    return species_list


class atomic_species:
    def __init__(self, name, A=None):
        """
        Arguments:
            name: spectroscopic_name
        """
        self.name = name
        self.Z, self.Ne = spectroscopy_to_atomic_notation(name)
        self.verner_data = None
        if self.Ne > 0:
            self.verner_data  = {
                var : val for (var, val) in
                zip(_Verner.columns, _Verner[(_Verner['Z']==self.Z)
                                             &(_Verner['Ne']==self.Ne)]
                    .values[0])
            }
        if A is None:
            if self.Z == 1:
                self.A = 1
            else:
                self.A = 2*self.Z
        self.ion_pot = self.verner_data['E_th']*const.eV
        E_arr = np.linspace(self.verner_data['E_th'], self.verner_data['E_max'],
                            100000)
        wl_arr = np.linspace(const.hc/(self.verner_data['E_max']*const.eV),
                             const.hc/(self.verner_data['E_th']*const.eV),
                             100000)
        self.sigma = Akima1DInterpolator(E_arr, self.cross_section(E_arr ))
        self.sigma_wl = (
            Akima1DInterpolator(wl_arr,
                                self.cross_section(const.hc/(wl_arr*const.eV)))
        )
        self.atomic_data = _atomic_masses.loc[(_atomic_masses['Z']==self.Z)
                                              &(_atomic_masses['A']==self.A)]
        self.mass = self.atomic_data['mass']*const.Da


    def sigma_find_E(self, sigma):
        E_arr = np.linspace(self.verner_data['E_th'], self.verner_data['E_max'],
                            100000)
        return Akima1DInterpolator(E_arr, self.sigma(E_arr)-sigma).roots()


    def cross_section(self, E, units='cm^2', valid_range=True):
        if self.verner_data is None:
            print('No Verner data due to faulty atomic_species.\n'
                  '  Try initializing.')
            return
        kshell_ionpots = { 6:308.67523704,  7:426.57706307,  8:563.23652347, 
                          12:1336.8046517 , 14:1877.59068171, 16:2512.3371317 }
        x = E/self.verner_data['E_0'] - self.verner_data['y_0']
        y = np.sqrt(x**2+self.verner_data['y_1']**2)
        Fy = (((x-1)**2+self.verner_data['y_w']**2)
              *(y**((self.verner_data['P']-11.)/2.)
                /((1+np.sqrt(y/self.verner_data['y_a']))
                  **(self.verner_data['P']))))
        sigma_0 = self.verner_data['sigma_0']
        sigma = Fy*sigma_0
        if (self.Z in [6,7,8,9,12,14,16]): #x-rays photoionize certain species' K shell electrons instead
            #Band et al. 1990
            a1,a2,a3,a4,a5 = [11.132,5.612,6.9668,8.9690,0.35327] #fit params
            b1,b2,b3,b4    = [1.6798e-2,3.6703e-4,9.7881e-4,8.0332e-3]

            lowZs = np.array([[2, 2, 112.88, 3.9492],
            [3, 3, 144.47, 9.3820],
            [3, 2, 201.22, 3.3064],
            [4, 4, 184.11, 13.931],
            [4, 3, 175.67, 15.889],
            [4, 2, 200.76, 10.343]]) #E_0 and sigma_0 for smaller elements

            if (self.Z in lowZs[:,0]) & (self.Ne in lowZs[:,1]):
                E_0,sig_0 = lowZs[(lowZs[:,0]==Z) & (lowZs[:,1]==N)][0][2:]
            else:
                E_0 = a1 + a2*self.Z + a3*self.Z**2 + (a4+a5*self.Z)*self.Ne #in eV
                sig_0 = 1/(b1 + b2*self.Z + b3*self.Z**2 + b4*self.Ne) #in Mb (x1e-18 to get to cm^2)

            y = E / E_0 
            Fy_xray = y**(-1.5) * (1+y**(0.5))**(-4)
            
            high_E = kshell_ionpots[self.Z] #setting Kshell ionization energy
            self.verner_data['E_max'] = np.max(E) #raises the ceiling included in the Verner data 
            if isinstance(sigma,float)==False: 
                sigma[E>=high_E] = Fy_xray[E>=high_E]*sig_0 #above the k shell ionization threshold, use K shell cross-sec
            else: #if monofrequency
                if E>=high_E:
                    sigma = Fy_xray*sig_0
                else:
                    sigma = sig_0
        
        if units == 'cm^2':
            sigma *= 1e-18
        elif units == 'Mb':
            pass
        if valid_range:
#             print(self.verner_data['E_max'])
            return np.where(E <= self.verner_data['E_max'],
                            np.where(E >= self.verner_data['E_th'],
                                    sigma, 0),
                            0)
        else:
            return np.where(E >= self.verner_data['E_th'], sigma, 0)
        
        
    def cross_section_derivative(self, E, units='cm^2', valid_range=True,
                                 wrt='E'):
        sigma = self.cross_section(E, units=units, valid_range=valid_range)
        if np.all(sigma == 0.0):
            return sigma
        x = E/self.verner_data['E_0'] - self.verner_data['y_0']
        y = np.sqrt(x**2+self.verner_data['y_1']**2)
        dy_dx = x/y
        sqrt_y_ya = np.sqrt(y/self.verner_data['y_a'])
        t1 = 2*(x-1)/((x-1)**2+self.verner_data['y_w']**2)
        t2 = (self.verner_data['P']-11.)/(2*y)
        t2 -= self.verner_data['P']/(2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))
        t2 *= dy_dx
        dsigma_dx = sigma*(t1+t2)
        if wrt == 'E':
            return dsigma_dx/self.verner_data['E_0']
        elif wrt == 'lambda':
            return -dsigma_dx*(E*const.eV/const.hc)*(E/self.verner_data['E_0'])
        elif wrt == 'nu':
            return dsigma_dx*(const.h/self.verner_data['E_0'])
        elif wrt == 'x':
            return dsigma_dx
    
    
    def cross_section_second_derivative(self, E, units='cm^2', valid_range=True,
                                        wrt='E'):
        sigma = self.cross_section(E, self.Z, units=units, valid_range=valid_range)
        if np.all(sigma == 0.0):
            return sigma
        dsigma_dx = self.cross_section_derivative(E, units=units,
                                                  valid_range=valid_range,
                                                  wrt='x')
        x = E/self.verner_data['E_0'] - self.verner_data['y_0']
        y = np.sqrt(x**2+self.verner_data['y_1']**2)
        dy_dx = x/y
        d2y_dx2 = (y-x*dy_dx)/y**2
        sqrt_y_ya = np.sqrt(y/self.verner_data['y_a'])

        t1 = 2*(self.verner_data['y_w']**2-(x-1)**2)
        t1 /= ((x-1)**2+self.verner_data['y_w']**2)**2
        t2 = (self.verner_data['P']-11.)/(2*y)
        t2 -= self.verner_data['P']/(2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))
        t2 *= d2y_dx2
        t3 = self.verner_data['P']*(1+2*sqrt_y_ya)
        t3 /= (2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))**2*sqrt_y_ya
        t3 -= (self.verner_data['P']-11.)/(2*y**2) 
        t3 *= dy_dx**2
        # d2sigma_dx2 = dsigma_dx**2/sigma + sigma*(t1+t2+t3))
        d2sigma_dx2 = np.divide(dsigma_dx**2, sigma,
                                out=np.zeros_like(dsigma_dx), where=sigma!=0)
        d2sigma_dx2 += sigma*(t1+t2+t3)
        if wrt == 'E':
            return d2sigma_dx2/self.verner_data['E_0']**2
        elif wrt == 'lambda':
            coeff = (E*const.eV/const.hc)*(E/self.verner_data['E_0'])
            return d2sigma_dx2*coeff**2 + 2*dsigma_dx*coeff*(E*const.eV/const.hc)
        elif wrt == 'nu':
            return d2sigma_dx2*(const.h/self.verner_data['E_0'])**2
        elif wrt == 'x':
            return d2sigma_dx2
        
        
    def cross_section_third_derivative(self, E, units='cm^2', valid_range=True,
                                       wrt='E'):
        sigma = self.cross_section(E, units=units, valid_range=valid_range)
        if np.all(sigma == 0.0):
            return sigma
        dsigma_dx = self.cross_section_derivative(E, units=units,
                                                  valid_range=valid_range,
                                                  wrt='x')
        d2sigma_dx2 = self.cross_section_second_derivative(E, units=units,
                                                    valid_range=valid_range,
                                                    wrt='x')
        x = E/self.verner_data['E_0'] - self.verner_data['y_0']
        y = np.sqrt(x**2+self.verner_data['y_1']**2)
        dy_dx = x/y
        d2y_dx2 = (y-x*dy_dx)/y**2
        d3y_dx3 = (2*x*dy_dx-y*(x*d2y_dx2+2*dy_dx))/y**3
        sqrt_y_ya = np.sqrt(y/self.verner_data['y_a'])
        
        t1 = 4*(x-1)*((x-1)**2-3*self.verner_data['y_w']**2)
        t1 /= ((x-1)**2+self.verner_data['y_w']**2)**3
        t2 = (self.verner_data['P']-11.)/(2*y)
        t2 -= self.verner_data['P']/(2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))
        t2 *= d3y_dx3
        t3 = self.verner_data['P']*(1+2*sqrt_y_ya)
        t3 /= (2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))**2*sqrt_y_ya
        t3 -= (self.verner_data['P']-11.)/(2*y**2) 
        t3 *= 3*dy_dx*d2y_dx2
        t4 = -self.verner_data['P']*(3+8*sqrt_y_ya**2+9*sqrt_y_ya)
        t4 /= (2*self.verner_data['y_a']*sqrt_y_ya*(1+sqrt_y_ya))**3*sqrt_y_ya**2
        t4 += (self.verner_data['P']-11.)/(y**3)
        t4 *= (dy_dx)**3
        
        d3sigma_dx3 = 3*d2sigma_dx2 - 2*dsigma_dx**2/sigma
        d3sigma_dx3 *= dsigma_dx/sigma
        d3sigma_dx3 += sigma*(t1+t2+t3+t4)

        if wrt == 'E':
            return d3sigma_dx3/self.verner_data['E_0']**3
        elif wrt == 'lambda':
            c1 = (E*const.eV/const.hc)
            c2 = (E/self.verner_data['E_0'])
            return -c1**3*c2*(6*dsigma_dx + 6*c2*d2sigma_dx2 + c2**2*d3sigma_dx3)
        elif wrt == 'nu':
            return d3sigma_dx3*(const.h/self.verner_data['E_0'])**3
        elif wrt == 'x':
            return d3sigma_dx3
