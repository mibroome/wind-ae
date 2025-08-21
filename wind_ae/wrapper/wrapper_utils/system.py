#!/bin/env python3

import copy

from . import constants as const

class system:
    """
    Object containing information about the planetary system
    """
    def __init__(self, Mp, Rp, Mstar, semimajor, Ftot, Lstar, name='My Planet'):
        """
        Arguments:
            Mp: Planetary mass (g)
            Rp: Transit oberservational radius (cm)
            Mstar: Stellar mass (g)
            semimajor: orbital seperation (cm)
            Ftot: Total integrated ionizing flux at planet (erg/s/cm^2)
        """
        # system parameters
        self.Mp = Mp
        self.Rp = Rp
        self.Mstar = Mstar
        self.semimajor = semimajor
        self.Ftot = Ftot
        self.Lstar = Lstar
        self.name = name


    def system_tuple(self):
        """
        Return variable values that describe the planetary system. Currently not
        including the given (unphysical) name of system in return.
        """
        return (self.Mp, self.Rp, self.Mstar, self.semimajor, self.Ftot, self.Lstar)


    def print_system(self, norm='Jupiter'):
        """Prints system variables in table.

        Args:
            norm (str): The normalization scheme to use ('Jupiter' or 'Earth').
        """
        if norm == 'Jupiter':
            mass_norm = const.Mjupiter
            rad_norm = const.Rjupiter
            letter = 'J'
        elif norm == 'Earth':
            mass_norm = const.Mearth
            rad_norm = const.Rearth
            letter = 'E'
        print(f'{self.name}:\n'
              '  System parameters (cgs)               # Normalized units\n'
              '    Mp:        {:e} g           # %8.2f M%s\n'
              '    Rp:        {:e} cm          # %8.2f R%s\n'
              '    Mstar:     {:e} g           # %8.2f Msun\n'
              '    semimajor: {:e} cm          # %8.2f au\n'
              '    Ftot:      {:e} erg/cm^2/s  # %8.2e FuvEarth\n'
              '    Lstar:     {:e} erg/s       # %8.2e Lsun\n'
              .format(self.Mp, self.Rp, self.Mstar, self.semimajor, self.Ftot, self.Lstar)
              % (self.Mp/mass_norm, letter, self.Rp/rad_norm, letter,
                 self.Mstar/const.Msun, self.semimajor/const.au, self.Ftot/const.FuvEarth,
                 self.Lstar/const.Lsun))
        return


    def value(self, var=None):
        if var is None:
            return (self.Mp, self.Rp, self.Mstar, self.semimajor, self.Ftot)
        elif var == "Mp":
            return self.Mp
        elif var == "Rp":
            return self.Rp
        elif var == "Mstar":
            return self.Mstar
        elif var == "semimajor":
            return self.semimajor
        elif var == "Ftot":
            return self.Ftot
        elif var == "Lstar":
            return self.Lstar
        elif var == "name":
            return self.name
        else:
            print("Don't recoginze var: %s" % var)
            return


    def assign(self, var, value):
        """
        Description:
            Assigns a new value to the given variable. Recognized variables:
            Mp, Rp, Mstar, semimajor, Ftot, name

        Arguments:
            var: name of variable to set
            value: numerical value to set variable to (in cgs)
        """
        if var == "Mp":
            self.Mp = value
        elif var == "Rp":
            self.Rp = value
        elif var == "Mstar":
            self.Mstar = value
        elif var == "semimajor":
            self.semimajor = value
        elif var == "Ftot":
            self.Ftot = value
        elif var == "Lstar":
            self.Lstar = value
        elif var == "name":
            self.name = value
        else:
            print("Don't recoginze var: %s" % var)
        return


    def deepcopy(self):
        new = self.__class__(copy.deepcopy(self.Mp),
                             copy.deepcopy(self.Rp),
                             copy.deepcopy(self.Mstar),
                             copy.deepcopy(self.semimajor),
                             copy.deepcopy(self.Ftot),
                             copy.deepcopy(self.Lstar),
                             copy.deepcopy(self.name))
        return new
