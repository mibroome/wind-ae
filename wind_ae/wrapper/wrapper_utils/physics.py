#!/bin/env python3

import copy
import pandas as pd
import numpy as np
import importlib.resources as pkg_resources
import wind_ae.wrapper.wrapper_utils.constants as const

class physics:
    """
    Object containing information about the physical parameters of the
    simualation.
    """
    def __init__(self, HX, species_list, molec_adjust, atomic_masses=np.array([0])):
        # phys_params.inp
        filepath = pkg_resources.files('wind_ae.McAstro.atoms').joinpath('atomic_table.txt')
        atomictable = pd.read_csv(filepath,comment='#')
        
        self.HX = HX
        self.species_list = species_list
        if atomic_masses[0] != 0:
            self.atomic_masses = atomic_masses
        else:
            verner_atomic_masses = np.zeros_like(HX)
            for j in range(len(species_list)):
                verner_atomic_masses[j] = atomictable['Atomic Mass'][atomictable['Name']==(species_list[j]).replace(' ','')] #Verner masses
            self.atomic_masses = verner_atomic_masses
            self.atomic_masses[0] = 1.67330000000000013e-24
            
        Z = np.zeros_like(HX)
        Ne = np.zeros_like(HX)
        for j in range(len(species_list)):
            Z[j] = atomictable['Z'][atomictable['Name']==(species_list[j]).replace(' ','')].iloc[0]
            Ne[j] = atomictable['N_electrons'][atomictable['Name']==(species_list[j]).replace(' ','')].iloc[0]
        self.Z = Z
        self.Ne = Ne 
        self.molec_adjust = molec_adjust
        
    def physics_tuple(self):
        return (self.HX, self.species_list, self.atomic_masses)
    
    def print_physics(self):
        print('Physics parameters\n'
              'mass_fraction:        '+str(self.HX)+'    # Mass fraction\n'
              'species_list:         '+self.species_list+'    #\n'
              'atomic_mass:          '+str(self.atomic_masses)+'   #g\n'
              'Z:                    '+str(self.Z)+'         #Atomic number\n'
              'Ne:                   '+str(self.Ne)+'        #Number of electrons per species\n',
              'Molec_adjust:         '+str(self.molec_adjust)+'        #Weighting factor to account for mols below windbase\n')
        
        return
    
    def value(self, var=None):
        if var is None:
            return (self.HX)
        elif var == "HX" or "mass_fraction":
            return self.HX
        elif var == "species_list":
            return self.species_list
        elif var == "atomic_masses":
            return self.atomic_masses
        elif var.upper() == "Z":
            return self.Z
        elif var.lower() == "ne":
            return self.Ne
        elif var == "molec_adjust":
            return self.molec_adjust
        else:
            print("Don't recoginze var: %s" % var)
            return

        
    def assign(self, var, value):
        if var == "HX" or "mass_fraction":
            self.HX = value
        elif var == "species_list":
            self.species_list = value
        elif var == "atomic_masses":
            self.atomic_masses = value
        elif var.upper() == "Z":
            self.Z = value
        elif var.lower() == "ne":
            self.Ne = value
        elif var == "molec_adjust":
            self.molec_adjust = value
        else:
            print("Don't recoginze var: %s" % var)
        return


    def deepcopy(self):
        new = self.__class__(copy.deepcopy(self.HX), 
                            copy.deepcopy(self.species_list),
                            copy.deepcopy(self.atomic_masses),
                            copy.deepcopy(self.Z),
                            copy.deepcopy(self.Ne),
                            copy.deepcopy(self.molec_adjust))
        return new
