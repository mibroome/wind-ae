# Common python packages
import os
import re
from scipy import integrate
import matplotlib
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
import warnings
from matplotlib.lines import Line2D
from cycler import cycler
import pandas as pd
import os
import sys
sys.path.append(sys.path[0] + "/..")

cc = plt.rcParams['axes.prop_cycle'].by_key()['color']

params = {
'xtick.direction': 'in',
'ytick.direction': 'in',
}
plt.rcParams.update(params)
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':15})
plt.rc('text', usetex=True)
warnings.filterwarnings("ignore")

# Check that pwd is relaxation directory and code is compiled
if os.path.exists(os.getcwd()+'/Quickstart.ipynb'):
    os.chdir('../') # May need to be set by hand if nonconventional setup
if not os.path.exists(os.getcwd()+'/bin/relaxed_ae'):
    print('ERROR: relaxation code not found!\n'
          '       Either not compiled, or not in relaxation code directory')
    assert(False)

# Wrapper modules
from wrapper.relax_wrapper import wind_simulation as wind_sim
from wrapper.wrapper_utils import constants as const
from wrapper.wrapper_utils.system import system
from wrapper.wrapper_utils.physics import physics
from wrapper.wrapper_utils.spectrum import spectrum
from wrapper.wrapper_utils.plots import *
import McAstro.atoms.atomic_species as McAtom