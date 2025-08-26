from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_bcs_base():
    assert sim.base_bcs() == (np.float64(1.0566958974807499), 6.381446068152, np.float64(21357.70226090895), 0.15348209185160616)

def test_ramp_bcs():
    assert sim.ramp_base_bcs() == 0

def test_ramp_Rmin():
    x = 1.0566958974807499
    assert sim.ramp_Rmin(x) == 0
    x = 40
    assert sim.ramp_Rmin(x) is None

def test_ramp_rho_rmin():
    assert sim.ramp_rho_rmin(21357.70226090895) == 0
    assert sim.ramp_rho_rmin(1) is None

def test_ramp_T_rmin():
    assert sim.ramp_T_rmin(0.15348209185160616) == 0
    assert sim.ramp_T_rmin(2) is None

