from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_bcs_base():
    Rmin,Rmax,rho_Rmin, T_Rmin = sim.base_bcs()
    assert np.floor(Rmin) == 1.0 #Rmin
    assert Rmax > Rmin
    assert rho_Rmin/1e5 <= 1.0 #rho(Rmin) in units of 1e-15 g/cm3
    assert (T_Rmin <= 1.0) & (T_Rmin > 0.1) #T(Rmin) in units of 1e4 K

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

