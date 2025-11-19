from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import numpy as np
import pytest

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'

def test_add_uservars():
    sim.load_planet(filename,calc_postfacto=False)
    sim.windsoln.add_user_vars()
    assert sim.windsoln.R_cori is not None

sim.load_planet(filename,calc_postfacto=False)
def test_calc_metallicity():
    assert sim.windsoln.current_metallicity() == 1.0

def test_calc_mdot():
    sim.windsoln.calc_massloss()
    assert np.round(sim.windsoln.Mdot) == 20135455210

def test_tau_array():
    assert np.round(sim.windsoln.tau_array(20)[0], 1) == 1628.8

def test_calc_tau1_radius():
    assert np.round(sim.windsoln.calc_tau1_radius(20)[0], 1) == 1.4

def test_calc_tau1_radius_error():
    with pytest.raises(Exception) as e_info:
        sim.windsoln.calc_tau1_radius(0)[0]

def test_calc_R_exo():
    assert np.round(sim.windsoln.calc_R_exo()[0], 1) == 3.8

def test_calc_Jeans():
    sim.windsoln.calc_R_exo()
    assert np.round(sim.windsoln.calc_Jeans(),1) == 391382772.9

def test_calc_Coriolis():
    sim.windsoln.calc_Coriolis()
    assert np.round(sim.windsoln.R_cori,1) == 6.4

def test_calc_ballistic():
    sim.windsoln.add_user_vars()
    sim.windsoln.calc_ballistic()
    assert sim.windsoln.n_ball == 31