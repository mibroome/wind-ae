from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_turn_off_tidal_grav():
    sim.turn_off_tidal_grav()
    assert sim.windsoln.flags_tuple[1] == 0

def test_raise_Ncol_sp():
    sim.raise_Ncol_sp(by_factor=2)
    assert np.round(sim.windsoln.Ncol_sp[0],3) == 0.041

def test_integrate_out():
    sim.run_wind(expedite=True)
    sim.integrate_out()
    assert sim.windsoln.flags_tuple[-1] == 1

def test_erf_velocity():
    sim.load_planet(filename)
    assert sim.erf_velocity(return_idx=True)[-1] == 103