from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
from wind_ae.wrapper.wrapper_utils.system import system
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_physics_tuple():
    assert np.round(sim.physics.physics_tuple()[0][0],2) == 0.80

def test_physics_value():
    assert np.round(sim.physics.value()[0],2) == 0.80

def test_physics_assign():
    sim.physics.assign("HX", [0.9,0.1])
    assert sim.physics.HX[0] == 0.9