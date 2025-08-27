from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
from wind_ae.wrapper.wrapper_utils.system import system
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_system_tuple():
    assert sim.system.system_tuple()[0] == 1.33e+30

def test_system_value():
    assert sim.system.value()[0] == 1.33e+30

def test_system_assign():
    sim.system.assign("Mp", 1e30)
    assert sim.system.Mp == 1e30