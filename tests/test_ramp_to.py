from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
from wind_ae.wrapper.wrapper_utils.system import system
import importlib.resources as pkg_resources

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_ramp_to():
    planet = sim.windsoln.planet_tuple
    result = sim.ramp_to(system(*planet),integrate_out=False)
    assert result == 0

def test_ramp_Ftot():
    planet = sim.windsoln.planet_tuple
    result = sim.ramp_Ftot(sim.windsoln.Ftot,integrate_out=False)
    assert result == 0

def test_ramp_var_error():
    assert sim.ramp_var("Ftot",var_class="test",var_end=100) == -100

def test_ramp_grav():
    planet = sim.windsoln.planet_tuple
    result = sim.ramp_grav(system(*planet),integrate_out=False)
    assert result == 0