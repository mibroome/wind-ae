from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

#loading current windsoln. Running "itself"
def test_run_wind():
    assert sim.run_wind() == 0

def test_run_wind_expedite():
    sim.run_wind(expedite=True)
    assert sim.windsoln.flags_tuple[-1] == 0
