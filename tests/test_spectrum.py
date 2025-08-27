from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

def test_load_spectrum():
    sim.load_spectrum()
    assert sim.spectrum is not None

def test_ramp_spectrum():
    sim.load_planet(path+'/saves/test.csv')
    result = sim.ramp_spectrum(Fnorm=0,norm_spec_range=[13.6,2000],goal_spec_range=[13.6,2000])
    assert result == 0

def test_ramp_to_user_spectrum():
    sim.load_planet(path+'/saves/test.csv')
    assert sim.ramp_to_user_spectrum('hd189733', updated_F=sim.windsoln.Ftot,plot=False) == 0