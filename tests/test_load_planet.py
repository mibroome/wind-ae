from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import os

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'saves/test.csv'

if os.path.exists(path+'saves/test_save.csv'):
    os.remove(path+'saves/test_save.csv')
if os.path.exists(path+'saves/test_easy_output.csv'):
    os.remove(path+'saves/test_easy_output.csv')

def test_load_planet():
    sim.load_planet(filename)
    assert sim.guess is not None

def test_load_uservars():
    sim.load_uservars(filename)
    assert sim.windsoln.R_cori is not None

def test_generate_rate_coeffs():
    sim.load_uservars(filename)
    sim.generate_rate_coeffs()
    with open(path+'src/rate_coeffs.h', 'r') as f:
        lines = f.readlines()
        output = lines[3][:15]
    assert output == "double R[6][59]"

sim.load_planet(filename)
def test_save_planet():
    sim.save_planet(path+'saves/test_save.csv')
    assert os.path.exists(path+'saves/test_save.csv') 

def test_easy_output():
    sim.easy_output_file(output_file=path+'saves/test_easy_output.csv')
    assert os.path.exists(path+'saves/test_easy_output.csv')

if os.path.exists(path+'saves/test_save.csv'):
    os.remove(path+'saves/test_save.csv')
if os.path.exists(path+'saves/test_easy_output.csv'):
    os.remove(path+'saves/test_easy_output.csv')