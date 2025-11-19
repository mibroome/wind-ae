from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
import importlib.resources as pkg_resources
import numpy as np

sim = wind_sim()
path = str(pkg_resources.files('wind_ae'))+'/'
filename = path+'/saves/test.csv'
sim.load_planet(filename)

#loading current windsoln. Running "itself"
def test_polish_bcs():
    assert sim.polish_bcs() == 0

def test_converge_mol_atomic_transition():
    assert sim.converge_mol_atomic_transition() == 0
    sim.windsoln.bolo_heat_cool = 0
    assert sim.converge_mol_atomic_transition() is None
    assert sim.converge_mol_atomic_transition(polish=True) == 0
    sim.windsoln.bolo_heat_cool = 1
    
def test_Rmax_convergence():
    sim.run_wind(expedite=True) #turns off outward integration
    assert sim.converge_Rmax() == 0

def test_Ncol_convergence():
    sim.run_wind(expedite=True)
    assert np.round(sim.self_consistent_Ncol()[1][1],5) == 0.00029
    sim.integrate_out()
    assert np.round(sim.self_consistent_Ncol()[1][1],4) == 0.0051

def test_turn_off_bolo():
    sim.turn_off_bolo()
    assert sim.windsoln.flags_tuple[-2] == 0 