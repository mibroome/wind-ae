from wrapper.relax_wrapper import wind_simulation as wind_sim
from wrapper.wrapper_utils import constants as const
from wrapper.wrapper_utils.system import system
import numpy as np
import time
import os

mass_grid = np.load('starts/mass_grid.npy')
radius_array = np.load('starts/radius_array.npy')

node=int(os.getcwd().split('/')[-1][4:])
# node=11

terminal_fail = 0
savepath='../relaxed-wind_good/saves/Grids/grid_9-17/'
progpath=''

sim = wind_sim()
open(progpath+f'progress{node}.txt', 'w').close()

keep_going=True
r_idx=0
while keep_going == True:
    try:
        sim.load_planet(f'starts/hi_{mass_grid[node][r_idx]:.2f}Me_{radius_array[r_idx]:.2f}Re.csv',
                       print_atmo_composition=False)
        keep_going=False
    except FileNotFoundError:
        r_idx+=1

for i in range(1,len(radius_array)):
    start = time.time()
    Mp = mass_grid[node][i]
    Rp = radius_array[i]
    
    local = time.localtime(start)
    f = open(progpath+f'progress{node}.txt','a+')
    f.write(f'{i} | {Mp:.2f}, {Rp:.2f}, Start time: {time.asctime(local)}\n')
    f.close()
    
    #ramping goes here    
    planet_tuple = np.copy(sim.windsoln.planet_tuple)
    planet_tuple[:2] = np.array([Mp*const.Mearth,Rp*const.Rearth])
    result = sim.ramp_grav(system(*planet_tuple),integrate_out=False,make_plot=False,
                           converge=False)
    if result == 0:
        sim.polish_bcs()
        actual_Mp = sim.windsoln.Mp/const.Mearth
        actual_Rp = sim.windsoln.Rp/const.Rearth
        sim.save_planet(savepath+f'hi_{actual_Mp:.2f}Me_{actual_Rp:.2f}Re.csv',overwrite=True)
        sim.save_planet(f'saves/last_good.csv',overwrite=True) #save locally
        fail_result=''
        terminal_fail=0 #reset the counter in case that was just a bad point
    else:
        terminal_fail+=1
        fail_result='FAILED '+str(result)
        sim.load_planet('saves/last_good.csv') #revert to last stable soln
    
    end = time.time()
    total = (end-start)
    
    f = open(progpath+f'progress{node}.txt','a+')
    f.write(f'       {fail_result}          Total time: '+time.strftime('%H:%M:%S', time.gmtime(total))+'\n')
    f.close()
    
    #if it has failed on two solutions in a row
    if terminal_fail>=2:
        f = open(progpath+f'progress{node}.txt','a+')
        f.write(f'--Stopping--')
        f.close()
        break 
