import numpy as np
import builtins

from ..radiative_transfer import rt_ray

def pop_yt_ray(i, j, ds, IP):
    """
    Description:
        Returns a populated rt_ray given the ray and data structure
        containing the simulation results. This is tailored towards results
        being analyzed with the yt project.
    
    Arguments:
        i: ith column of pixel
        j: jth row of pixel
        ds: data structure that contains simulation result
        IP: image plane object that contains ray's enter and exit points
        
    Returns:
        A populated rt_ray
    """
    # Grab enter and exit location of ray
    en_x, en_y, en_z = IP.ray_enter[j][i]
    ex_x, ex_y, ex_z = IP.ray_exit[j][i]
    # Turn ray into a vector
    ray_vec = [ex_x-en_x, ex_y-en_y, ex_z-en_z]
    ray_len = np.sqrt((ex_x-en_x)**2+(ex_y-en_y)**2+(ex_z-en_z)**2)
    ray_unt = [vec/ray_len for vec in ray_vec]
    ray_sgn = [np.sign(ray_vec[0]),
               np.sign(ray_vec[1]),
               np.sign(ray_vec[2])]
    # Grab data along ray 
    ds_ray = ds.ray((en_x,en_y,en_z), (ex_x,ex_y,ex_z))
    ray_srt = np.argsort(ds_ray['t'])
    # Generate rt_ray object
    n_cells = len(ds_ray['x'])
    radtran_ray = rt_ray(n_cells)
    # Find the cell faces in the direction of the ray
    cell_face = [[None]*n_cells]*3
    cell_face[0] = ds_ray['x'][ray_srt]+ray_sgn[0]*ds_ray['dx'][ray_srt]/2.
    cell_face[1] = ds_ray['y'][ray_srt]+ray_sgn[1]*ds_ray['dy'][ray_srt]/2.
    cell_face[2] = ds_ray['z'][ray_srt]+ray_sgn[2]*ds_ray['dz'][ray_srt]/2.
    # We are going to walk along the ray, so set out entry position and
    # calculate an upper bound on the maximum ds: ds=np.sqrt(3)*max(dx)
    pos = [en_x, en_y, en_z]
    over_max_ds = 2.*builtins.max(ds_ray['dx'][0],
                                  builtins.max(ds_ray['dy'][0],
                                               ds_ray['dz'][0]))
    # Populate values of ray
    for t in range(n_cells):
        # Calculate the line of sight velocity
        radtran_ray.v_LOS[t] = (ds_ray['velocity_x'][ray_srt][t].v*ray_unt[0]
                                +ds_ray['velocity_y'][ray_srt][t].v*ray_unt[1]
                                +ds_ray['velocity_z'][ray_srt][t].v*ray_unt[2])
        # Calculate the ray's neutral number density
        radtran_ray.n_abs[t] = ds_ray['nH'][ray_srt][t].v
        # Calculate the ray's ds of the cell
        ## To do so we step along the ray, and measure distance between hitting
        ## the next cell face. The distance is the smallest ds to hit cell_face
        ## in each dimension.
        cell_ds = over_max_ds
        for k in range(3):
            if ray_unt[k] != 0:
                cell_ds = builtins.min(cell_ds,
                                       (cell_face[k][t].v-pos[k])/ray_unt[k])
        radtran_ray.ds[t] = cell_ds
        # Calculate the ray's Temperature
        radtran_ray.T[t] = ds_ray['T'][ray_srt][t].v
        # Calculate the ray's mean molecular weight
        radtran_ray.mu[t] = ds_ray['mu'][ray_srt][t].v
        # Advance our position to next cell face
        for k in range(3):
            pos[k] = pos[k]+cell_ds*ray_unt[k]
    return radtran_ray