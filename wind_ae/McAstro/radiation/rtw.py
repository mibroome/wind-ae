#!/usr/bin/env python

from .radiative_transfer import radiative_transfer


class athena_ds:
    def __init__(self, filename):
        return


def rtask(n_cpus, rank, orb_elem, ds, alpha_func, r_star, I0_star, dx1, dx2,
         nu_lo, nu_hi, n_nu, ph_lo, ph_hi, n_ph,
         ix1_lo, ix1_hi, n_ix1, ix2_lo, ix2_hi, n_ix2):
    
    n_tot = n_oe*n_nu*n_ix1*n_ix2
    
    phs = np.linspace(ph_lo, ph_hi, n_ph)
    spec = np.linspace(nu_lo, nu_hi, n_nu)
    ix1 = np.linspace(ix1_lo, ix1_hi, n_ix1)
    ix2 = np.linspace(ix2_lo, ix2_hi, n_ix2)
    
    r_task = n_tot%n_cpus
    n_task = (n_tot-r_task)//n_cpus
    if rank < r_task:
        n_task += 1
        i_min = rank*n_task
        i_max = i_min+n_task
    else:
        i_min = r_task+rank*n_task
        i_max = i_min+n_task
        
    x1, x2, tau, I = (np.empty(n_task) for i in range(4))
    last_ix1, last_ix2 = None, None
    for i in range(i_min, i_max):
        i_nu = i%n_nu
        i_ph = ((i-i_nu)//n_nu)%n_oe
        i_ix1 = ((i-i_nu-i_ph)//(n_nu*n_ph))%n_ix1
        i_ix2 = ((i-i_nu-i_ph-i_ix1)//(n_nu*n_ph*n_ix1))%n_ix2
        
        orb_elem.TA = phs[i_ph]
        if ix1[i_ix1] != last_ix1 or ix2[i_ix2] != last_ix2:
            ray = fetch_ray(orb_elem, ds, ix1, ix2, dx1, dx2, r_star, I0_star)
            last_ix1 = ix1[i_ix1]
            last_ix2 = ix2[i_ix2]
        
        tau[i-i_min] = radiative_transfer(spec[i_nu], ray, alpha_func)
        I[i-i_min] = ray.I[-1]
        x1[i-i_min] = last_ix1*dx1
        x2[i-i_min] = last_ix2*dx2
        
    if savefile is not None:
        np.savez_compressed(savefile, x1=x1, x2=x2, tau=tau, I=I)
        
        
def fetch_ray(orb_elem, ds, ix1, ix2, dx1, dx2, r_star, I0_star):
    """
    Agruments:
        orb_elem: orbital elements of planet
        x_mn: left boundaries centered on planet (in centimeters)
        x_mx: right boundaries centered on planet (in centimeters)

    Notes:
        ^We use the decimal package to obtain maintain accuracy
        ^primary ray check not necessary
        ^could use more comments to explain method

    Source:
        1. https://en.wikipedia.org/wiki/Orbital_elements
        2. Section 2.8 of Solar System Dynamics by Murray & Dermott
    """
    # Convert inputs to decimal types to maintain accuracy
    adist = dm.Decimal(orb_elem.adist)
    eccen = dm.Decimal(orb_elem.eccen)
    Inc = dm.Decimal(orb_elem.Inc)
    LoAN = dm.Decimal(orb_elem.LoAN)
    AoP = dm.Decimal(orb_elem.AoP)
    TA = dm.Decimal(orb_elem.TA)
    dist = dm.Decimal(orb_elem.dist)
    dx1 = dm.Decimal(dx1)
    dx2 = dm.Decimal(dx2)
    # cosine and sine of orbital angular elements
    if LoAN != _dm0:
        cLoAN = dm.Decimal(math.cos(LoAN))
        sLoAN = dm.Decimal(math.sin(LoAN))
    else:
        cLoAN = _dm1
        sLoAN = _dm0
    if AoP != _dm0:
        cAoP = dm.Decimal(math.cos(AoP))
        sAoP = dm.Decimal(math.sin(AoP))
    else:
        cAoP = _dm1
        sAoP = _dm0
    if TA != _dm0:
        cTA = dm.Decimal(math.cos(TA))
        sTA = dm.Decimal(math.sin(TA))
    else:
        cTA = _dm1
        sTA = _dm0
    if Inc != _dm0:
        cInc = dm.Decimal(math.cos(Inc))
        sInc = dm.Decimal(math.sin(Inc))
    else:
        cInc = _dm1
        sInc = _dm0
    # Cross term in rotation matrix
    cAoPTA = cAoP*cTA - sAoP*sTA
    sAoPTA = cAoP*sTA + sAoP*cTA
    # Calculate distance needed to shift from planet to stellar frame
    frame_shift = [_dm0]*3
    frame_shift[0] = adist*(_dm1-eccen**2)/(_dm1+eccen*cTA)
    # Shift simulation boundaries into stellar frame
    for i in range(3):
        x_mx[i] = dm.Decimal(x_mx[i])+frame_shift[i]
        x_mn[i] = dm.Decimal(x_mn[i])+frame_shift[i]
    # unit vector in the star's reference frame that points to the observer
    ell = [0]*3
    ell[0] = cLoAN*cAoPTA - sLoAN*sAoPTA*cInc
    ell[1] = -cLoAN*sAoPTA - sLoAN*cAoPTA*cInc
    ell[2] = sLoAN*sInc

    # Allow for a small amount of tolerance when bound seeking
    abs_tol = dm.Decimal(1.0e-8)
    tol_mn = [_dm1-dm.Decimal(np.sign(x))*abs_tol for x in x_mn]
    tol_mx = [_dm1+dm.Decimal(np.sign(x))*abs_tol for x in x_mx]

    # Useful normalization factor: ell[0]**2 + ell[1]**2
    N1 = dm.Decimal(math.sqrt(_dm1-ell[2]**2))
    # parametric value of planet's location
    t_p = _dm1-(adist/dist)
    # Perform the same bound seeking method for where ray's enter and exit
    # the domain. However, now we are interested in secondary ray's that
    # create the image plane rather than the primary ray. These rays come
    # from the lab and interest locations (i*dx1, j*dx2) from the primary ray
    # at the location of the planet.
    ell_p = [_dm0]*3
    s_en, s_ex = (None for k in range(2))
    val = dist*t_p + j*dx2*ell[2]/N1
    ell_p[0] = dm.Decimal(-val*ell[0]-ix1*dx1*ell[1]/N1)
    ell_p[1] = dm.Decimal(-val*ell[1]+ix1*dx1*ell[0]/N1)
    ell_p[2] = dm.Decimal(-val*ell[2]+ix2*dx2/N1)
    for k in range(3):
        if ell_p[k] != 0:
            val2_mx = (x_mx[k]-dist*ell[k])/ell_p[k]
            if (tol_mn[0]*x_mn[0]
                <= val2_mx*ell_p[0]+dist*ell[0]
                <= x_mx[0]*tol_mx[0]
                and 
                tol_mn[1]*x_mn[1]
                <= val2_mx*ell_p[1]+dist*ell[1]
                <= x_mx[1]*tol_mx[1]
                and
                tol_mn[2]*x_mn[2]
                <= val2_mx*ell_p[2]+dist*ell[2]
                <= x_mx[2]*tol_mx[2]
               ):
                if s_en is None:
                    s_en = (x_mx[k]-dist*ell[k])/ell_p[k]
                else:
                    s_ex = (x_mx[k]-dist*ell[k])/ell_p[k]
            val2_mn = (x_mn[k]-dist*ell[k])/ell_p[k]
            if (tol_mn[0]*x_mn[0]
                <= val2_mn*ell_p[0]+dist*ell[0]
                <= x_mx[0]*tol_mx[0]
                and
                tol_mx[1]*x_mn[1]
                <= val2_mn*ell_p[1]+dist*ell[1]
                <= x_mx[1]*tol_mx[1]
                and
                tol_mx[2]*x_mn[2]
                <= val2_mn*ell_p[2]+dist*ell[2]
                <= x_mx[2]*tol_mx[2]
               ):
                if s_en is None:
                    s_en = (x_mn[k]-dist*ell[k])/ell_p[k]
                else:
                    s_ex = (x_mn[k]-dist*ell[k])/ell_p[k]
    if s_en is None or s_ex is None:
        print(s_en, s_ex)
        message = ('Chose a smaller image plane, as ray ({:d},{:d})'
                   ' did not pass through the simulation domain.\n'
                   .format(ix1,ix2))
        sys.stderr.write(message)
        sys.exit()
    if (s_en > s_ex):
        s_en, s_ex = s_ex, s_en
    # Store the (x,y,z) tuple of the ray enter and exit location
    en_x, en_y, en_z = tuple(map(float, (
        (s_en*ell_p[0]+dist*ell[0]-frame_shift[0]),
        (s_en*ell_p[1]+dist*ell[1]-frame_shift[1]),
        (s_en*ell_p[2]+dist*ell[2]-frame_shift[2]))))
    ex_x, ex_y, ex_z = tuple(map(float, (
        (s_ex*ell_p[0]+dist*ell[0]-frame_shift[0]),
        (s_ex*ell_p[1]+dist*ell[1]-frame_shift[1]),
        (s_ex*ell_p[2]+dist*ell[2]-frame_shift[2]))))
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
            
    # Set ray's intial flux
    if np.sqrt((ix2*dx2)**2+(ix1*dx1)**2) <= r_star:
        radtran_ray.I[0] = I0_star
    else:
        radtran_ray.I[0] = 0.0
        
    return radtran_ray