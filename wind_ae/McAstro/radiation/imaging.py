import sys
import math
import builtins
import numpy as np
import decimal as dm

from .radiative_transfer import rt_ray

_dm0 = dm.Decimal(0.0)
_dm1 = dm.Decimal(1.0)

class orbital_elements:
    def __init__(self, adist, eccen, Inc, LoAN, AoP, TA, dist):
        """
        Description:
            Class to contain all the orbital elements that define the
            planet. Also stores the distance to the star system.
            
        Arguments:
            adist: semimajor axis of orbit (in centimeters)
            eccen: eccentricity of orbit
            Inc: inclination of orbital plane (in radians)
            LoAN: orbital plane's longitude of ascending node (in radians)
            AoP: orbit's argument of periapsis (in radians)
            TA: true anomlay of p]lanet (in radians)
            dist: distance between lab and star (in centimeters)
        """
        self.adist = adist
        self.eccen = eccen
        self.Inc = Inc
        self.LoAN = LoAN
        self.AoP = AoP
        self.TA = TA
        self.dist = dist
        

class image_plane:
    def __init__(self, ix1, ix2, dx1, dx2, r_star, I0_star):
        """
        Description:
            Class that contains all information having to do with the image
            plane used for ray-tracing. Also stores the return values of
            doing ray tracing, i.e., the optical depth (tau) and specific
            spectral irradiance (I). Is capable of figuring out which sector
            of the image plane the current cpu is responsible for.
        
        Arguments:
            dx1: physical x distance at adist between pixels (in centimeters)
            dx2: physical y distance at adist between pixels (in centimeters)
            ilow: leftmost pixel, planet at 0
            ihigh: rightmost pixel, planet at 0
            jlow: bottommost pixel, planet at 0
            jhigh: topmost pixel, planet at 0
            
        Keyword arguments:
            rank: cpu rank to know what decomposed grid to do 
            size: number of cpus/grids to decompose image plane over
        """
        self.ix1 = ix1
        self.ix2 = ix2
        self.dx1 = dx1
        self.dx2 = dx2
        self.r_star = r_star
        self.I0_star = I0_star

    
    def generate_ray(self, ds, orb_elem, i, j):
        """
        Agruments:
            orb_elem: orbital elements of planet
            
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
        dx1 = dm.Decimal(self.dx1)
        dx2 = dm.Decimal(self.dx2)
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
        x_mn, x_mx = ([_dm0]*3 for k in range(2))
        for k in range(3):
            x_mn[k] = dm.Decimal(float(ds.domain_left_edge[k].v))+frame_shift[k]
            x_mx[k] = dm.Decimal(float(ds.domain_right_edge[k].v))+frame_shift[k]
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
        ell_p[0] = dm.Decimal(-val*ell[0]-i*dx1*ell[1]/N1)
        ell_p[1] = dm.Decimal(-val*ell[1]+i*dx1*ell[0]/N1)
        ell_p[2] = dm.Decimal(-val*ell[2]+j*dx2/N1)
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
                       .format(i,j))
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
            radtran_ray.n_abs[t] = ds_ray['n_HI'][ray_srt][t].v
            # Calculate the ray's ds of the cell
            ## To do so we step along the ray, and measure distance between hitting
            ## the next cell face. The distance is the smallest ds to hit cell_face
            ## in each dimension.
            cell_ds = float(over_max_ds.v)
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