import math
import numpy as np
import decimal as dm

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
    def __init__(self, dx1, dx2, ilow, ihigh, jlow, jhigh, rank=0, size=1):
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
        self.dx1 = dx1
        self.dx2 = dx2
        self.ilow = ilow
        self.ihigh = ihigh
        self.jlow = jlow
        self.jhigh = jhigh
        self.ray_enter = None
        self.ray_exit = None
        
        self.x1 = None
        self.x2 = None
        self.tau = None
        self.I = None
        
        self.rank = rank
        self.size = size
        self.rank_decompose()
        
    def decomposition(self, ngrids):
        for i in range(int(np.sqrt(ngrids)), 0, -1):
            if ngrids % i == 0:
                return i, int(ngrids/i)
            
            
    def rank_decompose(self):
        # Calculate number of decompositions in each direction
        nd_x2, nd_x1 = self.decomposition(self.size)
        # Calculate decomposition of a given rank index
        i_x1 = self.rank%nd_x1
        i_x2 = (self.rank-i_x1)//nd_x1
        # Calculate total number of pixels in each direction
        p_x1 = self.ihigh-self.ilow+1
        p_x2 = self.jhigh-self.jlow+1
        # Calculate remainding pixels after equally distributing pixels
        r_x1 = p_x1%nd_x1
        r_x2 = p_x2%nd_x2
        # Calculate base number of pixels done by a decomposition
        self.n_x1 = (p_x1-r_x1)//nd_x1
        self.n_x2 = (p_x2-r_x2)//nd_x2
        # Calculate min/max of decompostions x1 and x2
        if i_x1 < r_x1:
            self.n_x1 += 1
            self.x1_min = self.ilow+i_x1*self.n_x1
            self.x1_max = self.x1_min+self.n_x1
        else:
            self.x1_min = self.ilow+r_x1+i_x1*self.n_x1
            self.x1_max = self.x1_min+self.n_x1
        if i_x2 < r_x2:
            self.n_x2 += 1
            self.x2_min = self.jlow+i_x2*self.n_x2
            self.x2_max = self.x2_min+self.n_x2
        else:
            self.x2_min = self.jlow+r_x2+i_x2*self.n_x2
            self.x2_max = self.x2_min+self.n_x2
        # Set coordinates of each ray in image plane
        self.x1, self.x2 = (np.empty(self.n_x1*self.n_x2) for i in range(2))
        for j in range(self.x2_min, self.x2_max):
            for i in range(self.x1_min, self.x1_max):
                index = (i-self.x1_min)*self.n_x1+(j-self.x2_min)
                self.x1[index] = i*self.dx1
                self.x2[index] = j*self.dx2

    
    def calc_rays(self, orb_elem, x_mn, x_mx):
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
        if len(x_mx) != 3 or len(x_mn) != 3 or np.any(x_mx <= x_mn):
            print("ERROR: boundaries must be 3 dimensional AND \n"
                  "       every x_mx must be greater than the corresponding "
                  "x_mn.")
            return
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
        # Find the parameteric value where primary ray enters and exit domains
        #   (the primary ray connects the star and the lab)
        # A valid ray should enter and exit the domain, and thus there should be
        # two well defined parameteric values found. We store them as t_en and
        # t_ex, and check later that t_ex < t_en.
        t_en, t_ex = (None for j in range(2))
        for i in range(3):
            if ell[i] != 0:
                if (tol_mn[0]*x_mn[0]<=x_mx[i]*ell[0]/ell[i]<=x_mx[0]*tol_mx[0]
                    and
                    tol_mn[1]*x_mn[1]<=x_mx[i]*ell[1]/ell[i]<=x_mx[1]*tol_mx[1]
                    and
                    tol_mn[2]*x_mn[2]<=x_mx[i]*ell[2]/ell[i]<=x_mx[2]*tol_mx[2]
                   ):
                    if t_en is None:
                        t_en = x_mx[i]/(dist*ell[i])
                    else:
                        t_ex = x_mx[i]/(dist*ell[i])
                if (tol_mn[0]*x_mn[0]<=x_mn[i]*ell[0]/ell[i]<=x_mx[0]*tol_mx[0]
                    and
                    tol_mn[1]*x_mn[1]<=x_mn[i]*ell[1]/ell[i]<=x_mx[1]*tol_mx[1]
                    and
                    tol_mn[2]*x_mn[2]<=x_mn[i]*ell[2]/ell[i]<=x_mx[2]*tol_mx[2]
                   ):
                    if t_en is None:
                        t_en = x_mn[i]/(dist*ell[i])
                    else:
                        t_ex = x_mn[i]/(dist*ell[i])
        if t_en is None or t_ex is None:
            message=('The ray from the observer to the primary did not pass '
                     'through the domain.\nRun a larger simulation, or change '
                     'the orbital parameters to ensure the primary ray passes '
                     'through the domain.\n')
            sys.stderr.write(message)
            sys.exit()
        if (t_en < t_ex):
            t_en, t_ex = t_ex, t_en

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
        self.ray_enter, self.ray_exit = [[[(None,None,None)
                                           for i in range(self.n_x1)]
                                          for j in range(self.n_x2)]
                                         for k in range(2)]
        for j in range(self.x2_min, self.x2_max):
            for i in range(self.x1_min, self.x1_max):
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
                self.ray_enter[j-self.x2_min][i-self.x1_min] = tuple(map(float,(
                    (s_en*ell_p[0]+dist*ell[0]-frame_shift[0]),
                    (s_en*ell_p[1]+dist*ell[1]-frame_shift[1]),
                    (s_en*ell_p[2]+dist*ell[2]-frame_shift[2]))))
                self.ray_exit[j-self.x2_min][i-self.x1_min]  = tuple(map(float,(
                    (s_ex*ell_p[0]+dist*ell[0]-frame_shift[0]),
                    (s_ex*ell_p[1]+dist*ell[1]-frame_shift[1]),
                    (s_ex*ell_p[2]+dist*ell[2]-frame_shift[2]))))
        return