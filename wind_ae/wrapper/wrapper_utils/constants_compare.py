import math

## Celestial body masses/radii
M_jupiter = 1.898e30            # Jupiter's mass (g)
R_jupiter = 6.9911e9             # Jupiter's radius (cm)
M_earth = 5.9726e27       # Earth's mass (g)
R_earth = 6.3725e8        # Earth's radius (cm)
M_sun = 1.989e33          # Sun's mass (g)
R_sun = 6.96340e10        # Sun's radius (cm)
au = 1.495978707e13       # Astronomical unit (cm)

## Particle masses
amu = 1.66053906892e-24   # Atomic mass unit (g)
me = 9.1093837139e-28     # Electron mass (g)
mH = 1.007825e0*amu       # Hydrogen atomic mass (g)
mHe = 4.002603e0*amu      # Helium atomic mass (g)

## Ionization energies
e_HI = 13.5984345997e0    # Hydrogen ionization energy (eV)
e_HeI = 24.5873890110e0   # Helium first ionization energy (eV)
e_HeTR = 4.76777448600e0  # Helium triplet ionization energy (eV)
e_HeII = 54.4177655282e0  # Helium second ionization energy (eV)

## Physical constants
kb = 1.38e-16              # Boltzmann constant, CGS (erg/K)
kb_eV = 8.6167e-05         # Boltzmann constant, eV (eV/K)
G = 6.67259e-8             # Gravitational constant (dyn cm^2/g^2)
hp = 6.62620e-27           # Planck's constant, CGS (erg s)
hp_eV = 4.1357e-15         # Planck's constant, eV (eV s)
c_light = 2.99792458e10    # Speed of light (cm/s)


## Physical Constants (source: NIST)
elem_charge = 4.8032068e-10   # elementary charge (esu) [constants.py: e]
Ry = 2.1798741e-11            # Rydberg constant (erg)
mp = 1.672621923e-24          # mass of proton (g) [CODATA 2018]
eV = 1.6021772e-12            # electron volt (erg) -- unit conversion factor,
                               # distinct from the e_HI/e_HeI/... ionization
                               # energies above, which are already in eV

## Unit conversions
ns = 1.0e-9                 # nanosecond (s)
day = 8.64e4                # Julian day (s)
year = 3.15576e7            # Julian year (s)
srday = 8.616409054e4       # sidereal day (s)
sryear = srday*3.6525636e2  # sidereal year (s)
AA = 1.0e-8                 # Angstrom (cm)
pc = 3.086e18                # parsec (cm)
Mb = 1.0e-18                 # megabarn (cm**2)
bar = 1.0e6                  # bar (barye)
# Da (Dalton) is the same physical quantity as amu above; constants.py
# defines it separately (with slightly different CODATA precision) and
# aliases amu = Da. Kept here for completeness, not aliased to amu, so
# the two source precisions stay distinguishable.
Da = 1.66053906660e-24       # Dalton (g)

## Derived Constants
hc = hp*c_light
a0 = hp**2./(4.*math.pi**2.*me*elem_charge**2.)     # Bohr radius (cm)
sig_SB = 2*math.pi**5*kb**4/(15*c_light**2*hp**3)   # Stefan-Boltzmann constant (erg/cm**2/s/K**4)
RyH = mp/(me+mp)*Ry                                 # Hydrogen Rydberg constant (erg)

## Solar values not in the Excel table (source: IAU)
Lsun = 3.828e33          # Solar bolometric luminosity (erg/s/cm**2)
FuvEarth = 14/3          # Earth's recieved uv [10-91.2 nm] flux (erg/s/cm**2)

## Lyman Alpha (Doublet Mean, source: NIST ASD)
Lya_wl = 1215.67e-8   # wavelength (cm)
Lya_nu = c_light/Lya_wl  # frequency (1/s)
Lya_f12 = .4164       # oscillator strength
Lya_gamma = 6.265e8   # natural broadening frequency (1/s)
