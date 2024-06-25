from .mass_relation_models import Eker2018

__all__ = []

def Eker_mass_luminosity(mass, shift=False):
    return Eker2018.MLR(mass, shift=shift)
Eker_mass_luminosity.__doc__ = (Eker2018._MLR.__doc__)
__all__ += ['Eker_mass_luminosity']


def Eker_luminosity_mass(lum):
    return Eker2018.LMR(lum)
Eker_luminosity_mass.__doc__ = (Eker2018._LMR.__doc__)
__all__ += ['Eker_luminosity_mass']


def Eker_mass_radius(mass, shift=False):
    return Eker2018.MRR(mass, shift=shift)
Eker_mass_radius.__doc__ = (Eker2018._MRR.__doc__)
__all__ += ['Eker_mass_radius']


def Eker_mass_temperature(mass, shift=False):
    return Eker2018.MTR(mass, shift=shift)
Eker_mass_temperature.__doc__ = (Eker2018._MTR.__doc__)
__all__ += ['Eker_mass_temperature']