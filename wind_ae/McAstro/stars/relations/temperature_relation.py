from .temperature_relation_models import Sung2013
from .temperature_relation_models import Ballesteros2012

__all__ = []

def Sung_teff_to_BV(T_eff, lum_class='V'):
    return Sung2013.Temp2BV(T_eff, lum_class=lum_class)
Sung_teff_to_BV.__doc__ = (Sung2013._Temp2BV.__doc__)
__all__ += ['Sung_teff_to_BV']

def Ballesteros_teff_to_BV(T_eff):
    return Ballesteros2012.Temp2BV(T_eff)
Ballesteros_teff_to_BV.__doc__ = (Ballesteros2012.Temp2BV.__doc__)
__all__ += ['Ballesteros_teff_to_BV']