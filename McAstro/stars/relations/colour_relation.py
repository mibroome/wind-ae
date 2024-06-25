from .colour_relation_models import Mamajek2008

__all__ = []

def Mamajek_rotation_rate(BV0, age):
    return Mamajek2008.rotation_rate(BV0, age)
Mamajek_rotation_rate.__doc__ = (Mamajek2008.rotation_rate.__doc__)
__all__ += ['Mamajek_rotation_rate']

def Mamajek_stellar_age(BV0, P_rot):
    return Mamajek2008.stellar_age(BV0, P_rot)
Mamajek_stellar_age.__doc__ = (Mamajek2008.stellar_age.__doc__)
__all__ += ['Mamajek_stellar_age']