#!/usr/bin/env python3

# Sources: 
#   masses: https://en.wikipedia.org/wiki/Stellar_classification
#   colors: http://www.vendian.org/mncharity/dir3/starcolor/

class_lowerbound = {'M':0.08, 'K':0.45, 'G':0.80, 'F':1.04,
                    'A':1.40, 'B':2.10, 'O':16}


class_upperbound = {'M':0.45, 'K':0.80, 'G':1.04, 'F':1.40,
                    'A':2.10, 'B':16, 'O':999}


class_color = {'M':'#ffcc6f', 'K':'#ffd2a1', 'G':'#fff4ea', 'F':'#f8f7ff',
               'A':'#cad7ff', 'B':'#aabfff', 'O':'#9bb0ff'}


def mass_class(mass):
    if mass <= 0.08:
        return None
    if mass <= 0.45:
        return 'M'
    if mass <= 0.80:
        return 'K'
    if mass <= 1.04:
        return 'G'
    if mass <= 1.40:
        return 'F'
    if mass <= 2.10:
        return 'A'
    if mass <= 16:
        return 'B'
    else:
        return 'O'


def plot_classifications(ax, labels=True):
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    width = xlims[1]-xlims[0]
    mass = max(0.08, xlims[0])
    while mass < xlims[1]:
        next_mass = class_upperbound[mass_class(mass+0.01)]
        if next_mass > xlims[1]:
            next_mass = xlims[1]
        mass = (mass+next_mass)/2
        if next_mass < xlims[1]:
            ax.axvline(next_mass, ls='--',
                       color=class_color[mass_class(mass)], zorder=0)
        if labels:
            ax.text((mass-xlims[0])/width, 0.9, mass_class(mass), ha='center', 
                    color=class_color[mass_class(mass)], transform=ax.transAxes)
        mass = next_mass