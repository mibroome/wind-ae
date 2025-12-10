import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from cycler import cycler
import wind_ae.McAstro.atoms.atomic_species as McAtom
from matplotlib.lines import Line2D
from wind_ae.wrapper.wrapper_utils import constants as const
from wind_ae.wrapper.wrapper_utils.spectrum import spectrum

from IPython import display


_SI_prefix = {"y":1e-24, "z":1e-21, "a":1e-18, "f":1e-15, "p": 1e-12,
              "n":1e-9, "u":1e-6, "µ":1e-6, "m":1e-3, "c":1e-2, "d":0.1,
              "h":100, "k":1000, "M":1e6, "G":1e9, "T":1e12, "P":1e15,
              "E":1e18, "Z":1e21, "Y":1e14}
_clestial_prefix = {'Jupiter':7.1492e7, 'Earth':6.3781366e6} #radius in cm

def _convert_unit(unit_in, unit_out):
    if unit_in in _SI_prefix:
        if unit_out in _SI_prefix:
            return _SI_prefix[unit_in]/_SI_prefix[unit_out]
        elif unit_out in _clestial_prefix:
            return _SI_prefix[unit_in]/_clestial_prefix[unit_out]
    elif unit_in in _clestial_prefix:
        if unit_out in _SI_prefix:
            return _clestial_prefix[unit_in]/_SI_prefix[unit_out]
        elif unit_out in _clestial_prefix:
            return _clestial_prefix[unit_in]/_clestial_prefix[unit_out]
    else:
        raise SystemExit('unit prefix error')


def four_panel_plot(windsoln, ax, ax_Ys=None, label=None, alpha=0.8,
                    pvar=[['v','rho'],['T','Ys_HI']], norm=[[1e5,1e0],[1e0,1e0]],
                    radius_prefix=None, 
                    sub_sonic=False, past_rmin=False,
                    sonic_vert=True):
    """
    Produces a velocity, density, temperature, and neutral fraction plot from intermediate solutions while ramping. For a more aesthetic four-panel plot use `quickplot()`.

    Args:
        windsoln: windsoln object to be plotted.
        ax (np.ndarray): The ax array on which plots are placed (must have shape = (2,2)).
        ax_Ys (matplotlib.axes.Axes, optional): Axis for plotting ionization fraction. Defaults to None.
        label (str, optional): Label to be added to ax[0][0]. Defaults to None.
        alpha (float, optional): Alpha of lines plotted. Defaults to 0.8.
        pvar (list of lists of str, optional): Variables to plot. Defaults to [['v','rho'],['T','Ys_HI']].
        norm (list of lists of floats, optional): Values to norm (divide) pvars with. Defaults to [[1e5,1e0],[1e0,1e0]].
        radius_prefix (str, optional): SI prefix for x-axis (None defaults to units of Rp). Defaults to None.
        sub_sonic (bool, optional): Only plot sub-sonic region. Defaults to False.
        past_rmin (bool, optional): Only plot past Rmin region. Defaults to False.
        sonic_vert (bool, optional): Add vertical line at sonic point. Defaults to True.

    Returns:
        None
    """
    if ax.shape != (2,2):
        print("ERROR: four panel plot requires ax with shape (2,2)")
        return
    if radius_prefix in _SI_prefix or radius_prefix in _clestial_prefix:
        radius_norm = _convert_unit('c', radius_prefix)
        radius = radius_norm*windsoln.soln['r']
        radius_norm *= windsoln.Rp
    else:
        if radius_prefix is not None and radius_prefix != 'Rp':
            print('WARNING: radius_prefix not found, defaulting to units of Rp.')
        radius = windsoln.soln_norm['r']
        radius_norm = 1.0
    empty = (True if ax[0][0].lines == [] else False)
    # All plots
    for i in range(len(ax)):
        for j in range(len(ax[0])):
            if i == len(ax) - 1:
                if radius_prefix in _clestial_prefix:
                    ax[i][j].set_xlabel(rf'Radius ($R_\mathrm{{{radius_prefix}}}$)')
                elif radius_prefix in _SI_prefix:
                    ax[i][j].set_xlabel(f'Radius ({radius_prefix}m)')
                else:

                    ax[i][j].set_xlabel(r'Radius ($R_p$)') 
            # Setup x-axis
            xlims = (ax[i][j].get_xlim() if not empty else (np.inf, -np.inf))
            x_overplot = 1.+ 0.01
            xlo = radius[0]
            xhi = radius_norm*windsoln.Rmax
            if sub_sonic or not windsoln.integrate_outward:
                xhi = radius_norm*(windsoln.soln_norm['z'][1]+1.)
            if past_rmin:
                xlo = radius_norm*windsoln.Rmin
            ax[i][j].set_xlim([min(xlims[0], xlo)/x_overplot,
                               x_overplot*max(xlims[1], xhi)])
            ax[i][j].xaxis.set_major_locator(ticker.AutoLocator())
            ax[i][j].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax[i][j].ticklabel_format(axis='x', useMathText=True)
            if i != len(ax)-1:
                ax[i][j].tick_params(axis='x', labelbottom=False)
                ax[i][j].ticklabel_format(axis='x', style='plain')
            if j == len(ax[0])-1:
                ax[i][j].yaxis.set_label_position("right")
                ax[i][j].yaxis.tick_right()
                ax[i][j].yaxis.set_ticks_position("right")
            # Setup y-axis
            ylims = (ax[i][j].get_ylim() if not empty else (np.inf, -np.inf))
            xmask = (radius<xhi)&(radius>=xlo)
            ylo = windsoln.soln[pvar[i][j]][xmask].min()
            yhi = windsoln.soln[pvar[i][j]][xmask].max()
            ymin = min(ylims[0], 0.9*ylo/norm[i][j])
            # if ymin < 0:
            #     ymin = 1e-10
            # print(ymin, max(ylims[1], 1.1*yhi/norm[i][j]))
            # ax[i][j].set_ylim([ymin,
            #                    max(ylims[1], 1.1*yhi/norm[i][j])])
            if ax[i][j].get_yscale() == 'linear':
                ax[i][j].yaxis.set_major_locator(ticker.AutoLocator())
                ax[i][j].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            else:
                ax[i][j].yaxis.set_major_locator(ticker.LogLocator(numticks=80))
            # Do the plotting
            p = ax[i][j].plot(radius[xmask],
                              windsoln.soln[pvar[i][j]][xmask]/norm[i][j], '-',
                              alpha=alpha, zorder=99)
            if i==0 and j==0:
                line_color = p[0].get_color()
                p[0].set_label(label)
            else:
                p[0].set_color(line_color)
            # Plot sonic point as dashed red line (overtop line_color solid)
            vaf = 0.75 # vert_alpha_factor
            if sonic_vert:
                ax[i][j].axvline(radius_norm*(windsoln.soln_norm['z'][1]+1.),
                                 c='red', alpha=vaf*alpha, ls='--', zorder=3)
                ax[i][j].axvline(radius_norm*(windsoln.soln_norm['z'][1]+1.),
                                 c=line_color, alpha=vaf*alpha, ls='-', zorder=2)
            if j == len(ax[0])-1:
                ax[i][j].tick_params(axis='both', which='both', right=True)
            if pvar[i][j] == 'v':
                ax[i][j].set_ylabel(r'Velocity (km$\,$s$^{-1}$)')
                ax[i][j].set_ylim([0, ax[0][0].get_ylim()[1]])
            elif pvar[i][j] == 'rho':
                ax[i][j].set_yscale('log')
                ax[i][j].set_ylabel(r'Density (g$\,$cm$^{-3}$)')
            elif pvar[i][j] ==  'T':
                ax[i][j].set_ylabel(r'Temperature (K)')
                ax[i][j].set_yscale('log')
            elif pvar[i][j] == 'Ys_HI':
                ax[i][j].set_ylabel(r'H Neutral Fraction')
                ax[i][j].set_yscale('log')
            else:
                print("WARNING: pvar {:s} has no preset labeling."
                      .format(pvar[i][j]))
            if j == len(ax[0])-1:
                ax[i][j].set_ylabel(ax[i][j].yaxis.get_label_text(),
                                    rotation=270, va="bottom")
    return

def _custom_rc_params(line_color,nspecies):
    #setting fontsize based on number of species in legends
    if nspecies > 4:  #smaller legends
        fontsize=11
        columns = 2
    else:
        fontsize=14
        columns=1
    
    #choosing a colormap that is intuitive for the given line color. 
    #The colormap is a dynamic way to have and plot an indefinite number of species
    if line_color in ['black','k','tab:gray','dimgray','dimgrey','grey','gray','darkgray','darkgrey',
                      'silver','lightgrey','lightgray']:
        colormap=plt.cm.Greys_r(np.linspace(0,1,nspecies+1))
    if line_color in ['r','red','tab:red','lightcoral','indianred','brown',
                      'firebrick','maroon','darkred','salmon','tomato','darksalmon','coral',
                      'orangered','lightsalmon','crimson']:
        colormap=plt.cm.OrRd(np.linspace(0,1,nspecies+1))
    if line_color in ['tab:orange','darkorange','orange']:
        colormap=plt.cm.Oranges(np.linspace(0,1,nspecies+1))
    if line_color in ['y','yellow','gold','darkgoldenrod','goldenrod','khaki','tab:olive']:
        colormap=plt.cm.YlOrBr(np.linspace(0,1,nspecies+1))
    if line_color in ['tab:brown','sienna','chocolate','saddlebrown','sandybrown','peru']:
        colormap=plt.cm.copper_r(np.linspace(0,1,nspecies+1))
    if line_color in ['g','green','tab:green','olivedrab','yellowgreen','darkolivegreen',
                      'greenyellow','chartreuse','forestgreen','limegreen',
                      'darkgreen','lime','seagreen','mediumaquamarine',
                      'lawngreen','palegreen','mediumseagreen','springgreen','mediumspringgreen']:
        colormap=plt.cm.Greens(np.linspace(0,1,nspecies+1))
    if line_color in ['aquamarine','turquoise','lightseagreen',
                      'mediumturquoise','paleturquoise','darkslategrey',
                      'darkslategrey','teal','darkcyan','cyan','aqua',
                      'darkturquoise','cadetblue','tab:cyan','c']:
        colormap=plt.cm.GnBu(np.linspace(0,1,nspecies+1))
    if line_color in ['b','tab:blue','blue','powderblue',
                      'lightblue','deepskyblue','skyblue','lightskyblue',
                      'steelblue','dodgerblue','lightslategrey',
                      'lightslategray','slategray','slategrey',
                     'lightsteelblue','cornflowerblue',
                      'royalblue','midnightblue','midnightblue','navy','darkblue','mediumblue']:
        colormap=plt.cm.Blues(np.linspace(0,1,nspecies+1))
    if line_color in ['tab:purple','purple','slateblue','darkslateblue',
                      'mediumslateblue','mediumpurple','rebeccapurple',
                     'blueviolet','indigo','darkorchid','darkviolet',
                      'mediumorchid','thistle','plum','violet']:
        colormap=plt.cm.Purples(np.linspace(0,1,nspecies+1))
    if line_color in ['tab:pink','darkmagenta','fuschia','magenta',
                      'orchid','mediumvioletred','deeppink','hotpink',
                     'palevioletred','pink','lightpink','m']:
        colormap=plt.cm.spring_r(np.linspace(0,1,nspecies+1))
        
    custom_cycler = (cycler(linestyle=['-', '--', ':', '-.',
                                       (0, (1, 10)),
                                       (0, (1, 1)),
                                       (0, (3, 5, 1, 5, 1, 5)),
                                       (5, (10, 3)),
                                       (0, (3, 1, 1, 1)),
                                       (0, (1, 1)),
                                       (0, (5, 10)),
                                       (0, (3, 10, 1, 10)),
                                       (0, (5, 5)),
                                       (0, (3, 5, 1, 5)),
                                       (0, (3, 10, 1, 10, 1, 10)),
                                       (0, (5, 1)),
                                       (0, (3, 1, 1, 1, 1, 1))]))
    
    return colormap, custom_cycler, fontsize, columns



def quick_plot(soln, Mdot_legend=True, c='k', ls='-', label='',label_dim=[0,1.3,2],
             ion_label=True,first_plotted=True, ax=0): 
    """
    Plots density (g/cm³), temperature (K), velocity (10 km/s), and ionization fraction as a function of r (Rp).

    Args:
        soln (windsoln): Wind solution object (sim.windsoln).
        Mdot_legend (bool, optional): If True, put Mdot in legend of plot. Else, just prints. Defaults to True.
        c (str, optional): Line color. Defaults to 'k'.
        ls (str, optional): Line style. Defaults to '-'.
        label (str, optional): Label for the plot. Defaults to ''.
        label_dim (list, optional): Location of label and ncols [x, y, ncols]. Defaults to [0, 1.3, 2].
        ion_label (bool, optional): If True, show ionization legend. Defaults to True.
        first_plotted (bool, optional): True if this is the first of many plots on the same axes. Defaults to True.
        ax (matplotlib.axes.Axes, optional): Axes object to plot on. Defaults to 0.

    Returns:
        matplotlib.axes.Axes: Axes object (if first_plotted=True).
        str: Title string summarizing species and mass fractions.
    """
    try: # Check if R_cori has be calculated. If not, calculate all post-facto user variables
        soln.R_cori
    except AttributeError:
        if soln.integrate_outward == 0:
            soln.R_cori = 1e10
        else:
            soln.add_user_vars()
    
    radius_norm = 1.0
    alpha=0.5
    radius = soln.soln_norm['r']
    R_H = soln.semimajor*(soln.Mp/(3*soln.Mstar))**(1/3) / soln.Rp
    nspecies = soln.nspecies
    
    if len(label) == 0:
        Mstr = soln.Mp/const.Mearth 
        Rstr = soln.Rp/const.Rearth
        Fstr = str(int(soln.Ftot))
        spectrum = soln.spec_src_file
        label = f'{Mstr:.2f}Mearth_{Rstr:.2f}Rearth_{Fstr:s}F_{spectrum:s}'
    
    colormap,custom_cycler,fontsize,columns = _custom_rc_params(c,nspecies)

    stack=2
    if first_plotted==True:
        fig, ax = plt.subplots(stack,2,sharex=True,figsize=[11,6])
        fig.subplots_adjust(hspace=0)
        ax[1,0].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                        c=c, alpha=alpha, ls='--', zorder=3,label='Sonic Point')
        ax[1,0].axvline(R_H,ls=':',c=c,label='Hill Radius')
        if soln.Rmax > soln.R_cori:
            ax[1,0].fill_betweenx((0,np.max(soln.soln['rho'])),
                                  soln.R_cori,soln.Rmax,alpha=0.3,color=c,
                                  label=r'R$_{cori}$ = %.2fR$_p$' %soln.R_cori)
        ax[1,0].legend()
    #Density (+Hill Radius and sonic point for legend purposes)
    ax[0,0].semilogy(soln.soln_norm['r'], soln.soln['rho'],c=c,ls=ls,
                     label=label)
    ax[0,0].legend(loc='upper left', 
                   bbox_to_anchor=(label_dim[0],label_dim[1]),ncol=label_dim[2],frameon=False)    
    ax[0,0].set_ylabel(r'Density (g/cm$^3$)')
 
    #Velocity 
    ax[1,0].plot(soln.soln_norm['r'], soln.soln['v']/1e6,c=c,ls=ls)
    ax[1,0].set_ylabel(r'Velocity (10 km/s)')
    
    #Temperature
    mdot = (soln.Mdot)#/const.Msun)*3.154e+7 #in Msun/year
    print(f'*****{label:s} Mdot = {mdot:.2e} g/s ******')
    ax[0,1].plot(soln.soln_norm['r'], soln.soln['T']/1000,c=c,
                     ls=ls,label=r'$\dot{M}$ = %.1e g/s' %(mdot))
    ax[0,1].set_ylim(np.min(soln.soln['T']/1000)*0.97,np.max(soln.soln['T']/1000)*1.03)
    ax[0,1].set_ylabel(r'Temperature (1000 K)')
    ax[0,1].set_yscale('log')
    ax[0,1].get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax[0,1].get_yaxis().set_minor_formatter(ticker.ScalarFormatter())
    if Mdot_legend==True:
        leg = ax[0,1].legend(edgecolor='k')
        leg.get_frame().set_linewidth(1.5)    #Ionization fraction
    title = ''
    ax[1,1].set_prop_cycle(custom_cycler)
    try:
        soln.soln['n_H']
    except KeyError:
        soln.add_user_vars()
    for j,spname in enumerate(soln.species_list):
        spname = spname.replace(' ','')
        ax[1,1].semilogy(soln.soln_norm['r'],
                         1-soln.soln['Ys_'+spname],label=spname,c=c)
        title += spname+': %.2f, '%(soln.HX[j])
    ax[1,1].set_ylabel(r'Ionization Fraction')
    ax[1,1].set_ylim((1e-2,0.99999999))
#     ax[1,1].set_yticks([1e-3,1e-2,1e-1],labels=[r'$10^{\text{-}3}$',r'$10^{-2}$',r'$10^{-1}$'])
    if ion_label==True:
        ax[1,1].legend(fontsize=fontsize,ncol=columns)
    #Plotting shared lines
    if soln.Rmax > soln.R_cori:
        if first_plotted==False:
            ax[0,0].fill_betweenx((0,np.max(soln.soln['rho'])),
                                  soln.R_cori,soln.Rmax,alpha=0.3,color=c)
        ax[0,1].fill_betweenx((0,np.max(soln.soln['T'])),
                              soln.R_cori,soln.Rmax,alpha=0.3,color=c)
        ax[1,0].fill_betweenx((0,np.max(soln.soln['v']/(10*100*1000))),
                              soln.R_cori,soln.Rmax,alpha=0.3,color=c)
        ax[1,1].fill_betweenx((0,1),soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    for k in range(2):
        ax[1,k].set_xlabel(r'Radius (R$_p$)')
        for m in range(stack):
            ax[m,k].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                             c=c, alpha=alpha, ls='--', zorder=3)
            ax[m,k].axvline(R_H,ls=':',c=c)
    plt.gca().set_xlim(left=soln.Rmin)
    if first_plotted==True:
        return ax


    
def six_panel_plot(soln,Mdot_legend=True,c='k',ls='-',label='',
                   label_dim=[0,1.3,2],ion_label=True,first_plotted=True,ax=0): 
    '''
    Plots density (g/cm3), temperature (K), velocity (10 km/s), ionization fraction, column density (g/cm2), and number density (1/cm2), as a function of r (Rp).
        
    Args:
        soln - windsoln object (sim.windsoln)
        Mdot_legend - Bool; if True, put Mdot in legend of plot. Else, just prints.
        c - str; line color 
        ls - str; line style
        label - str; line label 
        label_dim - list; default=[0,1.3,2]. Location of label and ncols [x,y,ncols]. 
        first_plotted - Bool; True if this the first of many OR the ONLY SixPlot 
                        to be plotted on the same axes. 
        ax - matplotlib axis obj; if first_plotted=False, provide axis object so this 
             will be be plotted on desired figure with other simulations for comparison
    Returns:
        ax - axes object (if first_plotted=True)
        
    Example:
        ax1 = SixPlot(sim1.windsoln, first_plotted=True)
        SixPlot(sim2.windsoln, ax=ax1)
        SixPlot(sim3.windsoln, ax=ax1)
    '''
    try: # Check if R_cori has be calculated
        soln.R_cori
    except AttributeError:
        if soln.integrate_outward == 0:
            soln.R_cori = 1e10
        else:
            soln.add_user_vars()
    
    if len(label) == 0:
        Mstr = soln.Mp/const.Mearth 
        Rstr = soln.Rp/const.Rearth
        Fstr = str(int(soln.Ftot))
        spectrum = soln.spec_src_file
        label = f'{Mstr:.2f}Mearth_{Rstr:.2f}Rearth_{Fstr:s}F_{spectrum:s}'
    radius_norm = 1.0
    alpha=0.5
    radius = soln.soln_norm['r']
    # minrad = float(radius[soln.soln['T'] == np.min(soln.soln['T'])])
    minrad = float(radius.iloc[soln.soln['T'].idxmin()])
    R_H = soln.semimajor*(soln.Mp/(3*soln.Mstar))**(1/3) / soln.Rp
    nspecies = soln.nspecies
    colormap,custom_cycler,fontsize,columns = _custom_rc_params(c,nspecies)

    stack=3
    if first_plotted==True:
        fig, ax = plt.subplots(stack,2,sharex=True,figsize=[11,9])
        fig.subplots_adjust(hspace=0)
        ax[1,0].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                     c=c, alpha=alpha, ls='--', zorder=3,label='Sonic Point')
        ax[1,0].axvline(R_H,ls=':',c=c,label='Hill Radius')
        if soln.Rmax > soln.R_cori:    
            ax[1,0].fill_betweenx((0,np.max(soln.soln['rho'])),
                                  soln.R_cori,soln.Rmax,alpha=0.3,color=c,
                                  label=r'R$_{cori}$ = %.2fR$_p$' %soln.R_cori)
    
    #Density (+Hill Radius and sonic point for legend purposes)
    if first_plotted == False and soln.Rmax > soln.R_cori:
        ax[0,0].fill_betweenx((0,np.max(soln.soln['rho'])),
                              soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    ax[0,0].semilogy(soln.soln_norm['r'], soln.soln['rho'],c=c,
                     ls=ls,label=label)
    ax[0,0].legend(loc='upper left', 
                   bbox_to_anchor=(label_dim[0],label_dim[1]),ncol=label_dim[2],frameon=False)    
    ax[0,0].set_ylabel(r'Density (g/cm$^3$)')
    #Velocity 
    ax[1,0].plot(soln.soln_norm['r'], soln.soln['v']/1e6,c=c,ls=ls)
    ax[1,0].legend()
    
    if soln.Rmax > soln.R_cori:
        ax[1,0].fill_betweenx((0,np.max(soln.soln['v']/1e6)),
                              soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    ax[1,0].set_ylabel(r'Velocity (10 km/s)')

    #Temperature
    mdot = (soln.Mdot)#/const.Msun)*3.154e+7 #in Msun/year
    print(f'***** {label:s}: Mdot = {mdot:.2e} g/s *****')
    ax[2,0].plot(soln.soln_norm['r'], soln.soln['T']/1e3,c=c,
                 ls=ls,label=r'$\dot{M}$ = %.1e g/s' %(mdot))
    ax[2,0].set_ylabel(r'Temperature (1000 K)')
    ax[2,0].set_ylim(np.min(soln.soln['T']/1000)*0.97,np.max(soln.soln['T']/1000)*1.03)
    ax[2,0].set_yscale('log')
    ax[2,0].get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax[2,0].get_yaxis().set_minor_formatter(ticker.ScalarFormatter())
    

    if soln.Rmax > soln.R_cori:
        ax[2,0].fill_betweenx((0,np.max(soln.soln['T'])),
                              soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    if Mdot_legend==True:
        leg = ax[2,0].legend(edgecolor='k')
        leg.get_frame().set_linewidth(1.5)

    #Number density
    try:
        soln.soln['n_H']
    except KeyError:
        soln.add_user_vars()
    species_copy = soln.species_list.copy()
    spaced = McAtom.formatting_species_list(species_copy)
    ax[0,1].set_prop_cycle(custom_cycler)
    max_n = []
    for j in range(soln.nspecies):
        element_name = ((spaced[j]).split())[0]
        lowest_state = ((spaced[j]).split())[1] #will have to adapt for elements with more than 1 ionization state
        #converting to arabic numbers to make future multiple-ionization-state version of the code easier
        highest_state = McAtom.arabic_to_roman(McAtom.roman_to_arabic(lowest_state)+1)
        total = 'n_'+element_name
        neutral = 'n_'+element_name+lowest_state
        ionized = 'n_'+element_name+highest_state
        ax[0,1].semilogy(soln.soln_norm['r'],soln.soln[total],ls='-',c=colormap[j])
        ax[0,1].semilogy(soln.soln_norm['r'],soln.soln[neutral],ls='--',
                         c=colormap[j],label=r'$%s_{%s}$'%(element_name,lowest_state))
        ax[0,1].semilogy(soln.soln_norm['r'],soln.soln[ionized],
                         ls=':',c=colormap[j],label=r'$%s_{%s}$'%(element_name,highest_state))
        max_n = np.append(max_n,np.max(soln.soln[total]))
    if soln.Rmax > soln.R_cori:
        ax[0,1].fill_betweenx((0,np.max(max_n)),soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    ax[0,1].legend(ncol=2,fontsize=fontsize) 
    ax[0,1].set_ylim(bottom=1e3)
    ax[0,1].set_ylabel(r'n (cm$^{-3}$)')

    #Ionization fraction
    title = ''
    ax[1,1].set_prop_cycle(custom_cycler)
    for j,spname in enumerate(soln.species_list):
        spname = spname.replace(' ','')
        ax[1,1].semilogy(soln.soln_norm['r'],1-soln.soln['Ys_'+spname],label=spname,c=c)
        title += spname+': %.2f, '%(soln.HX[j])
    ax[1,1].set_ylabel(r'Ionization Fraction')
    if soln.Rmax > soln.R_cori:
        ax[1,1].fill_betweenx((0,1),soln.R_cori,soln.Rmax,alpha=0.3,color=c)
    ax[1,1].set_ylim((1e-3,0.99999999))
    if ion_label==True:
        ax[1,1].legend(fontsize=fontsize,ncol=columns)

    #Column Density
    ax[2,1].set_prop_cycle(custom_cycler)
    max_ncol = []
    for spname in soln.species_list:
        spname = spname.replace(' ','')
        ax[2,1].semilogy(soln.soln_norm['r'],soln.soln['Ncol_'+spname],label=spname,c=c)
        max_ncol = np.append(max_ncol,np.max(soln.soln['Ncol_'+spname]))
    ax[2,1].set_ylabel(r'Ncol (cm$^{-2}$)')
    ax[2,1].legend(fontsize=fontsize,ncol=columns)
    if soln.Rmax > soln.R_cori:
        ax[2,1].fill_betweenx((0,np.max(max_ncol)),soln.R_cori,soln.Rmax,alpha=0.3,color=c)

    #Plotting shared lines
    for k in range(2):
        ax[2,k].set_xlabel(r'Radius (R$_p$)')
        for m in range(stack):
            ax[m,k].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                             c=c, alpha=alpha, ls='--', zorder=3)
            ax[m,k].axvline(R_H,ls=':',c=c)
#             if soln.Rmax > soln.R_cori:
#             if soln.Rmax > 10:
#                 ax[m,k].fill_betweenx((0,1e100),soln.R_cori,soln.Rmax,alpha=0.3,color='k')
#                 custom_lines = [Line2D([0], [0], color='tab:gray', lw=4)]
#                 ax[0,0].legend(custom_lines, [r'R$_{cori}$ = %.2f' %soln.R_cori])
    plt.gca().set_xlim(left=soln.Rmin)
    plt.gca().set_xlim(right=soln.Rmax)
    if first_plotted==True:
        return ax    
    

def energy_plot(windsoln, ax=0, alpha=0.8, all_terms=False, 
                CII_line_cool=False, CIII_line_cool=False, OII_line_cool=False,OIII_line_cool=False, 
                legend=True,sub_sonic=True):
    """ Plots energy balance terms used in the energy equation (Broome et al. 2025)

    Args: 
        windsoln: The wind solution object containing the simulation data. ("sim.windsoln")
        ax: The axis to plot on (default is 0, which creates a new figure)
        alpha: Transparency level for the plot lines (default is 0.8). 
                Useful when overplotting multiple on same axes
        all_terms: If True, plot terms not included in Wind-AE, e.g., free-free cooling, (default is False)
        CII_line_cool: If True, include CII line cooling terms (default is False)
        CIII_line_cool: If True, include CIII line cooling terms (default is False)
        OII_line_cool: If True, include OII line cooling terms (default is False)
        OIII_line_cool: If True, include OIII line cooling terms (default is False)
        legend: If True, display the legend (default is True)
        sub_sonic: If True, sets x-axis upper limit at sonic point radius

    Returns:
        None
    """
    try:
        windsoln.soln['heat_ion']
    except KeyError:
        windsoln.add_user_vars()
    ## Uncommenting this block allows for live updating of plots in some versions of Jupyter notebooks
#     display.display(plt.gcf())
#     display.clear_output(wait=True)
#     plt.clf()
#     if ax==0:
#         fig,ax = plt.subplots()
    if ax==0:
        fig,ax = plt.subplots()
    ncols=1
    fontsize=14
    r = windsoln.soln_norm['r'][1:]
    ax.plot(r, windsoln.soln['heat_ion'][1:], '-',
            alpha=alpha, c='red',label='ionization heating')
    if windsoln.lyacool != 0:
        try:
            windsoln.soln['cool_lyman']
        except KeyError:
            windsoln.add_user_vars()
        ax.plot(r, -windsoln.soln['cool_lyman'][1:], '--',
                alpha=alpha, c='tab:cyan', label='Lyman-alpha cooling')
    ax.plot(r, -windsoln.soln['cool_PdV'][1:], '--',
            alpha=alpha, c='darkblue',label='PdV cooling')
    if windsoln.bolo_heat_cool != 0:
        ax.plot(r, windsoln.soln['boloheat'][1:], ls=(0, (3, 1, 1, 1, 1, 1)),
            alpha=alpha, c='lightcoral',label='bolometric heating')
        ax.plot(r, -windsoln.soln['bolocool'][1:], ls=(0, (3, 1, 1, 1, 1, 1)),
            alpha=alpha, c='lightblue',label='bolometric cooling')
    ax.plot(r, windsoln.soln['heat_advect'][1:], '-.',
            alpha=alpha, c='maroon', label='advective heating')
    ax.plot(r, -windsoln.soln['heat_advect'][1:], '-.',
            alpha=alpha, c='dodgerblue', label='advective cooling')
    # if all_terms:
    ncols=2
    fontsize=12
    ax.plot(r, -windsoln.soln['cool_rec'][1:], ':',
        alpha=alpha, c='teal',label='recombination cooling')
    if all_terms:
        ax.plot(r, windsoln.soln['cool_cond'][1:], ':',
            alpha=alpha, c='navy',label='conductive cooling')
    ax.plot(r, -windsoln.soln['cool_cond'][1:], ':',
        alpha=alpha, c='lightsalmon',label='conductive heating')
#         ax.plot(r, -windsoln.soln['cool_free'][1:], ':',
#             alpha=alpha, c='paleturquoise',label='free-free cooling')
    if CII_line_cool:
        ncols=2
        fontsize=12
        line_str = ['1570000','2326','1334']
        for i,line in enumerate(line_str):
            alpha = ((len(line_str))-i)/len(line_str)
            ax.plot(r, -windsoln.soln['cool_CII_'+line+'A'][1:], '-',c='indigo',
                    alpha=alpha, label=r'CII %s$\mathrm{\AA}$'%line)
    if CIII_line_cool:
        ncols=2
        fontsize=12
        line_str = ['1910','977']
        for i,line in enumerate(line_str):
            alpha = ((len(line_str))-i)/len(line_str)
            ax.plot(r, -windsoln.soln['cool_CIII_'+line+'A'][1:], '-',c='tab:purple',
                    alpha=alpha, label=r'CIII %s$\mathrm{\AA}$'%line)
    if OII_line_cool:
        ncols=2
        fontsize=12
        line_str = ['7320','3727','2741','834']
        for i,line in enumerate(line_str):
            alpha = ((len(line_str))-i)/len(line_str)
            ax.plot(r, -windsoln.soln['cool_OII_'+line+'A'][1:], '-',c='slategrey',
                    alpha=alpha, label=r'OII %s$\mathrm{\AA}$'%line)
    if OIII_line_cool:
        ncols=2
        fontsize=12
        line_str = ['520000','5000','166','84']
        for i,line in enumerate(line_str):
            alpha = ((len(line_str))-i)/len(line_str)
            ax.plot(r, -windsoln.soln['cool_OIII_'+line+'A'][1:], '-',c='slateblue',
                    alpha=alpha, label=r'OIII %s$\mathrm{\AA}$'%line)
 
        
    peak = 10**np.ceil(np.log10(1.1*windsoln.soln['heat_ion'].max()))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'Energy Rates (erg s$^{-1}$ cm$^{-3}$)')
    ax.set_xlabel(r'Radius ($R_p$)')
    xlo, xhi = windsoln.Rmin-0.02, windsoln.Rmax
    if sub_sonic or not windsoln.integrate_outward:
        xhi = windsoln.R_sp
    ax.set_xlim([xlo, xhi])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(ticker.ScalarFormatter())
    ax.set_ylim([peak/1e6, peak])
    if legend:
        ax.legend(ncol=ncols,bbox_to_anchor=(1.05,0),loc='lower left',fontsize=fontsize)
    return


def coriolis_plot(windsoln, ax, var='position', label_system=None):
    """    
    Args:
        windsoln: The wind solution object containing the data to plot.
        ax: The matplotlib axis to plot on.
        var: The variable to plot ('position', 'velocity', 'pos_angle', 'vel_angle').
        label_system: Optional label for the color legend.
    """
    windsoln.add_user_vars()
    s_arr = np.linspace(windsoln.R_sp, windsoln.Rmax, 1000, endpoint=False)
    legend = ax.get_legend()
    labels = [] if legend is None else [str(x._text) for x in legend.texts]
    handles = [] if legend is None else legend.legendHandles
    cori_vert_color = 'm'
    if var == 'position':
        p1, = ax.plot(s_arr, windsoln.cori_pos(s_arr)[0])
        p2, = ax.plot(s_arr, windsoln.cori_pos(s_arr)[1])
        p3, = ax.plot(s_arr, np.sqrt(windsoln.cori_pos(s_arr)[0]**2
                               +windsoln.cori_pos(s_arr)[1]**2))
        p4, = ax.plot(s_arr, s_arr, ls=':')
        ax.legend(handles+[p1, p2, p3, p4],
                  labels+[r'Streamline $x(s)$ coord', r'Streamline $y(s)$ coord position', r'Streamline radial distance ($r(s)$)', r'Substellar radial distance ($s$)'])
        ax.set_xlabel(r'Streamline distance ($R_p$)')
        ax.set_ylabel(r'Distance from core ($R_p$)')
    elif var == 'velocity':
        p1, = ax.plot(s_arr, windsoln.cori_vel(s_arr)[0]/1e5)
        p2, = ax.plot(s_arr, windsoln.cori_vel(s_arr)[1]/1e5)
        p3, = ax.plot(s_arr, windsoln.v_fit(s_arr)/1e5)
        # Sanity check that v = sqrt(v_x^2+v_y^2)
        p4, = ax.plot(s_arr, np.sqrt(windsoln.cori_vel(s_arr)[0]**2
                               +windsoln.cori_vel(s_arr)[1]**2)/1e5, ls=':')

        ax.legend(handles+[p1, p2, (p3, p4)],
                  labels+[r'$v_x$', r'$v_y$', r'$v_{\mathrm{tot}}$'])
        ax.set_xlabel(r'Streamline distance ($R_p$)')
        ax.set_ylabel(r'Velocity (km/s)')
    elif var == 'pos_angle':
        # arctan2(y,x)
        phi = np.arctan2(windsoln.cori_pos(s_arr)[1],
                         windsoln.cori_pos(s_arr)[0])
        if phi[0] < - np.pi/2.:
            for i, p in enumerate(phi):
                if p > -np.pi/2.:
                    phi[i] -= 2.*np.pi
        p1, = ax.plot(s_arr, phi[0]-phi)
        ax.axhline(np.pi/4., c='r', ls='--', zorder=0)
        ax.legend(handles+[p1],
                  labels+[r'$-\Delta\phi_{\mathrm{pos.}}$'])
        ax.set_xlabel(r'Streamline distance ($R_p$)')
        ax.set_ylabel(r'Deflection angle (radians)')
    elif var == 'vel_angle':
        # arctan2(y,x)
        phi = np.arctan2(windsoln.cori_vel(s_arr)[1],
                         windsoln.cori_vel(s_arr)[0])
        if phi[0] < - np.pi/2.:
            for i, p in enumerate(phi):
                if p > -np.pi/2.:
                    phi[i] -= 2.*np.pi
        p1, = ax.plot(s_arr, phi[0]-phi)
        ax.axhline(np.pi/4., c='r', ls='--', zorder=0)
        if label_system is not None:
            ax.legend(handles+[p1], labels+[label_system])
            cori_vert_color = p1.get_color()
        else:
            ax.legend(handles+[p1],
                      labels+[r'$-\Delta\phi_{\mathrm{vel.}}$'])
        ax.set_xlabel(r'Streamline distance ($R_p$)')
        ax.set_ylabel(r'Deflection angle (radians)')
    else:
        print("ERROR: Do not recognize plotting variable: {:s}\n"
              "Aviable variables: 'position', 'velocity', "
              "'pos_angle', 'vel_angle'\n"
             .format(var))
    ax.axvline(windsoln.R_cori, c=cori_vert_color, ls='--', zorder=0)
    return


def ballistic_plot(windsoln, ax, show_star=False, t_end=None):
    """
    Plot ballistic trajectories. Windsoln must be integrated out (if not, run sim.integrate_out())
    There is currently no stellar wind term.

    Args:
        ax (matplotlib.axes.Axes): The axes to plot on.
        show_star (bool): Whether to show the star marker at the semimajor axis.
        t_end (float): The end time for the ballistic trajectories.

    Returns:
        None
    """
    windsoln.add_user_vars()
    if t_end is None:
        t_end = np.log10(windsoln.stream_time(windsoln.R_cori)
                         -windsoln.stream_time(windsoln.R_sp))
    if show_star:
        def semimajor_ball_time(t):
            return (windsoln.semimajor
                    +windsoln.ballistic_sols[int((windsoln.n_ball)/1.75)]
                    .sol(t)[0])
        t_end = np.log10(optimize.fsolve(semimajor_ball_time, 10**t_end)[0])
        ax.plot(-windsoln.semimajor/windsoln.Rp, 0, 'y*')
    time = np.logspace(2, 1.01*t_end, 2000)
    # Plot ballistic trajectories
    for i, bs in enumerate(windsoln.ballistic_sols):
        ax.plot(bs.sol(time)[0,:]/windsoln.Rp, bs.sol(time)[1,:]/windsoln.Rp,
                alpha=0.7, ls=':')
    # Plot estimated sub-stellar trajectory
    s = np.linspace(windsoln.R_sp, windsoln.Rmax, 10000)
    ax.plot(windsoln.cori_pos(s)[0], windsoln.cori_pos(s)[1], 'b-')
    # Overlay length scales
    circle1 = plt.Circle((0, 0), windsoln.R_sp, color='r', lw=5, fill=False,
                         zorder=0, alpha=0.8)
    circle2 = plt.Circle((0, 0), windsoln.R_cori, color='m', lw=5, fill=False,
                         zorder=0, alpha=0.8)
    circle3 = plt.Circle((0, 0), windsoln.Rmin, color='k', lw=5, zorder=10,
                         alpha=0.8)
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.add_artist(circle3)
    ax.set_xlabel(r'x ($R_p$)')
    ax.set_ylabel(r'y ($R_p$)')
    ax.set_aspect('equal')
    ax.set_xlim([*ax.get_ylim()])


def collisionality_plot(windsoln, ax, legend=True, same_colors=False,
                        sub_sonic=False, past_rmin=False,
                        sonic_vert=True, exo_vert=False):
    """Deprecated (requires updates to Knudsen number calculations).
    """
    empty = (True if ax.lines == [] else False)
    nlines = 3+sonic_vert+exo_vert
    nplotted = len(ax.lines)//nlines
    xlo = windsoln.soln_norm['r'][0]
    xhi = windsoln.Rmax
    if sub_sonic or not windsoln.integrate_outward:
        xhi = 0.997*(windsoln.soln_norm['z'][1]+1.)
    if past_rmin:
        xlo = windsoln.Rmin
    xlims = (ax.get_xlim() if not empty else (np.inf, -np.inf))
    x_overplot = 1.01
    ax.set_xlim([min(xlims[0], xlo)/x_overplot, x_overplot*max(xlims[1], xhi)])
    radius = windsoln.soln_norm['r']
    xmask = (radius<xhi)&(radius>=xlo)
    max_mask = (windsoln.soln['Kn_ion_HI'][xmask]<1e6)
    # plot Knudsen numbers
    p = [None,]*3
    p[0] = ax.plot(radius[xmask], windsoln.soln['Kn_hb'][xmask], ls='-.',
                   label='Hardbody', zorder=1*(nplotted+1))
    p[1] = ax.plot(radius[xmask], windsoln.soln['Kn_Co'][xmask], ls='--', lw=5,
                   label='Coulomb', zorder=2*(nplotted+1))
    p[2] = ax.plot(radius[xmask][max_mask],
                   windsoln.soln['Kn_ion_HI'][xmask][max_mask],
                   ls=':', label='Ionization ', zorder=3*(nplotted+1))
    if same_colors:
        samecolor = plt.rcParams['axes.prop_cycle'].by_key()['color'][nplotted]
        for i in range(3):
            p[i][0].set_color(samecolor)
    if sonic_vert:
        ax.axvline(windsoln.R_sp, ls='--', c='r', zorder=0)
    if exo_vert:
        ax.axvline(windsoln.r_exo, ls=':', c='c', zorder=0)
    if legend:
        ax.legend()
    ax.set_yscale('log')
    ax.set_ylabel('Knudsen Number')
    ax.set_xlabel(r'Radius ($R_p$)')
    # Gradientx
    cmap = plt.cm.Reds
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = colors.ListedColormap(my_cmap)
    x_lims, y_lims = ax.get_xlim(), ax.get_ylim()
    if (windsoln.soln['Kn_hb'][xmask].max() <= 1e2 and
        windsoln.soln['Kn_Co'][xmask].max() >= 1e2):
        y_lims = (y_lims[0], 1e2)
    if y_lims[1] > 1e-2:
        # aspect='auto' prevents the code from being a baby about the log scale
        ax.imshow(np.linspace(1e-2, 1e2, 10000).reshape(10000, 1), origin='lower',
                  cmap=my_cmap, interpolation='bicubic', aspect='auto',
                  norm=colors.LogNorm(vmin=5e-2, vmax=5e1),
                  extent=(x_lims[0],x_lims[1],1e-2,1e2), zorder=-1)
    if y_lims[1] > 1e2:
        ax.imshow([[1,1],[1,1]], cmap=plt.cm.Reds_r, aspect='auto',
                  extent=(x_lims[0],x_lims[1],1e2,y_lims[1]), zorder=-2)
    ax.set_ylim(y_lims)
    return

def spectrum_plot(windsoln, var='F_wl',xaxis='wl',semimajor_au=1.0,highlight_euv=True,
                  wl_norm=1e-7,print_warnings=False):    
    """
    Plot the spectrum. Displays the observations, smoothed, and spans.  

    Args:
        var (str): Which variable plotted, energy ('F_wl') or number ('Phi_wl')
        xaxis (str): 'wl' or 'energy'; x variable. Wavelength in nm or energy in eV 
        semimajor_au (float): semimajor axis in units of au to which to scale the spectrum
        highlight_euv (bool): default=True; highlights EUV and XUV range and prints APPROXIMATE fluxes in ergs/s/cm2 each range.
        wl_norm (float): default=1e-7; wavelength normalization factor to convert nm to cm.
        print_warnings (bool): default=False; prints any warnings generated when loading spectrum.

    Returns:
        (fig, ax): The figure and axis object plotted on
    """
    if windsoln.spec_src_file != 'scaled-solar':
        spec = spectrum(lisird=False,
                        spectrum_file=windsoln.spec_src_file,
                        wl_norm=wl_norm,print_warning=print_warnings)
    else:
        spec = spectrum(date=windsoln.spec_date, wl_norm=wl_norm,print_warning=print_warnings)
    for name in windsoln.species_list:
        spec.add_species(name)
    soln_resolved = windsoln.spec_resolved/wl_norm
    spec.set_resolved(*soln_resolved)
    soln_normalized = windsoln.spec_normalized/wl_norm
    spec.set_normalized(*soln_normalized)
    soln_window = windsoln.spec_window/wl_norm
    spec.set_window(*soln_window, kind=windsoln.spec_kind)
    return spec.plot(var,xaxis,semimajor_au,highlight_euv)


def spectrum_binning_plot(windsoln, var='F_wl',xaxis='wl',semimajor_au=1.0, plot_polys=False, 
                          wl_norm=1e-7,print_warnings=False):
    """Plot the spectrum. Displays the observations, smoothed, and spans.
    
    *Legend Key:*
    Subbin edges -  are automatically set at ionization edges (including K-shell ionization edges) 
    for species present in a given simulation for maximum accuracy in calculating ionization rates.
    Crits -  dashed lines are the critical points in the smoothing 
    Bin edges - wavelength range window edges

    Args:
        var (str): Which variable plotted, energy ('F_wl') or number ('Phi_wl')
        xaxis (str): 'wl' or 'energy'; x variable. Wavelength in nm or energy in eV
        semimajor_au (float): semimajor axis in units of au to which to scale the spectrum
        highlight_euv (bool): default=True; highlights EUV and XUV range and prints
                        APPROXIMATE fluxes in ergs/s/cm2 each range.
        wl_norm (float): default=1e-7; wavelength normalization factor to convert nm to cm.
        print_warnings (bool): default=False; prints any warnings generated when loading spectrum.

    Returns:
        (fig, ax): The figure and axis object plotted on
    """
    if windsoln.spec_src_file != 'scaled-solar':
        spec = spectrum(lisird=False,
                        spectrum_file=windsoln.spec_src_file,
                        wl_norm=wl_norm,print_warning=print_warnings)
    else:
        spec = spectrum(date=windsoln.spec_date, wl_norm=wl_norm,print_warning=print_warnings)
    for name in windsoln.species_list:
        spec.add_species(name)
    soln_resolved = windsoln.spec_resolved/wl_norm
    spec.set_resolved(*soln_resolved)
    soln_normalized = windsoln.spec_normalized/wl_norm
    spec.set_normalized(*soln_normalized)
    soln_window = windsoln.spec_window/wl_norm
    spec.set_window(*soln_window, kind=windsoln.spec_kind)    
    return spec.binning_plot(var=var,xaxis=xaxis,semimajor_au=semimajor_au,plot_polys=plot_polys)
