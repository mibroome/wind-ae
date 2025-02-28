import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from cycler import cycler

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
                    pvar=[['v','rho'],['T','Ys']], norm=[[1e5,1e0],[1e0,1e0]],
                    radius_prefix=None, 
                    sub_sonic=False, past_rmin=False,
                    sonic_vert=True):
    """
    Description:
        Produces a velocity, density, temperature, and neutral fraction plot
        from a given windsoln.
        
    Arguments:
        windsoln: windsoln object to be plotted
        ax: the ax array on which plots are placed (must have shape = (2,2))
        
    Keyword arguments:
        ax_Ys: ax for plotting ionization fraction on (e.g. ax_Ys=ax[0][1].twinx())
        label: label to be added to ax[0][0] (use fig.legend() to display)
        alpha: alpha of lines plotted
        pvar: variables to plot
        norm: values to norm (divide) pvars with
        radius_prefix: SI prefix for x-axis (None defaults to units of Rp)
                       (assuming windsoln.soln in cgs)
        sub_sonic: only plot sub-sonic region (boolean)
        past_rmin: only plot past Rmin region (boolean)
        sonic_vert: add verticle line at sonic point (boolean)
        
    Returns:
        Nothing (but plots windsoln in supplied axis array)
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
                    ax[i][j].set_xlabel(r'Radius ($R_\mathrm{p}$)') 
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
            ax[i][j].set_ylim([min(ylims[0], 0.9*ylo/norm[i][j]),
                               max(ylims[1], 1.1*yhi/norm[i][j])])
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
            elif pvar[i][j] == 'Ys':
                ax[i][j].set_ylabel(r'Neutral Fraction')
                ax[i][j].set_yscale('log')
            else:
                print("WARNING: pvar {:s} has no preset labeling."
                      .format(pvar[i][j]))
            if j == len(ax[0])-1:
                ax[i][j].set_ylabel(ax[i][j].yaxis.get_label_text(),
                                    rotation=270, va="bottom")
    return

def FourPlot(soln,Mdot_legend=True,line_color='k',line_style='-',first_of_many=False,ax=0): 
    '''
    Description: 
        Plots density (g/cm3), temperature (K), velocity (10 km/s), 
        ionization fraction as a function of r (Rp).
        
    Arguments:
        soln - windsoln object (sim.windsoln)
        Mdot_legend - Bool; if True, put Mdot in legend of plot. Else, just prints.
        line_color - str; line color 
        line_style - str; line style
        first_of_many - Bool; True if this the first of many Four plots 
                        to be plotted on the same axes. 
    Returns:
        ax - axes object (if first_of_many=True)
        title
    '''
    radius_norm = 1.0
    alpha=0.5
    radius = soln.soln_norm['r']
    minrad = float(radius[soln.soln['T'] == np.min(soln.soln['T'])])
    R_H = soln.semimajor*(soln.Mp/(3*soln.Mstar))**(1/3) / soln.Rp
    
    custom_cycler = (cycler(linestyle=['-', '--', ':', '-.',(0, (1, 10))]))

    stack=2
    if first_of_many==True:
        fig, ax = plt.subplots(stack,2,sharex=True,figsize=[11,6])
        fig.subplots_adjust(hspace=0)
        ax[0,0].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                     c=line_color, alpha=alpha, ls='--', zorder=3,label='Sonic Point')
        ax[0,0].axvline(R_H,ls=':',c=line_color,label='Hill Radius')
    #Density (+Hill Radius and sonic point for legend purposes)
    ax[0,0].semilogy(soln.soln_norm['r'], soln.soln['rho'],c=line_color)
    ax[0,0].legend()
    ax[0,0].set_ylabel(r'Density (g/cm$^3$)')
    #Velocity 
    ax[1,0].plot(soln.soln_norm['r'], soln.soln['v']/(10*100*1000),c=line_color)
    ax[1,0].set_ylabel(r'Velocity (10 km/s)')
    #Temperature
    mdot = (soln.Mdot)#/const.Msun)*3.154e+7 #in Msun/year
    print('*******Mdot = %.2e g/s*********' %mdot)
    ax[0,1].semilogy(soln.soln_norm['r'], soln.soln['T'],c=line_color,label=r'$\dot{M}$ = %.1e g/s' %(mdot))
    ax[0,1].set_ylabel(r'Temperature (K)')
    if Mdot_legend==True:
        ax[0,1].legend()
#   a  ax[2,0].set_ylim((5e2,5e4))
    #Ionization fraction
    title = ''
    ax[1,1].set_prop_cycle(custom_cycler)
    for j,spname in enumerate(soln.species_list):
        spname = spname.replace(' ','')
        ax[1,1].semilogy(soln.soln_norm['r'],1-soln.soln['Ys_'+spname],label=spname,c=line_color)
        title += spname+': %.2f, '%(soln.HX[j])
    ax[1,1].set_ylabel(r'Ionization Fraction')
    ax[1,1].set_ylim((1e-3,0.99999999))
    ax[1,1].set_xticks([1e-3,1e-2,1e-1],labels=[r'$10^{\text{-}3}$',r'$10^{-2}$',r'$10^{-1}$'])
    ax[1,1].legend()
    #Plotting shared lines
    for k in range(2):
        ax[1,k].set_xlabel(r'Radius (R$_p$)')
        for m in range(stack):
            ax[m,k].axvline(radius_norm*(soln.soln_norm['z'][1]+1.),
                             c=line_color, alpha=alpha, ls='--', zorder=3)
            ax[m,k].axvline(R_H,ls=':',c=line_color)
    plt.gca().set_xlim(left=soln.Rmin)
    if first_of_many==True:
        return ax,title


def energy_plot(windsoln, ax, alpha=0.8, all_terms=False, 
                CII_line_cool=False, CIII_line_cool=False, OII_line_cool=False,OIII_line_cool=False, 
                legend=True,sub_sonic=True):
    ax.plot(windsoln.soln_norm['r'], windsoln.soln['heat_ion'], '-',
            alpha=alpha, c='maroon',label='ionization heating')
    ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_lyman'], '--',
            alpha=alpha, c='tab:cyan', label='Lyman-alpha cooling')
    ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_PdV'], '--',
            alpha=alpha, c='tab:blue',label='PdV cooling')
    ax.plot(windsoln.soln_norm['r'], windsoln.soln['boloheat'], ls=(0, (3, 1, 1, 1, 1, 1)),
        alpha=alpha, c='tab:red',label='bolometric heating')
    ax.plot(windsoln.soln_norm['r'], -windsoln.soln['bolocool'], ls=(0, (3, 1, 1, 1, 1, 1)),
        alpha=alpha, c='lightblue',label='bolometric cooling')
    ax.plot(windsoln.soln_norm['r'], windsoln.soln['heat_advect'], '-.',
            alpha=alpha, c='tab:orange', label='advective heating')
    ax.plot(windsoln.soln_norm['r'], -windsoln.soln['heat_advect'], '-.',
            alpha=alpha, c='tab:gray', label='advective cooling')
    if all_terms:
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_rec'], ':',
            alpha=alpha, c='slateblue',label='recombination cooling')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_cond'], ':',
            alpha=alpha, c='navy',label='conductive cooling')
        ax.plot(windsoln.soln_norm['r'], windsoln.soln['cool_cond'], ':',
            alpha=alpha, c='lightsalmon',label='conductive heating')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_free'], ':',
            alpha=alpha, c='paleturquoise',label='free-free cooling')
    if CII_line_cool:
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_CII_1570000A'], '-',
            alpha=0.5, label=r'CII 15700000$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_CII_2326A'], '-',
            alpha=0.5, label=r'CII 2326$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_CII_1334A'], '-',
            alpha=0.5, label=r'CII 1334$\mathrm{\AA}$')
    if CIII_line_cool:
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_CIII_1910A'], '-',
            alpha=0.5, label=r'CIII 1910$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_CIII_977A'], '-',
            alpha=0.5, label=r'CII 977$\mathrm{\AA}$')
    if OII_line_cool:
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OII_834A'], '-',
            alpha=0.5, label=r'OII 834$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OII_2741A'], '-',
            alpha=0.5, label=r'OII 2741$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OII_3727A'], '-',
            alpha=0.5, label=r'OII 3727$\mathrm{\AA}$')    
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OII_7320A'], '-',
            alpha=0.5, label=r'OII 7320$\mathrm{\AA}$') 
    if OIII_line_cool:
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OIII_520000A'], '-',
            alpha=0.5, label=r'OIII 520000$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OIII_5000A'], '-',
            alpha=0.5, label=r'OIII 5000$\mathrm{\AA}$')
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OIII_166A'], '-',
            alpha=0.5, label=r'OIII 166$\mathrm{\AA}$')    
        ax.plot(windsoln.soln_norm['r'], -windsoln.soln['cool_OIII_84A'], '-',
            alpha=0.5, label=r'OIII 84$\mathrm{\AA}$') 
        
    peak = 10**np.ceil(np.log10(windsoln.soln['heat_ion'].max()))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(r'Energy Rates (erg s$^{-1}$ cm$^{-3}$)')
    ax.set_xlabel(r'Radius ($R_p$)')
    xlo, xhi = 0.9, windsoln.Rmax
    if sub_sonic or not windsoln.integrate_outward:
        xhi = windsoln.r_crit
    ax.set_xlim([xlo, xhi])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(ticker.ScalarFormatter())
    ax.set_ylim([peak/1e6, peak])
    if legend:
        ax.legend()
    return


def coriolis_plot(windsoln, ax, var='position', label_system=None):
    s_arr = np.linspace(windsoln.r_crit, windsoln.Rmax, 1000, endpoint=False)
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
                  labels+[r'$x$', r'$y$', r'$r$', r'$s$'])
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
        ax.set_ylabel('velocity (km/s)')
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
    ax.axvline(windsoln.r_cori, c=cori_vert_color, ls='--', zorder=0)
    return


def ballistic_plot(windsoln, ax, show_star=False, t_end=None):
    if t_end is None:
        t_end = np.log10(windsoln.stream_time(windsoln.r_cori)
                         -windsoln.stream_time(windsoln.r_crit))
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
    s = np.linspace(windsoln.r_crit, windsoln.Rmax, 10000)
    ax.plot(windsoln.cori_pos(s)[0], windsoln.cori_pos(s)[1], 'b-')
    # Overlay length scales
    circle1 = plt.Circle((0, 0), windsoln.r_crit, color='r', lw=5, fill=False,
                         zorder=0, alpha=0.8)
    circle2 = plt.Circle((0, 0), windsoln.r_cori, color='m', lw=5, fill=False,
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
        ax.axvline(windsoln.r_crit, ls='--', c='r', zorder=0)
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