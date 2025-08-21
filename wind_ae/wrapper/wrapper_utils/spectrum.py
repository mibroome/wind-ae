#!/bin/env python3

import copy
import numpy as np
import matplotlib.pyplot as plt

from . import constants as const
from wind_ae.McAstro.planets.insolation.glq_spectrum import glq_spectrum
import wind_ae.McAstro.atoms.atomic_species as McAtom
import importlib.resources as pkg_resources


class spectrum:
    path = str(pkg_resources.files('wind_ae'))+'/'
    def __init__(self, lisird=True, date='2009-01-01', spectrum_file='', wl_norm=1e-7, print_warning=False):
        """
        Initialize a Spectrum object for solar or user-defined spectra.

        Args:
            lisird (bool, optional): If True, use a scaled FISM2 `LISIRD <https://lasp.colorado.edu/lisird/data/fism_daily_hr>`__ solar spectrum from the selected date. If False, use user input. Defaults to True.
            date (str, optional): Date of solar observation (format: '%YYYY-%mm-%dd'). Defaults to '2009-01-01'.
            spectrum_file (str, optional): If lisird is False, path to CSV file of user spectrum. File should contain headers: 'wl' (cm), 'F_wl' (erg/s/cm^2/cm), 'unc', 'nu' (1/s), 'F_nu'. Defaults to ''.
            wl_norm (float, optional): Converts cm to another wavelength unit via division. Defaults to 1e-7 (nanometers).
            print_warning (bool, optional): If True, print warnings. Defaults to False.

        Attributes:
            data (np.ndarray): Spectrum data array.
            data_norm (np.ndarray): Normalized spectrum data array.
            kind (str): Type of spectrum ('full').
            norm_span (np.ndarray): Normalized wavelength span.
            rslv_span (np.ndarray): Resolved wavelength span.
            spectrum_file (str): Path to spectrum file (if user input).
            kshell_ionpots (dict): K-shell ionization potentials for select elements.

        """
        #Scaled and smoothed solar spectrum from LISIRD mission
        if lisird == True:
            self.date = date
            self.wl_norm = wl_norm
            self.glq_spectrum = glq_spectrum(date=date, wl_norm=wl_norm)
        #or smoothed user-input spectrum
        else:
            self.date = date
            self.wl_norm = wl_norm
            self.glq_spectrum = glq_spectrum(lisird=False, filename=spectrum_file, wl_norm=wl_norm)
        self.data = self.glq_spectrum.data
        self.data_norm = self.glq_spectrum.data_norm
        self.kind = 'full'
        self.norm_span = np.array([self.glq_spectrum.spectrum.wl_min,
                                   self.glq_spectrum.spectrum.wl_max])
        # Set resolve before set normalized (as smoothing requires resolved)
        self.set_resolved(*self.norm_span)
        self.set_normalized(*self.norm_span)
        self.set_window(*self.rslv_span)
        self.spectrum_file = spectrum_file
        self.kshell_ionpots = { 'C':308.67523704,  'N':426.57706307,  'O':563.23652347, 
                    'Mg':1336.8046517 , 'Si':1877.59068171, 'S':2512.3371317 }
        return


    def set_resolved(self, lob, upb, match_norm=False):
        """Wavelength of upper and lower boundaries for the desired range of the input spectrum. 
        Values are in NANOMETERS. Suggested to set to same as window.
        
        Args:
            lob (float): Lower boundary of the window (in NANOMETERS)
            upb (float): Upper boundary of the window (in NANOMETERS)

        Returns:
            None
        """
#         if (lob*self.wl_norm) < 1e-4 or upb > 1e4:
#             print("WARNING: Bounds should be given in .")
        self.rslv_span = np.asarray([lob, upb])
        return
        # # Error check
        # if lob < self.norm_span[0] or upb > self.norm_span[1]:
        #     print('ERROR: Resolved should be a subspace of the normalization.\n'
        #           f'       {lob:.5g} < {self.norm_span[0]:.5g} or {upb:.5g} > '
        #           f'{self.norm_span[1]:.5g}.\n'
        #           '       Adjust normalized to contain the requested span.')
        #     return
        # # Set resolved span and normalize (calculate resolved quantities)
        # self.rslv_span = np.asarray([lob, upb])
        # self.glq_spectrum.normalize(self.norm_span, rslv_span=self.rslv_span)
        # # Update dataframes after normalizing
        # self.data = self.glq_spectrum.data
        # self.data_norm = self.glq_spectrum.data_norm
        # return


    def set_normalized(self, lob, upb):
        """Wavelength of upper and lower boundaries for the desired range of the input spectrum. 
        Values are in NANOMETERS. Suggested to set to same as window.
        
        Args:
            lob (float): Lower boundary of the window (in NANOMETERS)
            upb (float): Upper boundary of the window (in NANOMETERS)

        Returns:
            None
        """
        # Error check
        if lob > self.rslv_span[0] or upb < self.rslv_span[1]:
            print('ERROR: Normalized should contain the resolved.\n'
                  '       Adjust resolved to be a subspace of requested span.')
            return
        # Set normalized span, which define the smoothed domain
        self.norm_span = np.asarray([lob, upb])
#         print("spectrum.py: print norm_span",self.norm_span)
        self.glq_spectrum.bin_breaks = self.norm_span
        self.glq_spectrum.truncate_spectrum(*self.norm_span)
        # Normalize the spectrum and smooth
        self.glq_spectrum.normalize(self.norm_span, rslv_span=self.rslv_span)
        self.glq_spectrum.smooth_spectrum(crit_dist=.1e-7/self.wl_norm,
                                          conserve_phi=True)
        # Subbin spectrum based on fitting a polynomial
        self.glq_spectrum.fit_polynomial(lsq_err=1e-5, crits_enclosed=3)
        # Update dataframes after normalizing and smoothing
        self.data = self.glq_spectrum.data
        self.data_norm = self.glq_spectrum.data_norm
        return


    def set_window(self, lob, upb=None, kind='full', verbose=True):
        """Wavelength of upper and lower boundaries for the desired range of the input spectrum. 
        Window values are in NANOMETERS.

        Args:
            lob (float): Lower boundary of the window (in NANOMETERS)
            upb (float): Upper boundary of the window (in NANOMETERS)
            kind (str): Kind of window ('full', 'fixed', etc.).
            verbose (bool): Whether to print verbose output.
        """
        if upb is None or lob == upb:
            if kind == 'full':
                print("ERROR: Kind cannot be full for monochromatic source.")
                return
            elif kind == 'fixed':
                self.mono_wl = lob
            else:
                self.mono_wl = np.zeros(self.glq_spectrum.n_bins)
                for b in range(self.glq_spectrum.n_bins):
                    mean_E = (np.nansum(self.glq_spectrum.F_smth_rslv[b])
                              /np.nansum(self.glq_spectrum.Phi_smth_rslv[b]))
                    self.mono_wl[b] = const.hc/(mean_E*self.wl_norm)
            self.wndw_span = np.asarray([self.mono_wl[0], self.mono_wl[0]])
        else:
            if kind != 'full':
                print("ERROR: Kind must be full when window has non-zero span.")
                return
            self.wndw_span = np.asarray([lob, upb])

            if (verbose and
                np.any(abs(self.wndw_span-self.rslv_span)/self.wndw_span>=
                1e-5)):
                print('WARNING: Window and resolved do not match.\n'
                      f'         Window: [{self.wndw_span[0]:.5g}, '
                      f'{self.wndw_span[1]:.5g}] != [{self.rslv_span[0]:.5g}, '
                      f'{self.rslv_span[1]:.5g}]: resolved')
            self.mono_wl = (
                const.hc/(self.glq_spectrum.F_tot/self.glq_spectrum.Phi_tot)
                /self.wl_norm
            )
        self.glq_spectrum.set_abscissas(self.wndw_span, trans_poly_deg=2)
        self.kind = kind
        return


    def add_species(self, species):
        """Adds species to spectrum smoothing algorithm in McAstro/planets/insolation/glq_spectrum.py.

        Args:
            species (str or list): The species to add to the spectrum smoothing.
        """
        if type(species) is str:
            spaced_species = McAtom.formatting_species_list([species])[0]
            self.glq_spectrum.add_species(spaced_species)
            if spaced_species[0] in self.kshell_ionpots:
                self.glq_spectrum.add_species(spaced_species,kshell=True)
        else:
            for s in species:
                spaced_species = McAtom.formatting_species_list([s])[0]
                self.glq_spectrum.add_species(spaced_species)
                if spaced_species[0] in self.kshell_ionpots:
                    print("triggered")
                    self.glq_spectrum.add_species(spaced_species,kshell=True)
        self.species_list = self.glq_spectrum.species_list
        return


    def remove_species(self, species):
        """Removes species from spectrum smoothing algorithm in McAstro/planets/insolation/glq_spectrum.py.

        Args:
            species (str or list): The species to remove from the spectrum smoothing. 
        """
        if type(species) is str:
            species = McAtom.formatting_species_list([species]).replace(' ','')
            self.glq_spectrum.remove_species(species)
        else:
            for s in species:
                s = McAtom.formatting_species_list([s]).replace(' ','')
                self.glq_spectrum.remove_species(s)
        self.species_list = self.glq_spectrum.species_list
        return


    def change_date(self, date):
        """DEPRECATED: Unlikely to need. Changes date if using solar FISM2 spectra downloaded from the Lisird database.
        Wind-AE no longer scrapes the database for spectra."""
        # Error check
        if date == self.date:
            print(f'WARNING: {date} is already the date loaded.\n'
                  '         Returning without alterations.')
            return
        # Generate new glq_spectrum
        self.glq_spectrum = glq_spectrum(date=date, wl_norm=self.wl_norm)
        # Give the new spectrum the user defined properities of the old specturm
        for species in self.species_list:
            self.add_species(species.atomic_data.name)
        # Set resolve before set normalized (as smoothing requires resolved)
        self.set_resolved(*self.norm_span)
        self.set_normalized(*self.norm_span)
        self.set_window(*self.rslv_span)
        return


    def generate(self, kind, savefile=path+'inputs/spectrum.inp'):
        """Generates the spectrum file in wind_ae/inputs/spectrum.inp.
        Will not write newly generated spectrum unless this function is called.

        Args:
            kind (str): The kind of spectrum to generate (e.g., 'full', 'mono', etc.) 
            savefile (str): Should be wind_ae/inputs/spectrum.inp
        """
        self.table = self.glq_spectrum.write_csv(savefile, kind=kind,mono_wl=None)
        return


    def binning_plot(self,var='F_wl',xaxis='wl',semimajor_au=1.0, plot_polys=False):
        """Plot the spectrum with binning lines for the spectrum loaded into the spectrum() object.

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

        Returns:
            (fig, ax): The figure and axis object plotted on
        """

        self.glq_spectrum.plot(var,xaxis,semimajor_au,plot_polys)
        return
        
    def plot(self, var='F_wl',xaxis='wl',semimajor_au=1.0,highlight_euv=True):
        """Plot the spectrum. Displays the observations, smoothed, and spans.  

        Arguments:
            var (str): Which variable plotted, energy ('F_wl') or number ('Phi_wl')
            xaxis (str): 'wl' or 'energy'; x variable. Wavelength in nm or energy in eV 
            semimajor_au (float): semimajor axis in units of au to which to scale the spectrum
            highlight_euv (bool): default=True; highlights EUV and XUV range and prints 
                           APPROXIMATE fluxes in ergs/s/cm2 each range.

        Returns:
            (fig, ax): The figure and axis object plotted on
        """
        fig, ax = plt.subplots()
        mk = ((self.data_norm['wl']>=self.glq_spectrum.bin_breaks[0])
              &(self.data_norm['wl']<=self.glq_spectrum.bin_breaks[-1]))
        if var == 'F_wl':
            smth = 'f_wl_smth'
            norm = self.glq_spectrum.F_tot/self.glq_spectrum.wl_norm * (1/semimajor_au)**2
            ax.set_ylabel(f'Spectral irradiance at {semimajor_au:.2f} au\n'
                          r'($F_{\lambda}$) [erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')
        elif var == 'Phi_wl':
            smth = 'phi_wl_smth'
            norm = self.glq_spectrum.Phi_tot/self.glq_spectrum.wl_norm * (1/semimajor_au)**2
            ax.set_ylabel(f'Spectral photon irradiance at {semimajor_au:.2f} au\n'
                      r'($\Phi_{\lambda}$) [cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')
        else:
            print(f"ERROR: var: {var}, not recognized. Use 'F_wl' or 'Phi_wl.")
            return
        if xaxis=='wl':
            l1, = ax.plot(self.data_norm['wl'][mk], 
                          self.data[var][mk]*(1/semimajor_au)**2, lw=1,
                          label='Spectrum')
            l2, = ax.plot(self.data_norm['wl'][mk], norm*self.data_norm[smth][mk],
                          lw=2, label='Smoothed')
            v1 = ax.axvline(self.wndw_span[0], zorder=0, ls='--', c='k', lw=3,
                            label=f'Window: {self.wndw_span[0]:.2f} - {self.wndw_span[1]:.1f} nm ')
            
            if highlight_euv == True:
                wl_range = np.array(self.data_norm['wl'][mk])
                flux = norm*self.data_norm[smth][mk]
                delta_wl = wl_range[1]-wl_range[0]
                flux *= np.diff(wl_range,prepend=wl_range[0]-delta_wl)*1e-7
                max_F = max(self.data[var][mk]*(1/semimajor_au)**2)*0.8

                euv_range = [max(12.4,min(wl_range)),
                             max(wl_range)]
                ax.axvspan(euv_range[0],euv_range[1],color='tab:purple',alpha=0.1)
                if max(self.data_norm['wl'][mk])>91:
                    xuv_range = [min(wl_range),12.4]
                    ax.axvspan(xuv_range[0],xuv_range[1],color='c',alpha=0.1)
                    F_xuv = sum(flux[(wl_range>xuv_range[0]) & (wl_range<xuv_range[1])])
                    ax.text(xuv_range[0]+np.diff(xuv_range)/4,max_F,
                            r'F$_{X\rm{-}ray}\sim$%.0f'%F_xuv,color='darkcyan',
                            fontsize=14,weight='bold')
                F_euv = sum(flux[(wl_range>euv_range[0])&(wl_range<euv_range[1])])
                ax.text(euv_range[0]+np.diff(euv_range)/3,max_F,
                         r'F$_{EUV}\sim$%.0f'%F_euv,color='tab:purple',
                        fontsize=14,weight='bold')

            
            ax.axvline(self.wndw_span[1], zorder=0, ls='--', c='k', lw=3)
            ax.set_xlabel(r'Wavelength ($\lambda$) [nm]')
        
        if xaxis=='energy':
            convert = const.hc/(self.wl_norm*const.eV)
            l1, = ax.plot(convert/self.data_norm['wl'][mk], 
                          self.data[var][mk]*(1/semimajor_au)**2, lw=1,
                          label='Spectrum')
            l2, = ax.plot(convert/self.data_norm['wl'][mk],
                          norm*self.data_norm[smth][mk],
                          lw=2, label='Smoothed')
            v1 = ax.axvline(convert/self.wndw_span[0], zorder=0, ls='--', c='k', lw=3,
                            label=f'Window: {const.hc/(self.wndw_span[1]*self.glq_spectrum.wl_norm)/const.eV:.2f} - {const.hc/(self.wndw_span[0]*self.glq_spectrum.wl_norm)/const.eV:.0f} eV')
            
            if highlight_euv == True:
                wl_range = np.array(self.data_norm['wl'][mk])
                flux = norm*self.data_norm[smth][mk]
                delta_wl = wl_range[1]-wl_range[0]
                flux *= np.diff(wl_range,prepend=wl_range[0]-delta_wl)*1e-7
                wl_range = convert/wl_range #converting to eV
                max_F = max(self.data[var][mk]*(1/semimajor_au)**2)*0.8

                euv_range = [max(10,min(wl_range)),min(100,max(wl_range))]
                ax.axvspan(euv_range[0],euv_range[1],color='tab:purple',alpha=0.1)
                if max(self.data_norm['wl'][mk])>91:
                    xuv_range = [100,max(wl_range)]
                    ax.axvspan(xuv_range[0],xuv_range[1],color='c',alpha=0.1)
                    F_xuv = sum(flux[(wl_range>xuv_range[0]) & (wl_range<xuv_range[1])])
                    ax.text(xuv_range[0]+np.diff(xuv_range)/4,max_F,
                            r'F$_{Xray}\sim$%.0f'%F_xuv,color='darkcyan',
                            fontsize=14,weight='bold')
                F_euv = sum(flux[(wl_range>euv_range[0])&(wl_range<euv_range[1])])
                ax.text(euv_range[0]+np.diff(euv_range)/5,max_F,
                         r'F$_{EUV}\sim$%.0f'%F_euv,color='tab:purple',
                        fontsize=14,weight='bold')


            ax.axvline(convert/self.wndw_span[1], zorder=0, ls='--', c='k', lw=3)
            ax.set_xlabel(r'Energy ($E$) [eV]')
            ax.set_xticks([20,60],labels=['20','60'])

            ax.set_xscale('log')
        if self.spectrum_file == '':
            self.spectrum_file = f'Scaled_solar: {self.date}'
        ax.set_title(self.spectrum_file)
        ax.set_yscale('log')
        fig.tight_layout(pad=0.3)
        fig.subplots_adjust(bottom=0.3, top=0.9)
        lines = [l1, l2, v1]
        labels = [l.get_label() for l in lines]
        fig.legend(lines, labels, bbox_to_anchor=(0.5, -0.00), loc='lower center',
                   ncol=3, fontsize=14)
        return fig, ax