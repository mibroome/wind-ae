#!/usr/bin/env python3

import os
import numpy as np
from ftplib import FTP
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy import interpolate, integrate

from wind_ae.McAstro.stars.spectrum.black_body import BlackBody, Planck_irradiance_wl

class PHOENIX():
    """
    Description:
        Class for analysising steller PHOENIX model, uses the Gottingen
        Spectral Library (GSL).
    
    Keyword arguments:
        Teff: stellar effective temperature (in Kelvins)
        log10g: log base 10 of stellar surface gravity
        Z: stellar metalicity (in [Fe/H])
        alpha: stellar metalicity of named metals (in [alpha/Fe])
               Named metals: O, Ne, Mg, Si, S, Ar, Ca, and Ti
               Note: [alpha/H] = [Fe/H] + [alpha/Fe]
        res: GSL spectrum resolution ('high', 'medium')
        verbose: debugging output (boolean)
        
    Methods:
        spec_plot: plot the stellar spectrum (spectral radiance vs
                                              wavelength)
        spec_integral:
        spec_wgt_integral:
        
    Source paper:
        Husser et al. 2013 (2013A&A...553A...6H)
    """
    fits_directory = os.path.dirname(os.path.abspath(__file__))+'/'
    wl_fits = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'

    def __init__(self, Teff=10170, log10g=4.09, Z=-0.03, alpha=0.0,
                 folder='GSL_spectrum_fits/', res='high', verbose=False):
        """ 
        Default:
            based off Kelt-9 (HD 195689) from Gaudi et al. 2017
            (2017Natur.546..514G)
                Teff     = 10170 (+/-450)
                log10(g) = np.log10(G*2.51*Msun/(2.361*Rsun)**2)
                         = 4.093 (+/-0.014)
                Z        = -0.03 (+/- 0.20)
                alpha    = 0.0
            Similarish to Fossati et al. 2018 (2018ApJ...868L..30F)
        """
        self.fits_directory += folder
        if self.fits_directory[-1] != '/':
            self.fits_directory += '/'
        self.wl_fits = self.fits_directory+self.wl_fits
        self.res = res
        self.verbose = verbose
        self.Teff = Teff
        self.log10g = log10g
        self.Z = Z
        self.alpha = alpha
        if self._range_check():
            print("ERROR: FAILED TO FULLY INITALIZE PHOENIX OBJECT")
            return
        self._round_input()
        self.PHOENIX_fits = None
        self._fetch_PHOENIX_fit()

        if self.res == 'medium':
            self.wl = np.linspace(3, 10, 70000)*1e-5
        else:
            with fits.open(self.wl_fits) as hdul:
                self.wl_hdr = hdul[0].header
                if self.wl_hdr['UNIT'] == 'Angstrom':
                    self.wl = hdul[0].data*1e-8
                    self.wl_hdr['UNIT'] = 'Centimeter'
        with fits.open(self.PHOENIX_fits) as hdul:
            self.spec_hdr = hdul[0].header
            self.spec = hdul[0].data
            if len(hdul) > 1:
                self.abun_hdr = hdul[1].header
                self.abun = hdul[1].data
        self.black_body = BlackBody(self.spec_hdr['PHXTEFF'])
        # Fit a black body to the EUV spectrum
        H_ionizing = self.wl < 912e-8
        if any(H_ionizing):
            popt, pcov = curve_fit(Planck_irradiance_wl, self.wl[H_ionizing],
                                   self.spec[H_ionizing],
                                   p0=[self.spec_hdr['PHXTEFF']])
            self.EUV_Teff = popt[0]
            self.EUV_black_body = BlackBody(self.EUV_Teff)

    def _range_check(self):
        """
        Description:
            Check that user supplied variables are within the GSL parameter
            space.
            
        Returns:
            if parameters are within GSL bounds (boolean)
        """
        # Z check
        Z_lo = -4.0
        Z_hi = 1.0
        if not (Z_lo <= self.Z <= Z_hi):
            print(r'Your Z (Z={:}) is outside the Gottingen Spectral '
                  'Library\'s Z range ({:} <= Z <= {:})'.format(self.Z,
                                                                Z_lo, Z_hi))
            return(True)
        # alpha check
        alpha_lo = -0.2
        alpha_hi = 1.2
        if not (alpha_lo <= self.alpha <= alpha_hi):
            print(r'Your $\alpha$ ($\alpha$={:}) is outside the Gottingen '
                  r'Spectral Library $\alpha$ range ({:} <= $\alpha$ <= {:})'
                  .format(self.alpha, alpha_lo, alpha_hi))
            return(True)
        # T_eff check
        Teff_lo = 2300
        Teff_hi = 15000
        if not (Teff_lo <= self.Teff <= Teff_hi):
            print(r'Your Teff (Teff={:}) is outside the Gottingen Spectral '
                  'Library\'s Teff range ({:} <= Teff <= {:})'
                  .format(self.Teff, Teff_lo, Teff_hi))
            return(True)
        # log_10(g) check
        log10g_lo = 0.0
        log10g_hi = 6.0
        if not (log10g_lo <= self.log10g <= log10g_hi):
            print(r'Your log10(g) (log10(g)={:}) is outside the Gottingen '
                  'Spectral Library\'s log10(g) range '
                  '({:} <= log10(g) <= {:})'.format(self.log10g,
                                                    log10g_lo, log10g_hi))
            return(True)
        return(False)

    def _round_input(self):
        """
        Description:
            Since the GSL is a grid of simulations we need to pick select
            the relevant grid points.
        
        Notes:
            ^Currently we round to nearest grid point
            ^In the future one could contemplate performing 4D interpolation
        """
        # Round T_eff
        if self.Teff <= 7000:
            round_Teff = 100*int(self.Teff/100+0.5)
        elif self.Teff <= 12000:
            round_Teff = 200*int(self.Teff/200+0.5)
        else:
            round_Teff = 500*int(self.Teff/500+0.5)
        if self.verbose:
            print("Rounding Teff from {:d} K to {:d} K".format(
                self.Teff, round_Teff))
        self.Teff = round_Teff
        # Round log_10(g)
        round_log10g = 0.5*int(2*self.log10g+0.5)
        if self.verbose:
            print("Rounding log10g from {:.1F} to {:.1F}".format(
                self.log10g, round_log10g))
        self.log10g = round_log10g
        # Round Z
        if self.Z <= -2.0:
            round_Z = -1*int(abs(self.Z)+0.5)
        else:
            round_Z = ((self.Z > 0)-(self.Z < 0))*0.5*int(2*abs(self.Z)+0.5)
        if round_Z == 0:
            round_Z = -0.0
        if self.verbose:
            print("Rounding Z from {:.1F} to {:.1F}".format(self.Z, round_Z))
        self.Z = round_Z
        # Round alpha
        round_alpha = 0.2*int(5*self.alpha+0.5)
        if self.verbose:
            print("Rounding alpha from {:.1F} to {:.1F}".format(
                self.alpha, round_alpha))
        self.alpha = round_alpha
        return

    def _fetch_PHOENIX_fit(self):
        # Make sure fits_directory exis
        if not os.path.isdir(self.fits_directory):
            os.makedirs(self.fits_directory)
        # Compose PHOENIX file name
        metal = '{:+.1F}'.format(self.Z)
        if self.alpha != 0:
            metal += '.Alpha={:+.2F}'.format(self.alpha)
        file = 'lte{:05d}-{:.2F}{:s}.'.format(
            int(self.Teff), self.log10g, metal)
        file += 'PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
        # Check that we have wavelength fits file
        if os.path.isfile(self.wl_fits) or self.res == 'medium':
            if self.verbose:
                print("We already fetched file.")
        else:
            ftp = FTP('phoenix.astro.physik.uni-goettingen.de')
            welcome = ftp.login()
            ftp.cwd('v2.0/HiResFITS/')
            wl_file = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
            if wl_file in ftp.nlst():
                if self.verbose:
                    print("Found {:s}".format(wl_file))
                    print("Downloading...")
                ftp.retrbinary("RETR " + wl_file,
                               open(self.fits_directory+wl_file,
                                    'wb').write)
                if self.verbose:
                    print("Finished")
            else:
                print('Did not find {:s} file in the Gottingen Spectral '
                      'Library'.format(wl_file))
                print(ftp.nlst())
            goodbye = ftp.quit()

        # Check if we already have file
        if os.path.isfile(self.fits_directory+file):
            if self.verbose:
                print("We already fetched file.")
        else:
            # Fetch fits form Gottingen Spectral Library ftp server
            ftp = FTP('phoenix.astro.physik.uni-goettingen.de')
            welcome = ftp.login()
            ftp.cwd('v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011')
            metal = 'Z'+metal
            if metal not in ftp.nlst():
                print('Did not find {:s} directory in the Gottingen Spectral '
                      'Library'.format(metal))
            else:
                ftp.cwd(metal)
                if file in ftp.nlst():
                    if self.verbose:
                        print("Found {:s}".format(file))
                        print("Downloading...")
                    ftp.retrbinary("RETR " + file,
                                   open(self.fits_directory+file, 'wb').write)
                    if self.verbose:
                        print("Finished")
                else:
                    print('Did not find {:s} file in the Gottingen Spectral '
                          'Library'.format(file))
            goodbye = ftp.quit()
        self.PHOENIX_fits = self.fits_directory+file

    def spec_plot(self, ax, **kwargs):
        ax.plot(self.wl, self.spec, **kwargs)

    def spec_integral(self, lob=1e-8, upb=912e-8):
        integrand = np.where((self.wl <= upb) & (self.wl >= lob))[0]
        result = np.trapz(self.spec[integrand], x=self.wl[integrand])
        if upb > self.wl[-1]:
            if self.verbose:
                print("Warning: Extending spectrum past data with black body")
            result += self.black_body.spec_integral(lob=self.wl[-1], upb=upb)
        if lob < self.wl[0]:
            if self.verbose:
                print("Warning: Extending spectrum past data with EUV black "
                      "body")
            result += self.EUV_black_body.spec_integral(lob=lob,
                                                         upb=self.wl[0])
        return result

    def spec_wgt_integral(self, x_wgt, wgt, lob=1e-8, upb=912e-8):
        wgt_lob = x_wgt[0]
        wgt_upb = x_wgt[-1]
        if lob < wgt_lob or upb > wgt_upb:
            print('ERROR: Will not interpolate past x_wgt bounds.')
            print('lob:{}, wgt_lob:{}, upb:{}, wgt_upg:{}'
                  .format(lob, wgt_lob, upb, wgt_upb))
            return
        f_wgt = interpolate.interp1d(x_wgt, wgt)
        integrand = np.where((self.wl <= upb) & (self.wl >= lob))[0]
        wgt = f_wgt(self.wl[integrand])
        return np.trapz(wgt*self.spec[integrand], x=self.wl[integrand])
