#!/usr/bin/env python3

import math, copy
import numpy as np
import pandas as pd
from scipy import integrate, interpolate
from scipy.signal import savgol_filter

from wind_ae.McAstro.utils import constants as const
from wind_ae.McAstro.stars.spectrum.lisird import lisird_spectrum
from wind_ae.McAstro.atoms.atomic_species import atomic_species


class species:
    def __init__(self):
        self.atomic_data = None
        self.Xfrac = None
        self.sigma = None
        self.I_ion = None

        
class glq_rates:
    """
    Description:
        Class for calculating the Gauss-Legendre quadrature dataponts of the
        ionization and heating rates of a stellar spectrum in an atmosphere
        comprised of the given absorbing species.
        
        Order of operation is to create a glq_rates object, then add the
        consitutients of the atmosphere with add_species(). Next the
        spectrum should be subbinned with subbin(). Following the max degree
        of the polynomial needed to approximate the rates is calculated by
        calling get_max_degree(). Then the GLQ datapoints are set by calling
        set_abscissas(). The data can be saved to a csv with write_csv().
    """
    def __init__(self, mission='fism2', date='2002-01-01'):
        # Stellar "Observation" info
        self.mission = mission
        self.date = date
        # Observed spectrum
        spectrum = lisird_spectrum(mission=mission, date=date, sort_values='wl')
        self.spectrum = spectrum
        self.data = spectrum.data
        self.data['Phi_wl'] = self.data['F_wl']/(const.hc/self.data['wl'])
        # Data in normalized units
        self.data_norm = self.data.copy()
        self.data_norm['wl'] = self.data_norm['wl']*1e7 
        self.data_norm['F_wl'] = self.data_norm['F_wl']*1e-7 
        # Composition of atmosphere
        self.n_species = 0
        self.species_list = []
        self.Xfrac_tot = 0.0
        # Binning of spectrum
        self.n_bins = 1
        self.bin_breaks = np.array([spectrum.wl_min, spectrum.wl_max])
        # Smoothing of spectrum
        self.F_polys, self.Phi_polys, self.log_F_polys, self.log_Phi_polys = (
            [None] for i in range(4)
        )
        self.F_poly_coeffs, self.Phi_poly_coeffs, self.rescale = (
            [None] for i in range(3)
        )
        self.log_F_poly_coeffs, self.log_Phi_poly_coeffs = (
            [None] for i in range(2)
        )
        # Gaussian Quadrature
        self.max_degree = None
        self.bin_absc, self.bin_wgts = ([None] for i in range(2))
        
        
    def bin_spectrum(self, bin_breaks):
        if type(bin_breaks) is float or len(bin_breaks) == 1: # monochromatic
            self.n_bins = 0
            if type(bin_breaks) is float:
                self.bin_breaks = [bin_breaks]
            else:
                self.bin_breaks = bin_breaks
            return
        # Enforce sorted list
        bin_breaks = np.array(bin_breaks)
        bin_breaks.sort()
        self.bin_breaks = bin_breaks
        self.n_bins = len(bin_breaks)-1
        # Check bins
        for ibin in range(self.n_bins):
            if (self.spectrum.wl_min > bin_breaks[ibin]
                or self.spectrum.wl_max < bin_breaks[ibin+1]):
                print("WARNING: Mission's spectrum does not span "
                      "requested bin's wavelengths.")
                print(f'        ({self.mission}) wavelegnth span: '
                      f'[{self.spectrum.wl_min*1e7:.2f}, '
                      f'{self.spectrum.wl_max*1e7:.2f}] nm')
                print('WARNING: Truncating binning to fit mission spectrum.')
                if bin_breaks[ibin+1] < self.spectrum.wl_min:
                    self.bin_breaks = self.bin_breaks[1:]
                    self.n_bins -= 1
                    continue
                elif bin_breaks[ibin] > self.spectrum.wl_max:
                    self.bin_breaks = self.bin_breaks[:ibin]
                    self.n_bins = len(self.bin_breaks)-1
                    break
                elif bin_breaks[ibin] < self.spectrum.wl_min:
                    self.bin_breaks[0] = self.spectrum.wl_min
                elif bin_breaks[ibin+1] > self.spectrum.wl_max:
                    dropped_bins = len(bin_breaks)-1-self.n_bins
                    last_bin = ibin+1-dropped_bins
                    self.bin_breaks[last_bin] = self.spectrum.wl_max
                    self.bin_breaks = self.bin_breaks[:last_bin+1]
                    self.n_bins = len(self.bin_breaks)-1
                    break
                else:
                    print('ERROR: This should be impossible... '
                          'Time for debugging!')
                    return
        self.F_polys, self.Phi_polys, self.log_F_polys, self.log_Phi_polys = (
            [None for i in range(self.n_bins)] for j in range(4)
        )
        self.F_poly_coeffs, self.Phi_poly_coeffs, self.rescale = (
            [None for i in range(self.n_bins)] for j in range(3)
        )
        self.log_F_poly_coeffs, self.log_Phi_poly_coeffs = (
            [None for i in range(self.n_bins)] for j in range(2)
        )
        
        
    def smooth_spectrum(self, window=49, degree=8, pdeg=3):
        """
        Description:
            Fits a polynomial to the spectrum. The spectrum of fism2 varies
            and is high resolution, making fitting polynomials difficult.
            First we smooth the peaks the troughs of the spectrum by passing
            the spectrum through a Savitzky-Golay filter. The filtered
            spectrum is then renormalized to conserve total energy in each
            bin. Next a 5th degree polynomial is fit to the filtered
            spectrum, again rescaling to preserver energy in each bin. The
            polynomial is calcualted by used a spline with with an infinite
            smoothing factor, which relaxes the spline to a single interval.
        
        Keyword arguments:
            window: datapoints used in Savitzky-Golay filter (default: 49)
            degree: Degree of Savitzky-Golay filter polynomial (default: 8)
            
        Returns:
            Nothing, but adds a filtered spectrum (F_savgol, Phi_savgol) and
            the a polynomial fit (F_poly, Phi_poly) to the object
        """
        if self.n_bins == 0: # monochromatic
            return
        wavelengths = self.data['wl']
        F_wl = self.data['F_wl']
        Phi_wl = self.data['Phi_wl']
        # pass spectrum through savgol filter
        F_hat = savgol_filter(np.log(F_wl), window, degree)
        Phi_hat = savgol_filter(np.log(Phi_wl), window, degree)
        F_tilda = copy.deepcopy(F_hat)
        Phi_tilda = copy.deepcopy(Phi_hat)
        F_Matilda = copy.deepcopy(F_hat)
        Phi_Matilda = copy.deepcopy(Phi_hat)
        F_expMatilda = copy.deepcopy(F_hat)
        Phi_expMatilda = copy.deepcopy(Phi_hat)
        for ibin in range(self.n_bins):
            # Masked data
            mk = ((self.data['wl']>=self.bin_breaks[ibin])
                    &(self.data['wl']<=self.bin_breaks[ibin+1]))
            wl = self.data['wl'][mk]
            F = self.data['F_wl'][mk]
            rescale = (integrate.simps(F, wl)
                       /integrate.simps(np.exp(F_hat[mk]), wl))
            F_tilda[mk] = np.exp(F_hat[mk])*rescale
            Phi_tilda[mk] = np.exp(Phi_hat[mk])*rescale
            # fit polynomial
            self.log_F_polys[ibin] = (
                interpolate.UnivariateSpline(wl, np.log(F_tilda[mk]),
                                             ext=2, k=5)
            )
            self.log_Phi_polys[ibin] = (
                interpolate.UnivariateSpline(wl, np.log(Phi_tilda[mk]),
                                             ext=2, k=5)
            )
            self.log_F_polys[ibin].set_smoothing_factor(np.inf)
            self.log_Phi_polys[ibin].set_smoothing_factor(np.inf)
            self.F_polys[ibin] = (
                interpolate.UnivariateSpline(wl*1e7, F_tilda[mk]*1e-7, ext=2, k=5)
            )
            self.Phi_polys[ibin] = (
                interpolate.UnivariateSpline(wl*1e7, Phi_tilda[mk]*1e-7, ext=2, k=5)
            )
            self.F_polys[ibin].set_smoothing_factor(np.inf)
            self.Phi_polys[ibin].set_smoothing_factor(np.inf)
            # np.poly
            self.F_poly_coeffs[ibin] = (
                np.polyfit(wl*1e7, F_tilda[mk]*1e-7, pdeg)
            )
            self.Phi_poly_coeffs[ibin] = (
                np.polyfit(wl*1e7, Phi_tilda[mk]*1e-7, pdeg)
            )
            self.log_F_poly_coeffs[ibin] = (
                np.polyfit(wl*1e7, np.log(F_tilda[mk]*1e-7), pdeg)
            )
            self.log_Phi_poly_coeffs[ibin] = (
                np.polyfit(wl*1e7, np.log(Phi_tilda[mk]*1e-7), pdeg)
            )
            # Rescale
            self.rescale[ibin] = (
                integrate.simps(F[mk], wl)
                /integrate.simps(np.exp(self.log_F_polys[ibin](wl)), wl)
            )
            F_Matilda[mk] = self.rescale[ibin]*self.F_polys[ibin](wl*1e7)
            Phi_Matilda[mk] = self.rescale[ibin]*self.Phi_polys[ibin](wl*1e7)
            F_expMatilda[mk] = np.exp(self.log_F_polys[ibin](wl))
            F_expMatilda[mk] *= self.rescale[ibin]
            Phi_expMatilda[mk] = np.exp(self.log_Phi_polys[ibin](wl))
            Phi_expMatilda[mk] *= self.rescale[ibin]
        self.data['F_savgol'] = F_tilda
        self.data['Phi_savgol'] = Phi_tilda
        self.data['F_poly'] = F_Matilda
        self.data['Phi_poly'] = Phi_Matilda
        self.data['F_exppoly'] = F_expMatilda
        self.data['Phi_exppoly'] = Phi_expMatilda
        
        
    def add_species(self, species_name, Xfrac=1.0,
                    adjust_bins=True, verbose=False):
        """
        Description:
            Add species to the atmosphere
            
        Keyword arguments:
            Xfrac: The number fraction of the species
            adjust_bins: Set a bin edge at the ionization edge (boolean)
            verbose: Report status of total Xfrac in atmosphere (boolean)
        """
        for species_obj in self.species_list:
            if species_obj.atomic_data.name == species_name:
                if species_obj.Xfrac == Xfrac:
                    print('WARNING: Species already in list, with same Xfrac.')
                    print('           Returning without doing anything.')
                    return
                else:
                    print('WARNING: Species already in list, '
                          'but with different Xfrac.\n'
                          '           Readjusting Xfrac...')
                    self.remove_species(species_name)
                    break
        new_atomic_data = atomic_species(species_name)
        if new_atomic_data.verner_data is None:
            return
        new_species = species()
        new_species.atomic_data = new_atomic_data
        new_species.sigma = (
            new_atomic_data.cross_section(const.hc/const.eV/self.data['wl'])
        )
        new_species.Xfrac = Xfrac
        new_species.I_ion = new_atomic_data.verner_data['E_th']*const.eV
        self.species_list.append(new_species)
        
        self.Xfrac_tot += Xfrac
        if verbose:
            self.check_Xfrac()
        if adjust_bins:
            new_bin_break = const.hc/new_species.I_ion
            if new_bin_break not in self.bin_breaks:
                self.bin_spectrum(list(self.bin_breaks)+[new_bin_break])
        self.n_species += 1
            

    def remove_species(self, species_name, adjust_bins=False, verbose=False):
        old_species = None
        for species_obj in self.species_list:
            if species_obj.atomic_data.name == species_name:
                old_species = species_obj
                continue
        if old_species is None:
            print('ERROR: Failed to find species in list')
            return
        self.Xfrac_tot -= old_species.Xfrac
        self.species_list.remove(old_species)
        if verbose:
            self.check_Xfrac()
        if adjust_bins:
            print('WARNING: Removing bin break associated with '
                  'ionization edge.')
            old_bin_break = const.hc/old_species.I_ion
            if old_bin_break in self.bin_breaks:
                self.bin_breaks = list(self.bin_breaks)
                self.bin_breaks.remove(old_bin_break)
                self.bin_spectrum(self.bin_breaks)
        self.n_species -= 1
        
        
    def normalize_Xfrac(self):
        for species_obj in self.species_list:
            species_obj.Xfrac /= self.Xfrac_tot
        self.Xfrac_tot /= self.Xfrac_tot 
        self.check_Xfrac(silent_pass=True)
        
        
    def check_Xfrac(self, silent_pass=False):
        if self.Xfrac_tot != 1:
            print('Total Xfrac still does not sum to 1 '
                  f'(Xfrac = {self.Xfrac_tot}).')
            return 0
        else:
            Xfrac_cum = 0.0
            for species_obj in self.species_list:
                Xfrac_cum += species_obj.Xfrac
            if Xfrac_cum != self.Xfrac_tot:
                print('ERROR: Xfrac_tot and cumaltive species Xfrac mismatch!')
                return 0
            if not silent_pass:
                print('Total Xfrac sums to 1!')
            return 1
        
        
    def finalize_initalization(self, normalize=False):
        # Check that Xfrac is 1.0
        if normalize:
            self.normalize_Xfrac()
        else:
            if self.check_Xfrac(silent_pass=True):
                return False
        # Truncate spectrum to ionizing region
        if len(self.species_list) < 1:
            print('ERROR: Can only truncate bins if given list of species.')
            return False
        min_I_ion = np.inf
        for species_obj in self.species_list:
            if species_obj.I_ion < min_I_ion:
                min_I_ion = species_obj.I_ion
        truncated_edge = const.hc/min_I_ion
        self.bin_breaks.sort()
        if self.bin_breaks[-1] != truncated_edge:
            if truncated_edge not in self.bin_breaks:
                new_breaks = list(self.bin_breaks)+[truncated_edge]
            else:
                new_breaks = list(self.bin_breaks)
            new_breaks.sort()
            edge_index = new_breaks.index(truncated_edge)
            new_breaks = new_breaks[:edge_index+1]
            self.bin_spectrum(new_breaks)
            
            
    def subbin(self, tau_cut=10, R2=1e-2):
        if self.n_bins == 0: # monochromatic
            return
        self.tau_cut = tau_cut
        self.R2 = R2
        self.original_bins = self.bin_breaks

        half_span = np.inf
        dwl = 0.0
        new_bins = [None]
        while new_bins != []:
            new_bins = []
            for ibin in range(self.n_bins):
                # Go from right
                r0_edge = 0.9999999*self.bin_breaks[ibin+1]
                l0_edge = 1.0000001*self.bin_breaks[ibin]
                max_l0 = l0_edge
                max_r0 = r0_edge
                last_sucessful0 = r0_edge
                last_fail0 = l0_edge
                # Loop checking bins: half_span > dwl
                while True:
                    half_span = (r0_edge-l0_edge)/2
                    wl = np.linspace(l0_edge, r0_edge, 4096)[1:-1]
                    hnu = const.hc/const.eV/wl
                    # Get relavant species in width
                    d3s_max = []
                    s_min = []
                    for s in self.species_list:
                        sigma = s.atomic_data.cross_section(const.hc/const.eV
                                                            /wl)
                        if sigma.any():
                            if not sigma.all():
                                print(f'WARNING: {s.atomic_data.name} '
                                      'has zero in array')
                            d3 = (s.atomic_data
                                  .cross_section_third_derivative(hnu,
                                                                  wrt='lambda'))
                            d3s_max.append(abs(d3).max())
                            s_min.append(sigma.min())
                    dwl0 = 0.0
                    for i in range(len(s_min)):
                        dwl0 += d3s_max[i]/s_min[i]
                    dwl0 = 1./dwl0
                    dwl0 *= 6.*R2/tau_cut
                    if dwl0 < 0:
                        print(s_min, d3s_max, dwl0)
                        sys.exit()
                    else:
                        dwl0 = (dwl0)**(1./3.)
                    if abs(1.-dwl0/half_span) < 1e-3:
                        break
                    if dwl0 > half_span:
                        if last_fail0 == l0_edge:
                            break
                        last_sucessful0 = l0_edge
                        l0_edge -= (l0_edge-last_fail0)/2
                    else:
                        last_fail0 = l0_edge
                        l0_edge += (last_sucessful0-l0_edge)/2
                # Go from left
                r1_edge = 0.999*self.bin_breaks[ibin+1]
                l1_edge = 1.00001*self.bin_breaks[ibin]
                max_l1 = l1_edge
                max_r1 = r1_edge
                last_sucessful1 = l1_edge
                last_fail1 = r1_edge
                # Loop checking bins: half_span > dwl
                while True:
                    half_span = (r1_edge-l1_edge)/2
                    wl = np.linspace(l1_edge, r1_edge, 4096)[1:-1]
                    hnu = const.hc/const.eV/wl
                    # Get relavant species in width
                    d3s_max = []
                    s_min = []
                    for s in self.species_list:
                        sigma = s.atomic_data.cross_section(const.hc/const.eV/wl)
                        if sigma.any():
                            if not sigma.all():
                                print(f'WARNING: {s.atomic_data.name} has zero in array')
                            d3 = (s.atomic_data
                                  .cross_section_third_derivative(hnu, wrt='lambda'))
                            d3s_max.append(abs(d3).max())
                            s_min.append(sigma.min())
                    dwl1 = 0.0
                    for i in range(len(s_min)):
                        dwl1 += d3s_max[i]/s_min[i]
                    dwl1 = 1./dwl1
                    dwl1 *= 6.*R2/tau_cut
                    if dwl1 < 0:
                        print(s_min, d3s_max, dwl1)
                        sys.exit()
                    else:
                        dwl1 = (dwl1)**(1./3.)
                    if abs(1.-dwl1/half_span) < 1e-3:
                        break
                    if dwl1 > half_span:
                        if last_fail1 == r1_edge:
                            break
                        last_sucessful1 = r1_edge
                        r1_edge -= (r1_edge-last_fail1)/2
                    else:
                        last_fail1 = r1_edge
                        r1_edge += (last_sucessful1-r1_edge)/2
                # hifsd
                add0 = True
                add1 = True
                if last_sucessful0 == r0_edge and last_fail0 == l0_edge:
                    add0 = False
                elif (r0_edge+l0_edge)/2-dwl0 <= max_l0:
                    add0 = False
                if last_sucessful1 == r1_edge and last_fail1 == l1_edge:
                    add1 = False
                elif (r1_edge+l1_edge)/2+dwl1 >= max_r1:
                    add1 = False
                loc0 = (r0_edge+l0_edge)/2-dwl0
                loc1 = (r1_edge+l1_edge)/2+dwl1
                width = r0_edge-l1_edge
                gap = width-2*(dwl0+dwl1)
                split = False
                if gap < dwl0 and gap < dwl1:
                    split = True
                    inc = width/3
                if add0 and add1:
                    if gap <= 0:
                        avg = ((r0_edge+l0_edge)/2-dwl0
                               +(r1_edge+l1_edge)/2+dwl1)/2
                        new_bins += [avg]
                    elif split:
                        new_bins += [max_l1+inc, max_l1+2*inc]
                    else:
                        new_bins += [(r1_edge+l1_edge)/2+dwl1,
                                     (r0_edge+l0_edge)/2-dwl0]
                elif add0:
                    new_bins += [(r0_edge+l0_edge)/2-dwl0]
                elif add1:
                    new_bins += [(r1_edge+l1_edge)/2+dwl1]
                else:
                    pass
            self.bin_spectrum(list(self.bin_breaks)+new_bins)
    
    
    def get_max_degree(self, Rn=1e-2):
        self.max_degree = [0,]*self.n_bins
        for ibin in range(self.n_bins):
            # Masked data
            mk = ((self.data['wl']>=self.bin_breaks[ibin])
                    &(self.data['wl']<=self.bin_breaks[ibin+1]))
            wl = self.data['wl'][mk]
            species = []
            for s in self.species_list:
                sigma = s.atomic_data.cross_section(const.hc/const.eV/wl)
                if sigma.any():
                    if not sigma.all():
                        print(f'WARNING: {s.atomic_data.name} '
                              'has zero in array')
                    species.append(s)
            term_A = 0.0
            term_B = 0.0
            delta_wl = (self.bin_breaks[ibin+1]-self.bin_breaks[ibin])
            wl0 = (self.bin_breaks[ibin+1]+self.bin_breaks[ibin])/2
            E0 = const.hc/wl0/const.eV
            for s in species:
                sigma = s.atomic_data.cross_section(const.hc/const.eV/wl)
                min_sigma = sigma.min()
                d_sigma = (s.atomic_data
                                  .cross_section_derivative(E0, wrt='lambda'))
                d2_sigma = (s.atomic_data
                                  .cross_section_second_derivative(E0,
                                                                   wrt='lambda'))
                term_A += d_sigma/min_sigma
                term_B += (delta_wl**2.)*d2_sigma/min_sigma/2.
            max_term = self.tau_cut*(term_A*delta_wl
                                     +term_B*delta_wl**2)
            i = 1
            while True:
                if i > 100:
                    # Abuse int being able to be raised to large powers, drop Rn
                    if int(max_term)**i <= math.factorial(i):
                        self.max_degree[ibin] = i#max(self.max_degree[ibin], i)
                        break
                elif i > 300:
                    # If still a problem just use Stirling's approximation
                    # max_degree[ibin] = int(max_term*np.exp(1))
                    self.max_degree[ibin] = int(max_term*np.exp(1))
                    break
                elif max_term**i <= math.factorial(i)*Rn:
                    self.max_degree[ibin] = i#max(self.max_degree[ibin], i)
                    break
                i += 1


    def set_abscissas(self):
        if self.max_degree is None:
            print('ERROR: First run get_max_degree()')
            return
        self.bin_absc = [None,]*self.n_bins
        self.bin_wgts = [None,]*self.n_bins
        for i in range(self.n_bins):
            n_points = int((self.max_degree[i]+1)/2+0.5)
            abscissas, weights = np.polynomial.legendre.leggauss(n_points)
            mask = ((self.data['wl']>=self.bin_breaks[i])
                    &(self.data['wl']<=self.bin_breaks[i+1]))
            masked_data = self.data[mask]
            left = masked_data.iloc[0]['wl']
            right = masked_data.iloc[-1]['wl']
            diff = (right-left)/2
            avg = (right+left)/2
            self.bin_absc[i] = diff*abscissas+avg
            self.bin_wgts[i] = diff*weights


    def write_csv(self, filename):
        if self.log_F_polys is None or self.log_Phi_polys is None:
            print('ERROR: First run smooth_spectrum()')
            return
        if self.max_degree is None:
            print('ERROR: First run get_max_degree()')
            return
        if self.bin_absc is None:
            print('ERROR: First run set_abscissas()')
            return
        headers = [r'$hc/\lambda_i$']
        headers += [r'$w_i\phi_{\lambda_i}$']
        for i, s in enumerate(self.species_list):
            headers += [rf'$\sigma_{{\lambda_i,{s.atomic_data.name}}}$']
        if self.n_bins == 0: # monochromatic
            FtoPHI = self.bin_breaks[0]/const.hc
            binrows = np.array([[1./FtoPHI, 1.]])
            for i, s in enumerate(self.species_list):
                cols = s.atomic_data.cross_section(1/(const.eV*FtoPHI))
                binrows = np.column_stack([binrows, cols])
            df = pd.DataFrame(binrows, columns=headers)
            df = df.convert_dtypes()
            comment = f'# NPTS: 1\n'
        else:
            table = []
            Phi_tot = 0.
            mk = ((self.data['wl']>=self.bin_breaks[0])
                    &(self.data['wl']<=self.bin_breaks[-1]))
            wl = self.data['wl'][mk]
            F = self.data['F_wl'][mk]
            Phi = F*wl/const.hc
            Phi_tot = integrate.trapz(Phi, wl)
            FtoPHI = 1./integrate.trapz(Phi/Phi_tot*const.hc/wl, wl)
            npts = 0
            for ibin in range(self.n_bins):
                col1 = const.hc/self.bin_absc[ibin]
                hnu = const.hc/const.eV/(self.bin_absc[ibin])
                col2 = (self.bin_wgts[ibin]*self.rescale[ibin]/Phi_tot
                        *np.exp(self.log_Phi_polys[ibin](self.bin_absc[ibin]))
                       )
                binrows = np.column_stack([col1, col2])
                for i, s in enumerate(self.species_list):
                    cols = s.atomic_data.cross_section(hnu)
                    binrows = np.column_stack([binrows, cols])
                table.append(binrows)
                npts += len(self.bin_absc[ibin])
            table = np.vstack(table)
            df = pd.DataFrame(table, columns=headers)
            df = df.convert_dtypes()
            comment = f'# NPTS: {npts}\n'
        comment += f'# NSPECIES: {self.n_species}\n# FtoPHI: {FtoPHI:.12e}\n'
        comment += f'# DATE: {self.date}\n'
        comment += f'# WINDOW: {self.bin_breaks[0]},{self.bin_breaks[-1]}\n'
        comment += '# IONPOTS: '
        for s in self.species_list:
            comment += f'{s.I_ion:.12e}, '
        comment = comment[:-2]+'\n# '  
        for h in headers:
            comment += h+', '
        comment = comment[:-2]+'\n'
        with open(filename, 'w') as file_:
            file_.write(comment)
            df.to_csv(file_, index=False, header=False, float_format='%.17e')
            
        return df