.. _spectrum:

Spectrum Module
=====================
.. Note::
   We avoid the cost of running a high resolution spectrum by fitting a polynomial to the input stellar spectrum.
   Any observed spectrum---such as the FISM2 solar spectrum from the LISIRD database that is the defualt solar spectrum---or
   realistic simulated XUV spectrum will vary widely in flux and shape across the spectral range, as well as be high resolution,
   making fitting polynomials difficult. If the the spectral qualities can be well approximated by low degree polynomial(s), though,
   it is inexpensive to accurate perform numerical integrations using Gauss-Legendre quadrature. Thus, we use the smoothing and binning
   algorithm (discussed in more detail in `McCann (2021) <https://escholarship.org/uc/item/0qc1m9zn>`_).

   Logarithmic fits and/or the least squares method would not locally (or, potentially, even globally) conserve energy along the 
   spectrum, thus we employ a Savitzky-Golay filter, which smooths evenly-spaced noisy data with a rolling polynomial.
   First we smooth the peaks the troughs of the spectrum by multi-passing the spectrum through the Savitzky-Golay filter.  
   The effect of running a Savitzky-Golay filter on small segment is similar to running a single pass filter on a larger wavelength 
   range, but it distorts the data less than a standard, larger single pass filter and better preserves the area under the smoothed 
   spectrum. The filtered spectrum is then renormalized to conserve total energy in each bin. Next a 5th degree polynomial is fit to 
   the filtered spectrum, again rescaling to preserve energy in each bin. The polynomial is calculated by used a spline with with an 
   infinite smoothing factor, which relaxes the spline to a single bin interval.

   Binning the spectrum allows us to run the multipass filtering in fewer smoothing passes and allows us to more accurately preserve 
   the spectrum shape, especially at the ionization wavelengths of species present in the model. Subbinning, in particular, allows us 
   to fit the spectrum with more low order polynomials, as opposed to fewer polynomials that would have to be higher order and would 
   be more difficulty to accurately and cheaply integrate using Gauss-Legendre quadrature. For our bins, we choose bin width 2:math:`r` 
   centered at wavelength :math:`\lambda_0` such that the error over :math:`\lambda\in[\lambda_0-r,\lambda_0+r]`  is less than 
   :math:`\epsilon`. 

   Bin edges are also placed at the ionization energies and K-shell ionization energies of the species present, unless one of the 
   ionization energies of an existing species is within 2 nm of an existing species' ionization energy.

   The physical effects of smoothing a spectrum are also mitigated by using the above method. 
   Ionizing energy is conserved since the peaks of the spectrum are smoothed and distributed locally---
   meaning there will be an equal amount of higher and lower than the peak ionization energy photons in the wind. 
   That being said, at the edges of the spectrum, where there are not necessarily symmetric peaks over which to smooth, this 
   method may over or underestimate the number of higher or lower energy photons. However, we take this smoothed approximation 
   for a high resolution spectrum to be adequate for **Wind-AE**.

Some of the functions in this module are redundant, but we retain them to allow for future updates.

********************************************
How to manually load/generate a spectrum
*********************************************
.. Note:: Changing the wavelength range can be more easily done via :meth:`wind_ae.wrapper.relax_wind.ramp_spectrum()`

This is not often necessary, but can be tricky, so an example is provided here.

Nominally `set_resolved()`, `set_normalized()`, and `set_window()` can all accept different wavelength ranges, but **Wind-A** has been simplified to accept the same wavelength range for all three. 

This ensures that `F_tot` is truely the total flux over the given range.

.. code-block:: python
   
   from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
   from wind_ae.wrapper.wrapper_utils.spectrum import spectrum  
   
   # For a scaled solar spectrum, set lisird=True
   spec = spectrum(lisird=False,spectrum_file='hd189733',wl_norm=1e-7)
   # Add species for spectrum binning/smoothing
   for name in ['HI','HeI','CI']:
      spec.add_species(name)
   # Setting the spectrum range in NANOMETERS
   # 13.6eV = 91nm, 2000eV = 0.6nm
   soln_resolved = [0.6, 91]
   spec.set_resolved(*soln_resolved)
   spec.set_normalized(*soln_resolved)
   kind = 'full' 
   spec.set_window(*soln_resolved, kind=kind)
   # Can plot via spec.plot() or spec.binning_plot() to confirm before generating
   spec.generate(kind=kind, savefile='wind_ae/inputs/spectrum.inp')

.. automodule:: wrapper.wrapper_utils.spectrum
   :members:
   :undoc-members:
   :show-inheritance: