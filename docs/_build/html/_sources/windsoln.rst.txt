.. _windsoln:

Analysis Module
=====================

.. automodule:: wrapper.wrapper_utils.windsoln

Perform postfacto calculations
------------------------------------
Any values that can be computed from velocity, temperature, density, and ionization structure are (and can be, upon requested) included in ``add_user_vars()`` and are accessible via the dataframe ``sim.windsoln.soln[]``

- ``sim.windsoln.soln[]`` values are all as a function of radius (normalized values available via ``sim.windsoln.soln[]``):
   - Radius - 'r' (cm), 'z' (:math:`r = R_min + q*z`), 
   - Mass density - 'rho' (g/cm^3)
   - Velocity -'v' (cm/s)
   - Temperature - 'T', 
   - Per-element "Neutral" fraction - 'Ys_HI', etc. ("neutral" refers to the the lowest ionization state in the simulation as **Wind-AE** can currently only model one ionization state per metal)
   - Per-element Column density  -'Ncol_HI', etc.
   - Number densities (1/cm^3) - 
      - Total electron - 'n_e'
      - Element - Total: 'n_H', neutral: 'n_HI', ionized: 'n_HII'. Such that 'n_H' = 'n_HI' + 'n_HII'
      - Total of all species: 'n_tot'
   - Mean molecular/atomic weight - 'mu' (dimensionless)
   - Pressure - 'P' (dyne/cm^2 = barye = microbar = :math:`10^{-6}` bar)
   - Ram pressure - 'ram' (barye)
   - Sound speed - 'cs' (cm/s)
   - Mach number - 'Mach' (dimensionless)
   - Density scale height - 'Hsc' (cm)

   - Multispecies Ionization Balance:
      - Ionization rate (per unit volume)- 'ion_rate_HI', etc. (:math:`s^{-1} cm^{-3}`)
      - Recombination rate (per unit volume) - 'recomb_HI', etc. (:math:`s^{-1} cm^{-3}`)
      - Advection rate (per unit volume) - 'advec_HI', etc. (:math:`s^{-1} cm^{-3}`)

   - Heating and Cooling Rates (per unit volume):
      - Currently Included in **Wind-AE**:
         - Photoionization rate - 'heat_ion' (ergs :math:`s^{-1} cm^{-3}`)
         - Recombination cooling - 'cool_rec' (ergs :math:`s^{-1} cm^{-3}`)
         - Advective Heating/Cooling - 'heat_advect' (ergs :math:`s^{-1} cm^{-3}`)
         - Bolometric Heating/Cooling - 'boloheat', 'bolocool' (ergs :math:`s^{-1} cm^{-3}`)
         - Line cooling:
            - Lyman-alpha cooling - 'cool_lyman' (ergs :math:`s^{-1} cm^{-3}`)
            - Carbon line cooling - "cool_CII_1570000A" (15700000 :math:`\AA` line), "cool_CII_2326A","cool_CII_1334A", "cool_CIII_1910A", "cool_CIII_977A"
            - Oxygen line cooling - "cool_OII_834A" (834 :math:`\AA` line), "cool_OII_2741A", "cool_OII_3727A", "cool_OII_7320A", "cool_OIII_520000A", "cool_OIII_5000A", "cool_OIII_166A", "cool_OIII_84A"
         - PdV cooling (work done by expanding gas) - 'cool_PdV' (ergs :math:`s^{-1} cm^{-3}`)
      - Not included in **Wind-AE**:
         - Conductive cooling - 'cool_cond' (ergs :math:`s^{-1} cm^{-3}`) **[PLANNED]**
         - Gravitational heating/cooling - 'cool_grav' (ergs :math:`s^{-1} cm^{-3}`)
      - Cumulative differential heating - 'cum_heat' (ergs :math:`s^{-1} cm^{-3}`)

   - Bernoulli Constants - 'bern', 
   - Hydrogen Knudsen number - 'Kn_hb_HI' (hardbody), 'Kn_Co_HI' (Coloumb), 'Kn_mx_HI' (mix of hardbody and Coloumb) [multispecies PLANNED] 

- A selection of helpful ``sim.windsoln`` variables:
   - ``soln`` - DataFrame containing all solution variables that vary as a function of radius.
   - ``soln_norm`` - DataFrame containing all solution variables normalized by their scales (then ``sim.windsoln.soln['rho'] = sim.windsoln.soln_norm['rho'] * sim.windsoln.scales_dict['rho']``)
   - ``scales_dict`` - Dictionary containing the scale factors for temperature, density, radius, etc.
   - Planet parameters: ``Mp`` (g), ``Rp`` (cm), ``semimajor`` (cm), ``Mstar`` (g), bolometric luminosity ``Lstar`` (ergs/s), flux at semimajor axis ``Ftot`` (ergs :math:`s^{-1} cm^{-2}`)
   - Physics parameters:  ``nspecies`` (number of unique elements), ``HX`` (mass fraction array of elements), ``species_list``, ``atomic_masses`` (g), ``molec_adjust`` (:math:`m_H`, average mean molecular weight of region of the atmosphere below the wind)
   - Spectrum parameters:
      - ``npts`` - number of wavelength points in smoothed spectrum
      - ``E_wl`` - energy per wavelength bin (ergs)
      - ``wPhi_wl`` - photon density normalized by total flux (:math:`=\phi e^{-\tau} / F_{tot}`) ()
      - ``F_wl`` - flux per wavelength bin (ergs :math:`cm^{-2} s^{-1}`) (``sum(F_wl) = Ftot``)
      - Generating spectrum: ``spec_src_file``, ``spec_kind``, ``spec_window`` (nm), etc. 
   - Boundary Condition parameters: 
      - Lower: ``Rmin`` (Rp), ``rho_rmin`` (RHO0 = sim.windsoln.scales_dict['rho']), ``T_rmin`` (T0), per-element neutral fraction ``Ys_rmin_HI``, etc. 
      - Upper: ``Rmax`` (Rp), per-element ``Ncol_sp_HI``, etc. (g/cm^2)
   - Flags: 
      - ``lyacool`` - Line cooling on/off.
      - ``tidalforce`` - Tidal force from star on/off.
      - ``bolo_heat_cool`` - Complementary error function that governs the transition from layer below wind where molecules may be present and bolometric heating / cooling dominate the energy budget on/off.
      - ``integrate_outward`` - Integrate outward from the lower boundary to the Coriolis radius (``Rmax = R_cori``) on/off.
   - Useful Radii:
      - Sonic point radius - ``R_sp`` (Rp)
      - XUV radius / wind launch radius / photoionization front - ``R_XUV`` (Rp), lowest radial extent of substantial photoionization energy deposition (i.e., where wind "launches")
      - Coriolis radius - ``R_cori`` (Rp)
      - Hill sphere - ``R_hill`` (Rp)
      - Exosphere - ``R_exo`` (Rp)
   - Optical depth - ``taus`` (dimensionless) 2D array of size [len(radii),len(spectrum bins)] (Not to be confused with the function ``tau_array``)


Assorted Useful Functions
----------------------------------
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_mu
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_massloss
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.tau_array
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_tau1_radius

Rarely useful and somewhat deprecated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_R_exo
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_Jeans
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_Coriolis
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_ballistic
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.regrid