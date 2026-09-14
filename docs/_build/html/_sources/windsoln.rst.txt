.. _windsoln:

Analysis Module
=====================

.. automodule:: wind_ae.wrapper.wrapper_utils.windsoln

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
      - Included in the energy equation when their corresponding flag is on (see :ref:`relax_wind`'s :ref:`heating-and-cooling-flags` and the ``Flags`` variables below):
         - Photoionization heating - 'heat_ion' (ergs :math:`s^{-1} cm^{-3}`) (always on)
         - Line cooling (governed by ``linecool`` flag, formerly ``lyacool``):
            - Lyman-alpha cooling - 'cool_lyman' (ergs :math:`s^{-1} cm^{-3}`)
            - Carbon line cooling - "cool_CII_1570000A" (15700000 :math:`\AA` line), "cool_CII_2326A","cool_CII_1334A", "cool_CIII_1910A", "cool_CIII_977A"
            - Oxygen line cooling - "cool_OII_834A" (834 :math:`\AA` line), "cool_OII_2741A", "cool_OII_3727A", "cool_OII_7320A", "cool_OIII_520000A", "cool_OIII_5000A", "cool_OIII_166A", "cool_OIII_84A"
         - Recombination cooling - 'cool_rec' (ergs :math:`s^{-1} cm^{-3}`) (governed by ``recombo_cool`` flag; on by default)
         - Free-free (bremsstrahlung) cooling - 'cool_free' (ergs :math:`s^{-1} cm^{-3}`) (governed by ``free_free_cool`` flag; on by default) **[NEW in v2.0]**
         - Conductive cooling - 'cool_cond' (ergs :math:`s^{-1} cm^{-3}`) (governed by ``conduction`` flag; **off** by default --- see warning below) **[NEW in v2.0]**. The associated total (electron + neutral) thermal conductivity coefficient is stored in 'cond_coeff' (erg :math:`s^{-1} cm^{-1} K^{-1}`).
         - Advective Heating/Cooling - 'heat_advect' (ergs :math:`s^{-1} cm^{-3}`)
         - Bolometric Heating/Cooling - 'boloheat'/'heat_bolo' (heating), 'bolocool'/'cool_bolo' (cooling) (ergs :math:`s^{-1} cm^{-3}`) (governed by the ``bolo_heat_cool`` scaling flag; the 'heat_bolo'/'cool_bolo' names are aliases added in v2.0 for readability, identical values to 'boloheat'/'bolocool')
         - PdV cooling (work done by expanding gas) - 'cool_PdV' (ergs :math:`s^{-1} cm^{-3}`)
      - Diagnostic only (always computed, but never included self-consistently in the energy equation):
         - Gravitational heating/cooling - 'cool_grav' (ergs :math:`s^{-1} cm^{-3}`)
      - Cumulative differential heating - 'cum_heat' (ergs :math:`s^{-1} cm^{-3}`)

   - Bernoulli Constants - 'bern', 
   - Hydrogen Knudsen number - 'Kn_hb_HI' (hardbody), 'Kn_Co_HI' (Coloumb), 'Kn_mx_HI' (mix of hardbody and Coloumb) [multispecies PLANNED] 

.. warning::
   As of v2.0, conductive, recombination, and free-free cooling can each be toggled
   independently via flags (see below). Turning on conduction (``sim.turn_on_conduction()``)
   makes the ODEs very stiff --- only do so after all other ramping is complete. The
   temperature dependence of free-free and recombination cooling can make Wind-AE
   numerically unstable, and ramping to higher metallicities with them on may fail.

- A selection of helpful ``sim.windsoln`` variables:
   - ``soln`` - DataFrame containing all solution variables that vary as a function of radius.
   - ``soln_norm`` - DataFrame containing all solution variables normalized by their scales (then ``sim.windsoln.soln['rho'] = sim.windsoln.soln_norm['rho'] * sim.windsoln.scales_dict['rho']``)
   - ``scales_dict`` - Dictionary containing the scale factors for temperature, density, radius, etc.
   - Planet parameters: ``Mp`` (g), ``Rp`` (cm), ``semimajor`` (cm), ``Mstar`` (g), bolometric luminosity ``Lstar`` (ergs/s), flux at semimajor axis ``Ftot`` (ergs :math:`s^{-1} cm^{-2}`)
   - Physics parameters:  ``nspecies`` (number of unique elements), ``HX`` (mass fraction array of elements), ``species_list``, ``atomic_masses`` (g), ``molec_adjust`` (:math:`m_H`, average mean molecular weight of region of the atmosphere below the wind), ``kappa_opt`` (optical opacity used in the molecular/bolometric layer, default 4e-3), ``kappa_IR`` (IR opacity used in the molecular/bolometric layer, default 1e-2), ``gamma`` (adiabatic index, default 5/3) **[kappa_opt, kappa_IR, and gamma are now per-planet, user-settable/rampable attributes as of v2.0; previously fixed]**, ``chain_parent`` (list, length ``nspecies``; ``chain_parent[j]`` gives the index of the parent species if species ``j`` is one ionization state above another species in ``species_list``, or -1 if independent --- see :ref:`metals`) **[NEW in v2.0]**
   - Spectrum parameters:
      - ``npts`` - number of wavelength points in smoothed spectrum
      - ``E_wl`` - energy per wavelength bin (ergs)
      - ``wPhi_wl`` - photon density normalized by total flux (:math:`=\phi e^{-\tau} / F_{tot}`) ()
      - ``F_wl`` - flux per wavelength bin (ergs :math:`cm^{-2} s^{-1}`) (``sum(F_wl) = Ftot``)
      - Generating spectrum: ``spec_src_file``, ``spec_kind``, ``spec_window`` (nm), etc. 
   - Boundary Condition parameters: 
      - Lower: ``Rmin`` (Rp), ``rho_rmin`` (RHO0 = sim.windsoln.scales_dict['rho']), ``T_rmin`` (T0), per-element neutral fraction ``Ys_rmin_HI``, etc. 
      - Upper: ``Rmax`` (Rp), per-element ``Ncol_sp_HI``, etc. (g/cm^2)
   - Flags (``flags_tuple``, in order): 
      - ``integrate_outward`` - Integrate outward from the lower boundary to the Coriolis radius (``Rmax = R_cori``) on/off.
      - ``tidalforce`` - Tidal force from star on/off.
      - ``linecool`` - Line cooling (Lyman-alpha and, if present, metal line cooling) on/off. **Renamed from** ``lyacool`` **in v2.0.**
      - ``bolo_heat_cool`` - Turns on (1.0) and off (0.0) and ramps (value between 0.0 and 1.0) in the complementary error function that governs the transition from layer below wind where molecules may be present and bolometric heating / cooling dominate the energy budget on/off.
      - ``conduction`` - Self-consistent conductive heat flux on/off. **[NEW in v2.0]** Off by default; introduces numerically stiff ODEs, so should only be turned on after all other ramping is complete (``sim.turn_on_conduction()`` / ``sim.turn_off_conduction()``).
      - ``recombo_cool`` - Recombination cooling on/off. **[NEW in v2.0]** Off by default.
      - ``free_free_cool`` - Free-free (bremsstrahlung) cooling on/off. **[NEW in v2.0]** Off by default. May cause numerical instability at high metallicity.
      - ``molec_layer`` - Turns on (1.0) and off (0.0) and ramps (value between 0.0 and 1.0) in the complementary error function that governs the mean-molecular-weight (mu) molecular-to-atomic transition, independent of ``bolo_heat_cool``. **[NEW in v2.0]** (``sim.turn_on_molecular_layer()`` / ``sim.turn_off_molecular_layer()``). Defaults to the value of ``bolo_heat_cool`` when loading older solution files that predate this flag.
   - Useful Radii:
      - Sonic point radius - ``R_sp`` (Rp), found inherently as a part of the relaxation method (i.e., C source code).
      - Wind launch radius - ``R_launch`` (Rp), computed as radius where photoionization heating dominates over other heating/cooling processes, which corresponds to where wind begins to accelerate.
      - Photoionization front / absorption radius - ``R_absorb`` (Rp), flux-per-wavlength-bin weighted mean :math:`\tau(\lambda)=1` surface.
      - XUV radius - ``R_XUV``=``R_launch`` (Rp), legacy.
      - Coriolis radius - ``R_cori`` (Rp), radius at which Coriolis force deflects wind 1 radian.
      - Hill sphere - ``R_hill`` (Rp)
      - Exosphere - ``R_exo`` (Rp), currently computed from H-H Knudsen number.
      - Roche lobe (L1 point) - ``R_roche`` (Rp), set after calling :meth:`~wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_roche_lobe` (computed automatically as part of ``add_user_vars()`` whenever ``integrate_outward=True``)
   - Optical depth - ``taus`` (dimensionless) 2D array of size [len(radii),len(spectrum bins)] (Not to be confused with the function ``tau_array``)


Assorted Useful Functions
----------------------------------
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.add_user_vars
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_mu
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_massloss
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.current_metallicity
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.tau_array
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_tau1_radius
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.alpha_rec

Rarely useful and somewhat deprecated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_R_exo
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_Jeans
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_Coriolis
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_ballistic
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.calc_roche_lobe
.. automethod:: wind_ae.wrapper.wrapper_utils.windsoln.wind_solution.regrid
