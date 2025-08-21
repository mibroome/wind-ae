.. _usage:

Is **Wind-AE** the right tool for me?
=========================================

**Wind-AE** is well-suited for users interested in **quickly estimating mass loss rates** or outflow structure. 
Outflow structure includes bulk temperature and per-species ionization fractions as a function of radius, so can be 
easily translated into approximating and **predicting observables and transits**, including metastable helium 
(He 10830:math:`\AA`) transits, though a He transit module is not yet included. Precise modeling of lower atmosphere 
(:math:`\lesssim 100` microbar) is considered necessary for highly accurate transit models, but Wind-AE can be easily 
coupled to lower atmosphere photochemistry models whose outputs (e.g., radius, temperature, abundances, ionization fractions,
etc. at 1 microbar) can be fed into Wind-AE as inputs.

.. note::
   *If you are interested in outflow structure:* Past the Coriolis turning radius (a few planetary radii) 3D physics 
   dominates, so **Wind-AE** does not integrate past that point. **Wind-AE** also makes simplifying assumptions about 
   the region below the wind-launch radius (~10 nanobars).

Because **Wind-AE** runs on the order of seconds to minutes, it can be (and has been) used to **model planet evolution**.

**Wind-AE** can model:
----------------------

- Multiple atomic species
- X-ray physics (secondary ionizations and K-shell ionization cross-sections for relevant metals)
- Both low and high stellar XUV flux
- **Heating & Cooling**: Ionization heating, bolometric heating & cooling (negligible in wind), PdV cooling (work done due to expansion of gas), radiative / atomic line cooling (Lyman-\alpha, OI, OII, OIII, CI, CII), recombination cooling

**Wind-AE** does not (currently) include:
-----------------------------------------

- **Magnetic fields**
- **Time dependence**
- **Diffusion/drag** - the atomic species set by the user are assumed to be entrained in the outflow and in thermal equilibrium. This is an appropriate assumption for species below the `crossover mass <https://ui.adsabs.harvard.edu/abs/1987Icar...69..532H>`_ and a warning will be raised.
- **Heating & Cooling**: 
   - Conduction (warning raised if relevant, **planned**)
   - H3+ line cooling (not planned)
   - Fe & Ca line cooling (relevant at high Z only, **planned**)
   - Free-free cooling (warning raised if relevant, not planned)
- Multiple ionization states of the same species (**planned**)

See Broome et al. (submitted) for more details.

Other tools and models:
-----------------------

- Want a rapid H/He model with power-law approximated XUV spectra? Check out `ATES <https://github.com/AndreaCaldiroli/ATES-Code>`_ (`Caldiroli et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..30C/abstract>`_)
- Do you want to set the mass loss rate (:math:`\dot{M}`) yourself or want an EUV isothermal Parker wind outflow model that runs in nanoseconds? Check out `p-winds <https://github.com/ladsantos/p-winds>`_ (`Dos Santos et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26A...659A..62D/abstract>`_)
- Do you want to use p-winds and get transit models for metals via Cloudy? Check out `Sunbather <https://github.com/antonpannekoek/sunbather>`_ (`Linssen et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..43L/abstract>`_)
- Want to leverage Cloudy and the hydrodynamic code PLUTO for more thorough XUV-irradiated, but slightly more expensive calculations? Check out `TPCI <https://ui.adsabs.harvard.edu/abs/2015A%26A...576A..21S/abstract>`_ (`Salz et al. 2015`)
- That sounds great, but you prefer to code in Python over C/C++? Check out `pyTPCI <https://ui.adsabs.harvard.edu/abs/2025ApJ...980...34R/abstract>`_ [Riley, Zhang, & Bean 2025]
- Do you care about diffusion throughout the wind? Check out `AIOLIS <https://github.com/Schulik/aiolos>`_ (`Schulik & Booth, 2022 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.523..286S/abstract>`_)
- Want to model the lower atmosphere in more detail? Check out CETIMB (Koskinen et al. 2022)
- Just want a grid of mass loss rates for pure-Hydrogen, low-flux-EUV-irradiated planets? See `Kubyshkina & Fossati <https://ui.adsabs.harvard.edu/abs/2021RNAAS...5...74K/abstract>`_
- Want a grid of mass loss rates for pure-Hydrogen, high-flux-XUV-irradiated planets? See `Owen & Jackson (2012) <https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.2931O/abstract>`_

.. note::
   Want your model added to this list or to update the short bio? Email mabroome@ucsc.edu