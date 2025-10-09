---
title: 'Wind-AE: A Python package for atmospheric escape with metals and X-rays'
tags:
  - Python
  - astronomy
  - exoplanets
  - atmospheric escape
  - Parker wind
  - relaxation method
authors:
  - name: Madelyn I. Broome
    orcid: 0000-0002-7520-5663
    equal-contrib: true
    affiliation: 1 
  - name: John H. McCann
    equal-contrib: true
    affiliation: 1
  - name: Ruth Murray-Clay
    orchid: 0000-0001-5061-0462
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: University of California Santa Cruz, Dept. of Astronomy & Astrophysics, Santa Cruz, US
   index: 1
date: 25 September 2025
bibliography: paper.bib

aas-journal: Astrophysical Journal
---

# Summary
Throughout their lives, close-in exoplanets are strongly irradiated by their host stars. In the first few hundred million years of life, though, atomic species in a planet’s upper atmosphere are photoionized by the young hot stars' strong X-ray and extreme UV radiation. This ionizing radiation from the star heats the planets’ upper atmospheres to up to 10,000 K and causes the gas to expand, accelerate from sub- to supersonic speeds (a.k.a., a Parker wind), and outflow from the planet – essentially photoevaporating part or all of the planet’s atmosphere in a process called atmospheric escape. Understanding the atmospheric mass loss histories of these planets is essential to understanding the evolution of exoplanets close to their host stars. However, mass loss rates are not directly observable. They can only be inferred from models.
To that end, we have developed `Wind-AE`: a fast, 1D, steady-state Parker wind forward relaxation model based on [@rmc:2009] that can model irradiation by multifrequency X-ray/extreme-UV (XUV) stellar spectra and the interactions of these high energy photons with "metals" (atomic species).


# Statement of Need

`Wind-AE` is a Python wrapper around an update of the widely-cited [@rmc:2009] C code. The update to the model physics primarily includes the ability to model multifrequency XUV photons interacting with metals in the planets' upper atmospehres. Modeling impact of the presence of metals and X-rays is important because the two have the ability to dramatically change the mass loss rates and wind radial structure for certain planets. `Wind-AE` is unique for its combination of X-ray/metal physics and speed. Speed is of particular interest to users in the broader exoplanet communities for whom `Wind-AE` presents a variety of uses including as a simple plug-and-play model that quickly returns mass loss rates and as a model that can be easily integrated with more sophisticated models to: simulate observability of atmospheric escape (as in [@anusha:2025]), simulate the evolution of planets' atmospheres over time (as in Tang et al. 2025 (submitted)), to investigate open questions about the demographics of exoplanets (as in [@parke:2025]), and make grids of mass loss rates (as in Broome et al. (2025)).

Because of the inherent sensitivity of the relaxation method (it requires a good initial guess that is close in parameter space to the goal solution to converge), a Python wrapper that contains algorithms for ramping smoothly through parameter space makes `Wind-AE` far more reliable, flexible, and user friendly. Python-based user interfaces are far more accessible than C to the majority astronomers, from researchers to undergraduate students who have used this code.
<!-- `Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. -->

<!-- # Verification and Documentation -->

<!-- `orbitize!` implements a full stack of automated testing and documentation building 
practices. We use GitHub Actions to automatically run a suite of unit tests, maintained in [orbitize/tests](https://github.com/sblunt/orbitize/tree/main/tests),
each time code is committed to the public repository or a pull request is opened. The Jupyter notebook
tutorials, maintained in [orbitize/docs/tutorials](https://github.com/sblunt/orbitize/tree/main/docs/tutorials), are also run automatically when a 
pull request to the `main` branch is opened. Documentation is built using `sphinx`, and hosted
on readthedocs.org at [orbitize.info](https://orbitize.readthedocs.io/en/latest/). We also
maintain a set of longer-running tests in [orbitize/tests/end-to-end-tests](https://github.com/sblunt/orbitize/tree/main/tests/end-to-end-tests) that show real
scientific use cases of the code. These tests are not automatically run.

`orbitize!` is documented through API docstrings describing individual functions, which are accessible on [our readthedocs page](https://orbitize.readthedocs.io/en/latest/api.html), a set of [Jupyter notebook tutorials](https://orbitize.readthedocs.io/en/latest/tutorials.html) walking the user through a particular application, a set of [frequently asked questions](https://orbitize.readthedocs.io/en/latest/faq.html),
and an in-progress ["manual"](https://orbitize.readthedocs.io/en/orbitize-manual/manual.html) describing orbit fitting with `orbitize!` from first principles. -->

# Comparison to Similar Packages
For detailed scientific comparisons of other 1D open source atmospheric escape models, see Appendix A of [our paper](). Broadly speaking, `Wind-AE` runs faster than many models that perform more detailed physical calculations and makes what we consider to be more physically realistic assumptions than many models of comparable or faster speed.

- [ATES](https://github.com/AndreaCaldiroli/ATES-Code) is a comparably rapid model that includes X-ray radiation from the star, but differs from `Wind-AE` in that it approximates the stellar spectrum as a power law, only models atomic hydrogen and helium in the atmosphere, and calculates the ionization balance post facto [@ates:2021].
- [p-winds](https://github.com/ladsantos/p-winds) is an isothermal Parker wind model which is approximately 10 times faster than `Wind-AE`, but serves a slightly different purpose as users set the mass loss rate to obtain the wind's structure. `p-winds` models extreme-UV radiation (EUV), but does not model X-rays. By solving for the temperature as a function of radius, `Wind-AE` is able to obtain more realistic outflow profiles, but at the cost of speed. ([@pwinds:2022]).
- [Sunbather](https://github.com/antonpannekoek/sunbather) wraps `p-winds` and Cloudy (a spectral synthesis code [@cloudy:2017]) to create transit models that include metals [@sunbather:2024].
- [TPCI](https://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/576/A21) leverages Cloudy and the hydrodynamic code PLUTO ([@pluto:2012] ) for more thorough X-ray and EUV-irradiated simulations, but slightly more expensive calculations ([@salz:2015]).
- [pyTPCI](https://ascl.net/2506.012) is a Python wrapper for TPCI [@pytpci:2025].
- [AIOLIS](https://github.com/Schulik/aiolos) models diffusion (a physical process whose impact on atmospheric escape is still being investigated) and the lower atmosphere in more detail and in a time-dependent manner, so it takes of order 100 times longer [schulik:2022].
- CETIMB includes more detailed lower atmosphere physics and is therefore about 500 times slower [@koskinen:2022].
- CHAIN integrates CLOUDY and [@kubyshkina:2018]. Like TPCI, via Cloudy, CHAIN has slightly more nuanced heating and cooling calculations and does model X-rays (though the published grid is EUV only [@grid:2021]), so is also more expensive.


# Acknowledgements
MIB, RMC, and JHM acknowledge funding support from NSF CAREER grant 1663706. MIB and RMC acknowledge support from NASA'S Interdisciplinary Consortia for Astrobiology Research (NNH19ZDA001N-ICAR) under grant number 80NSSC21K0597.

We would like to thank Jorge Fernandez-Fernandez, Emerson Tao, and Artem Aguichine for beta testing `Wind-AE`.

# References