# Wind-AE (BETA)

<h1 style="text-align: center;"><a href="https://wind-ae.readthedocs.io/en/latest/">Read the Docs</a></h1>

`Wind-AE` (pronounced /windy/) stands for "wind atmospheric escape" and is a relatively fast 1D, steady-state, hydrodynamic, non-isothermal, Parker wind relaxation code for modeling atmospheric escape based on [Murray-Clay et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...693...23M/abstract). `Wind-AE` is a forward model that solves the energy conservation, momentum conservation, and ionization equilibrium equations at the substellar point using the ["Numerical Recipes in C"](https://ui.adsabs.harvard.edu/abs/1992nrca.book.....P/abstract) relaxation method. This allows `Wind-AE` to quickly compute the atmosphere mass loss rate as well as upper atmosphere (<~100 microbar) velocity, temperature, density, and ionization structure as a function of altitude. 

`Wind-AE` updates [Murray-Clay et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...693...23M/abstract) to allow for the modeling of atomic metals and multifrequency XUV stellar spectra (Broome et al. submitted). If you use `Wind-AE`, please consider citing Broome et al. (submitted). 


We appreciate your patience while the docs are developed. In the meantime, take advantage of `Notebooks/Quickstart.ipynb` to get a quick orientation to `Wind-AE` and please report any bugs via [Github](https://github.com/mabroome/wind-ae/issues) or via email to mabroome@ucsc.edu.

![Build Status](https://github.com/sblunt/orbitize/actions/workflows/python-package.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/mibroome/wind-ae/badge.svg?branch=main)](https://coveralls.io/github/mibroome/wind-ae?branch=main)
![PyPI - Version](https://img.shields.io/pypi/v/wind_ae)
![Read the Docs (version)](https://img.shields.io/readthedocs/wind-ae/latest)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)](https://semaphorep.github.io/codeastro/)

Is Wind-AE the right tool for me?
----------------------
`Wind-AE` is well-suited for users interested in **quickly estimating mass loss rates** or outflow structure. Outflow structure includes bulk temperature and per-species ionization fractions as a function of radius, so can be easily translated into approximating and **predicting observables and transits**, including metastable helium (He 10830$\AA$) transits, though a He transit module is not yet included. Precise modeling of lower atmosphere ($\lssim 100 \mu$bar) is considered necessary for highly accurate transit models, but Wind-AE can be easily coupled to lower atmosphere photochemistry models whose outputs can (e.g., radius, temperature, abundances, ionization fractions, etc. at 1 $\mu$bar) can be fed into Wind-AE as inputs.  

>*If you are interested in outflow structure:* Past the Coriolis turning radius (a few planetary radii) 3D physics dominates, so `Wind-AE` does not integrate past that point. `Wind-AE` also makes simplifying assumptions about the region below the region below the wind-launch radius (~10 nanobars). 

Because `Wind-AE` runs on the order of seconds to minutes, it can be (and has been) used to **model planet evolution**.

#### `Wind-AE` can model:
- Multiple atomic species
- X-ray physics (secondary ionizations and K-shell ionization cross-sections for relevant metals)
- Both low and high stellar XUV flux
- **Heating & Cooling**: Ionization heating, bolometric heating & cooling (negligible in wind), PdV cooling (work done due to expansion of gas), radiative / atomic line cooling (Lyman-$\alpha$, OI, OII, OIII, CI, CII), recombination cooling
#### `Wind-AE` does not (currently) include:
- **Magnetic fields**
- **Time dependence**
- **Diffusion/drag** - the atomic species set by the user are assumed to be entrained in the outflow and in thermal equilibrium. This is an appropriate assumption for species below the [crossover mass](https://ui.adsabs.harvard.edu/abs/1987Icar...69..532H) and a warning will be raised.
- **Heating & Cooling**: Conduction (warning raised if relevant, **planned**), H3+ line cooling (not planned), Fe & Ca  line cooling (relevant at high Z only, **planned**), free-free cooling (warning raised if relevant, not planned) 
- Multiple ionization states of the same species (**planned**)
See Broome et al. (submitted) for more details.

- Want a rapid H/He model with power-law approximated XUV spectra? Check out [ATES](https://github.com/AndreaCaldiroli/ATES-Code) ([Caldiroli et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..30C/abstract))
- Do you want to set the mass loss rate ($\dot{M}$) yourself or want an EUV isothermal Parker wind outflow model that runs in nanoseconds? Check out [p-winds](https://github.com/ladsantos/p-winds) ([Dos Santos et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A..62D/abstract)).
- Do you want to use p-winds and get transit models for metals via Cloudy? Check out [Sunbather](https://github.com/antonpannekoek/sunbather) ([Linssen et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..43L/abstract))
- Want to leverage Cloudy and the hydrodynamic code PLUTO for more thorough XUV-irradiated, but slightly more expensive calculations? Check out TPCI ([Salz et al. 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...576A..21S/abstract))
- That sound great, but you prefer to code in Python over C/C++? Check out [pyTPCI](https://ascl.net/2506.012). ([Riley, Zhang, & Bean 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...980...34R/abstract))
- Do you care about diffusion throughout the wind? Check out [AIOLIS](https://github.com/Schulik/aiolos) ([Schulik & Booth, 2022](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523..286S/abstract)).
- Want to model the lower atmosphere in more detail? Check out CETIMB (Koskinen et al. 2022)
- Just want a grid of mass loss rates for pure-Hydrogen, low-flux-EUV-irradiated planets? See [Kubyshkina & Fossati](https://ui.adsabs.harvard.edu/abs/2021RNAAS...5...74K/abstract) 
- Want a grid of mass loss rates for pure-Hydrogen, high-flux-XUV-irradiated planets? See [Owen & Jackson (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.2931O/abstract)

>Want your model added to this list or to update the short bio? Email mabroome@ucsc.edu

Requirements
------------

`Wind-AE` requires the following packages and will pip install them automatically by following the Installation guide below.

* `python`>3.13
* `numpy` 
* `scipy`
* `astropy`
* `pandas`>=2.2.3
* `matplotlib` 
* `datetime`
* `pyarrow` 
* `fastparquet`
* `requests` 
* `ChiantiPy`

Installation
------------
Until `Wind-AE` is dockerized, it is recommended to use a Python environment to avoid dependency issues. However, if your system meets the above requirements, there is no need to create an environment and you can skip to the compilation step.

To create an environment use either
```angular2html
python3 -m venv venv_name.venv
source venv_name.venv/bin/activate
```
or using `conda`
```angular2html
conda create -n venv_name
conda activate venv_name
conda install pip
```
### Pip install

Recommended:
```angular2html
pip install --upgrade pip
```

Then 
```angular2html
pip install wind_ae
```


### OR Compile from source (BETA)

Clone the repository using
```angular2html
git clone https://github.com/mibroome/wind-ae/
```
or navigate to [github.com/mibroome/wind-ae/](https://github.com/mibroome/wind-ae/) and download and unzip the zip file.

To compile from the source,
```angular2html
pip install -r requirements.txt
pip install -e .
```

### Confirming the import was successful

Run tests (optional). Estimated time: 4 minutes.
```angular2html
cd wind-ae && pytest
```

Otherwise, you can test the install by running
```angular2html
python -c "import wind_ae"
```

Now you can run `Wind-AE` from anywhere! As seen in the tutorial (`Notebooks/Quickstart.ipynb`), the following imports are helpful for most purposes. 
```angular2html
from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
from wind_ae.wrapper.wrapper_utils.plots import energy_plot six_panel_plot quick_plot 
from wind_ae.wrapper.wrapper_utils import constants as const
from wind_ae.wrapper.wrapper_utils.system import system
from wind_ae.wrapper.wrapper_utils.spectrum import spectrum
```

> **Note**: If you ever need to interface directly with the `C` code, it lives in `wind_ae/src/` and can be excecuted from within the `wind_ae/` folder via `./bin/relaxed_ae`. The solution generated will be for a planet with the parameters detailed in the input files in the `Inputs/` folder. There is generally no need to interface with the `C` code and most standard tasks can be accomplished by using the Python wrapper.

## Future features and known problems:

- Computation of the complementary error function that governs the drop-off of bolometric heating/cooling is not truly self-consistent (`converge_mol_atomic_transition(polish=True, width=)`) and may require visual confirmation via `energy_plot()` (checking whether bolometric heating/cooling impede too far into photoionization heating or fall too short) and manual adjustment of the `width` parameter or:

```python
sim.load_planet('path/to/planet/file')
bcs = np.copy(sim.windsoln.bcs_tuple)
# erf_loc - normalized velocity value at radius where you want the erf to drop
# erf_rate - how quickly the erf drops off in units of Hsc at erf_loc
# To get initial estimate, run sim.erf_velocity(polish=True)
bcs[-1] = np.array([erf_loc, erf_rate])
sim.inputs.write_bcs(*bcs)
sim.run_wind()
```

- Knudsen number calculations currently only contain H-H collisions.
- Converting spectrum ``kind`` from ``'full'`` to ``'mono'`` occasionally has issues.

--------
### Check out the [open issues](https://github.com/mabroome/wind-ae/issues).