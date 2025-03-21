# Wind-AE (BETA)

`Wind-AE` (/windy/) stands for "wind atmospheric escape" and is a relatively fast 1D, steady-state, hydrodynamic, non-isothermal, Parker wind relaxation code for modeling atmospheric escape based on [Murray-Clay et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...693...23M/abstract). `Wind-AE` is a forward model that solves the energy conservation, momentum conservation, and ionization equilibrium equations at the substellar point using the ["Numerical Recipes in C"](https://ui.adsabs.harvard.edu/abs/1992nrca.book.....P/abstract) relaxation method. This allows `Wind-AE` to quickly compute the atmosphere mass loss rate as well as upper atmosphere (> ~100 microbar) velocity, temperature, density, and ionization structure as a function of altitude. 

`Wind-AE` updates [Murray-Clay et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...693...23M/abstract) to allow for the modeling of atomic metals and multifrequency XUV stellar spectra (Broome et al. submitted). If you use `Wind-AE`, please consider citing Broome et al. (in prep). 


We appreciate your patience while the docs are developed. In the meantime, take advantage of `Notebooks/Quickstart.ipynb` to get a quick orientation to `Wind-AE` and please report any bugs via [Github](https://github.com/mabroome/wind-ae/issues) or via email to mabroome@ucsc.edu.

Is Wind-AE the right tool for me?
----------------------
`Wind-AE` is well-suited for users interested in **quickly estimating mass loss rates** or outflow structure. Outflow structure includes bulk temperature and per-species ionization fractions as a function of radius, so can be easily translated into approximating and **predicting observables and transits**, including metastable helium (He 10830 Å) transits, though a He transit module is not yet included. Precise modeling of lower atmosphere (< ~100 microbar) is considered necessary for highly accurate transit models, but Wind-AE can be easily coupled to lower atmosphere photochemistry models whose outputs can (e.g., radius, temperature, abundances, ionization fractions, etc. at 1 microbar) can be fed into Wind-AE as inputs.  

>**If you are interested in outflow structure:** Past the Coriolis turning radius (a few planetary radii) 3D physics dominates, so `Wind-AE` does not integrate past that point. `Wind-AE` also makes simplifying assumptions about the region below the region below the wind-launch radius (~10 nanobars). 

Because `Wind-AE` runs on the order of seconds to minutes, it also can be (and has been) used to **model planet evolution**.

#### `Wind-AE` can model:
- Multiple atomic species
- X-ray physics (secondary ionizations and K-shell ionization cross-sections for relevant metals)
- Both low and high stellar XUV flux
- **Heating & Cooling**: Ionization heating, bolometric heating & cooling (negligible in wind), PdV cooling (work done due to expansion of gas), radiative / atomic line cooling (Lyman-$\alpha$, OI, OII, OIII, CI, CII), recombination cooling
#### `Wind-AE` does not (currently) include:
- **Magnetic fields**
- **Time dependence**
- **Diffusion/drag** - the atomic species set by the user are assumed to be entrained in the outflow and in thermal equilibrium. This is an appropriate assumption for species below the [crossover mass](https://ui.adsabs.harvard.edu/abs/1987Icar...69..532H) and a warning will be raised.
- **Heating & Cooling**: Conduction (warning raised if relevant, **planned**), H3+ line cooling (not planned), Fe & Ca  line cooling (relevant at high Z only, **planned**), free-free cooling (warning raised if relevant, not planned).
- Multiple ionization states of the same species (**planned**) - i.e., currently only OII->OIII not OII->OIII->OIV can be modeled.

### Models for other use cases
- Want a rapid H/He model with power-law approximated XUV spectra? Check out [ATES](https://github.com/AndreaCaldiroli/ATES-Code) ([Caldiroli et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..30C/abstract))
- Do you want to set the mass loss rate ($\dot{M}$) yourself or want an EUV isothermal Parker wind outflow model that runs in nanoseconds? Check out [p-winds](https://github.com/ladsantos/p-winds) ([Dos Santos et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...659A..62D/abstract)).
- Do you want to use p-winds and get transit models for metals via Cloudy? Check out [Sunbather](https://github.com/antonpannekoek/sunbather) ([Linssen et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..43L/abstract))
- Want to leverage Cloudy and the hydrodynamic code PLUTO for more thorough XUV-irradiated, but slightly more expensive calculations? Check out TCPI ([Salz et al. 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...576A..21S/abstract))
- Do you care about diffusion throughout the wind? Check out [AIOLIS](https://github.com/Schulik/aiolos) ([Schulik & Booth, 2022](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523..286S/abstract)).
- Want to model the lower atmosphere in more detail? Check out CETIMB (Koskinen et al. 2022)
- Just want a grid of mass loss rates for pure-Hydrogen, low-flux-EUV-irradiated planets? See [Kubyshkina & Fossati](https://ui.adsabs.harvard.edu/abs/2021RNAAS...5...74K/abstract) 
- Want a grid of mass loss rates for pure-Hydrogen, high-flux-XUV-irradiated planets? See [Owen & Jackson (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.2931O/abstract)

>Want your model added to this list or to update the short bio given above? Email mabroome@ucsc.edu

Requirements
------------

`Wind-AE` requires the following packages and will pip install them automatically by following the Installation guide below.

* `python`>3.9.12 
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
### Compile from source (only available in beta)

Clone the repository using
```angular2html
git clone https://github.com/mibroome/wind-ae/
```
or navigate to [github.com/mibroome/wind-ae/](https://github.com/mibroome/wind-ae/) and download and unzip the zip file.

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

To compile as an editable module from the source,
```angular2html
pip install --upgrade pip
python3 -m pip install -e .
```
(`pip` upgrade is recommended as `pip<25.0.1` has had issues installing dependencies. Deprecation warning about setup.py file can be ignored, fix in progress.)

<!-- REMEMBER to `make` the C code before running `Wind-AE` for the first time.
```angular2html
cd wind_ae/wind_ae/
make
```
You can ignore any warnings that pop up. -->

You can test the install by running
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

Future features and known problems
--------
Check out the [open issues](https://github.com/mabroome/wind-ae/issues).