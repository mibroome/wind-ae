.. _installation:

Installation
============

Requirements
______________

**Wind-AE** requires the following packages and will pip install them automatically by following the Installation guide below.

- python >= 3.13.0
- numpy
- scipy
- astropy
- pandas >= 2.2.3
- matplotlib
- datetime
- pyarrow
- fastparquet
- requests
- ChiantiPy

Installation Instructions
___________________________

Until ``Wind-AE`` is dockerized, it is recommended to use a Python environment to avoid dependency issues. However, if your system meets the above requirements, there is no need to create an environment and you can skip to the compilation step.

To create an environment use either

.. code-block:: bash

	python3 -m venv venv_name.venv
	source venv_name.venv/bin/activate

or using ``conda``

.. code-block:: bash

    conda create -n venv_name
    conda activate venv_name
    conda install pip

Pip install
-----------

Recommended to upgrade pip first:

.. code-block:: bash

    pip install --upgrade pip

Then

.. code-block:: bash

    pip install wind_ae

OR Compile from source (BETA)
-------------------------------

Clone the repository using

.. code-block:: bash

	git clone https://github.com/mibroome/wind-ae/

or navigate to `github.com/mibroome/wind-ae/ <https://github.com/mibroome/wind-ae/>`_ and download and unzip the zip file.

To compile from the source

.. code-block:: bash

    pip install -r requirements.txt
    pip install -e .

Confirming the import was successful
---------------------------------------

Run tests (optional). Estimated time: 4 minutes

.. code-block:: bash

    cd wind-ae && pytest

Otherwise, you can test the install by running

.. code-block:: bash

    python -c "import wind_ae"


Now you can run **Wind-AE** from anywhere! As seen in the tutorial, the following imports are helpful for most purposes:

.. code-block:: python

	from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
	from wind_ae.wrapper.wrapper_utils.plots import energy_plot, six_panel_plot, quick_plot
	from wind_ae.wrapper.wrapper_utils import constants as const
	from wind_ae.wrapper.wrapper_utils.system import system
	from wind_ae.wrapper.wrapper_utils.spectrum import spectrum

.. note::

    If you ever need to interface directly with the C code, it lives in ``wind_ae/src/`` and is compiled to
    ``wind_ae/bin/relaxed_ae``. The solution generated will be for a planet with the parameters detailed in the
    input files in the ``inputs/`` folder that ``relaxed_ae`` is run from. There is generally no need to interface
    with the C code and most standard tasks can be accomplished by using the Python wrapper.

    .. versionchanged:: 2.0
        Each ``sim = wind_sim()`` instance now writes its input files to its own isolated temporary working
        directory (``sim.workdir``, see :ref:`relax_wind`) instead of the installed ``wind_ae/`` folder, so that
        multiple simulations can run simultaneously without overwriting one another's inputs/outputs. To run
        ``relaxed_ae`` directly against a given instance's inputs, ``cd`` to *that instance's* workdir first:

        .. code-block:: bash

            cd sim.workdir       # e.g. print(sim.workdir) from Python first
            /Users/your/path/to/wind-ae/wind_ae/bin/relaxed_ae

Future features and known problems
___________________________________
- Computation of the complementary error function that governs the drop off of bolometric heating/cooling is not truly self-consistent (``converge_mol_atomic_transition(polish=True, width_factor=)``) and may require visual confirmation via ``energy_plot()`` (checking whether bolometric heating/cooling impede too far into photoionization heating or fall too short) and manual adjustment of the ``width_factor`` parameter:

.. code-block:: python

	sim.load_planet('path/to/planet/file')
	# width_factor widens the transition region (in scaleheights); start at 0 and
	# increase if energy_plot() shows the transition impeding too far into the
	# photoionization-heated region.
	sim.converge_mol_atomic_transition(polish=True, width_factor=0)

For full manual control over the erfc drop-off location/rate (rarely necessary), use
:meth:`~wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_molecular_erfc`:

.. code-block:: python

	sim.load_planet('path/to/planet/file')
	# erf_loc  - normalized velocity value at radius where you want the erf to drop
	# erf_rate - how quickly the erf drops off in units of Hsc at erf_loc
	# The estimator that converge_mol_atomic_transition() uses internally is available
	# as sim._erf_velocity(polish=True) if you need a starting estimate, but as an
	# underscore-prefixed method it is not a stable, supported part of the public API.
	sim.ramp_molecular_erfc(v_drop=erf_loc, rate=erf_rate)

- Knudsen number calculations currently only contain H-H collisions.
- Converting spectrum ``kind`` from ``'mono'`` to ``'full'`` occasionally has issues. Converting ``'full'`` to ``'mono'`` has no issues.
- **[v2.0]** Turning on conduction (``sim.turn_on_conduction()``) introduces numerically
  stiff ODEs; only do so after all other ramping is complete, and expect longer runtimes.
  Ramping to high metallicity with recombination and/or free-free cooling on can also be
  numerically unstable --- if a ramp fails, try turning those flags off during ramping and
  back on afterwards.

Check out the `open issues <https://github.com/mibroome/wind-ae/issues>`_.
