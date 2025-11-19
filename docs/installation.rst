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
________________

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

    If you ever need to interface directly with the C code, it lives in ``wind_ae/src/`` and can be executed from 
    within the ``wind_ae/`` folder via ``./bin/relaxed_ae``. The solution generated will be for a planet with the 
    parameters detailed in the input files in the ``Inputs/`` folder. There is generally no need to interface with the 
    C code and most standard tasks can be accomplished by using the Python wrapper.

Future features and known problems
___________________________________
- Computation of the complementary error function that governs the drop off of bolometric heating/cooling is not truly self-consistent (``converge_mol_atomic_transition(polish=True,width=)``) and may require visual confirmation via ``energy_plot()`` (checking whether bolometric heating/cooling impede too far into photoionization heating or fall too short) and manual adjustment of the ``width`` parameter 

.. code-block:: python

	sim.load_planet('path/to/planet/file')
	bcs = np.copy(sim.windsoln.bcs_tuple)
	# erf_loc - normalized velocity value at radius where you want the erf to drop
	# erf_rate - how quickly the erf drops off in units of Hsc at erf_loc
	# To get initial estimate, run sim.erf_velocity(polish=True)
	bcs[-1] = np.array([erf_loc,erf_rate])
	sim.inputs.write_bcs(*bcs)
	sim.run_wind()

- Knudsen number calculations currently only contain H-H collisions.
- Converting spectrum ``kind`` from ``'full'`` to ``'mono'`` occasionally has issues.

Check out the `open issues <https://github.com/mibroome/wind-ae/issues>`_.
