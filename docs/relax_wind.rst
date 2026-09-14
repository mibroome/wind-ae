.. _relax_wind:

Relaxation Module
=====================
For examples of how to use these functions, see `Quickstart <Quickstart.ipynb>`_.

.. automodule:: wind_ae.wrapper.relax_wrapper

Working directories & running simultaneous simulations
------------------------------------------------------
.. versionadded:: 2.0

Every instantiation of ``sim = wind_sim()`` creates its own isolated temporary working
directory (in your machine's ``tmp`` folder) containing its own ``inputs/``, ``saves/``,
and ``outputs/`` subfolders. This means multiple ``wind_simulation`` instances --- in
separate notebook cells, threads, or processes --- can be run **simultaneously** without
overwriting one another's input/output files, which was not possible in Wind-AE 1.x.

.. code-block:: python

   from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim

   sim1 = wind_sim()
   sim2 = wind_sim()  # independent workdir; safe to run alongside sim1

   sim1.load_planet('saves/planet_a.csv')
   sim2.load_planet('saves/planet_b.csv')
   # sim1.run_wind() and sim2.run_wind() can now be called concurrently
   # (e.g. from separate threads/processes/Jupyter notebooks) without input/output collisions.

The path to an instance's working directory is available via ``sim.workdir`` (see
:meth:`~wind_ae.wrapper.relax_wrapper.wind_simulation.where`). To execute the compiled C
relaxation code directly on a given instance's inputs, ``cd`` into that workdir and run
``/Users/your/path/to/wind-ae/wind_ae/bin/relaxed_ae`` from *within* it, e.g.:

.. code-block:: bash

   cd $(python -c "from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim; print(wind_sim().workdir)")
   /Users/your/path/to/wind-ae/wind_ae/bin/relaxed_ae

Temporary workdirs are cleaned up automatically when your machine restarts, or on demand
via :meth:`~wind_ae.wrapper.relax_wrapper.wind_simulation.cleanup`. ``wind_simulation``
also supports the context-manager protocol, which calls ``cleanup()`` for you on exit:

.. code-block:: python

   with wind_sim() as sim:
       sim.load_planet('saves/seed.csv')
       sim.run_wind()
       sim.save_planet('saves/output.csv')
   # sim.workdir has now been removed

If you pass an explicit ``workdir=`` to ``wind_sim()``, Wind-AE will use that directory
instead of creating a temporary one, and will **not** delete it on ``cleanup()`` --- you
are responsible for its lifetime in that case.

.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.where
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.cleanup

Loading and interacting with wind solutions
------------------------------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.load_planet
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.run_wind
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.direct_solve

.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.load_uservars
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.generate_rate_coeffs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.load_spectrum
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.save_planet
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.easy_output_file

Ramping functions
------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_to
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_Ftot
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_var
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_grav
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_star

Metals functions
-------------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.calc_metallicity
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.add_metals
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.remove_metals
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_metallicity

Heating & Cooling Flags
-------------------------------
.. versionadded:: 2.0
   Conductive cooling, recombination cooling, and free-free (bremsstrahlung) cooling are
   now independently toggleable, alongside line cooling and the molecular/bolometric
   layer. See ``sim.windsoln.flags_tuple`` and :ref:`windsoln` for the full flag list.

.. warning::
   Turning on conduction makes the ODEs very stiff. Only turn it on **after** all other
   ramping is complete (:meth:`~wind_ae.wrapper.relax_wrapper.wind_simulation.turn_on_conduction`).
   The temperature dependence of free-free and recombination cooling can make Wind-AE
   numerically unstable, and ramping to higher metallicities with these on may fail --- if
   so, try turning them off during ramping and back on afterwards.

.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_on_conduction
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_off_conduction
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_off_line_cool
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_on_molecular_layer
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_off_molecular_layer

Boundary Conditions 
-----------------------------------------------
Polishing (Converging to self-consistent BCs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.polish_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_base_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.converge_mol_atomic_transition
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_molecular_erfc
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.converge_Ncol_sp
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.converge_Rmax
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.integrate_out

Computing self consistent boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.find_base_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.find_self_consistent_Ncol
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_molec_adjust
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.raise_Ncol_sp

Explicitly Ramping Lower Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_T_rmin
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_Rmin
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_rho_rmin
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_off_tidal_grav

Ramping spectrum
-------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.flux_norm
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_spectrum
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_to_user_spectrum
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.format_user_spectrum

Plotting Aliases
-------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.energy_plot
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.six_panel_plot
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.quick_plot
