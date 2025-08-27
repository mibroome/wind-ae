.. _relax_wind:

Relaxation Module
=====================
For examples of how to use these functions, see `Quickstart <Quickstart.ipynb>`_.

.. automodule:: wind_ae.wrapper.relax_wrapper

Loading and interacting with wind solutions
------------------------------------------------
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.load_planet
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.run_wind

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
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.metallicity
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.add_metals
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.remove_metals
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_metallicity

Boundary Conditions 
-----------------------------------------------
Polishing (Converging to self-consistent BCs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.polish_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.ramp_base_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.run_isotherm
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.converge_Ncol_sp
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.converge_Rmax
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.integrate_out

Computing self consistent boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.base_bcs
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.self_consistent_Ncol
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.turn_off_bolo
.. automethod:: wind_ae.wrapper.relax_wrapper.wind_simulation.erf_velocity
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
