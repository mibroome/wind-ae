.. _plots:

Plotting
=====================
All plotting functions take windsoln objects as inputs.

.. code-block:: python
   
   from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim
   from wind_ae.wrapper.wrapper_utils.plots import spectrum_plot

   sim = wind_sim()
   sim.load_planet('planet.csv')
   spectrum_plot(sim.windsoln,)

The two most common plots for analysis (``energy_plot()`` and ``six_panel_plot()``) are 
also accessible via the functions in :ref:`relax_wind` as 
:meth:`wind_ae.wrapper.relax_wind.six_panel_plot()` and :meth:`wind_ae.wrapper.relax_wind.energy_plot()` (e.g., ``sim.six_panel_plot()``).

.. automodule:: wind_ae.wrapper.wrapper_utils.plots
   :members:
   :undoc-members:
   :show-inheritance: