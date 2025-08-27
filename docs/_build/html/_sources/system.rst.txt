.. _system:

System Module
=====================
Mostly for internal use, however it is useful for printing system parameters in a human readable format.

.. code-block:: python
   
   from wind_ae.wrapper.relax_wrapper import wind_simulation as wind_sim

   sim = wind_sim()
   sim.load_planet('planet.csv')
   sim.system.print_system()

.. automodule:: wind_ae.wrapper.wrapper_utils.system
   :members:
   :undoc-members:
   :show-inheritance: