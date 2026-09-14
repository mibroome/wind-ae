.. _metals:

Metals Module
=====================
Typically called via the functions in :ref:`relax_wind` (`sim.add_metals`, `sim.calc_metallicity`, etc.)

.. note::
   .. versionadded:: 2.0
      Multiple ionization states of the same element can now be modeled simultaneously
      as a linked ionization chain (e.g. ``sim.add_metals(['CI','CII'])`` tracks C I and
      C II together, with C II's population drawn self-consistently from C I's neutral
      fraction). Chain linkage is detected automatically from ``species_list`` --- any two
      species of the same element where one has exactly one fewer electron than the other
      are treated as parent/child. When a chain is present, the parent and child must be
      given the **same** mass fraction in ``phys_params.inp`` (handled automatically when
      going through ``add_metals()``).

.. automodule:: wind_ae.wrapper.wrapper_utils.metals
   :members:
   :undoc-members:
   :show-inheritance: