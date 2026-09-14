# Wind-AE v2.0 Release Notes

## Highlights

- **Simultaneous / parallel runs.** `wind_sim()` instances no longer share a single `inputs/`/`saves/` directory, so multiple simulations can now be run at the same time (e.g. across notebook cells, threads, or processes) without overwriting each other's files.
- **Multiple ionization states of the same species.** Wind-AE can now model, e.g., C I and C II simultaneously as a linked ionization chain — previously only one ionization state per element could be tracked at a time. See [New: multiple ionization states (species chains)](#new-multiple-ionization-states-species-chains) below.
- **Major speed-up for metals.** Adding metals and ramping to 10Z or 100Z used to take ~30-70 minutes; it now takes ~3 minutes (see [Performance](#performance) below for why).
- **Electron-counting bug fixes.** Errors in electron counting that affected recombination cooling, mean molecular/atomic weight (mu), and other derived quantities have been resolved.
- **New, independently toggleable cooling physics.** Conductive cooling, recombination cooling, and free-free (bremsstrahlung) cooling are now available as flags.

---

## Simultaneous runs & isolated working directories

Every instantiation of `sim = wind_sim()` now creates its own temporary working directory on your machine's `tmp` folder (`sim.workdir`), containing that instance's own `inputs/`, `saves/`, and `outputs/`. This eliminates the input/output overwrite conflicts that previously made it unsafe to run more than one Wind-AE simulation at a time.

```python
sim1 = wind_sim()
sim2 = wind_sim()  # fully independent from sim1 — safe to run concurrently
```

- Temp workdirs are deleted automatically when your machine restarts, or on demand via `sim.cleanup()`.
- `wind_simulation` now supports the context-manager protocol, which calls `cleanup()` for you: `with wind_sim() as sim: ...`.
- The new `sim.where()` prints the workdir path and its current contents.
- To run the compiled C binary directly against a given instance's inputs, `cd` into `sim.workdir` (not the installed `wind_ae/` package folder) and run `/path/to/wind-ae/wind_ae/bin/relaxed_ae` from there.
- Pass `workdir=<path>` to `wind_sim()` to opt out of the auto-created temp directory and manage the directory's lifetime yourself (`cleanup()` will not touch it).

---

## New: multiple ionization states (species chains)

Wind-AE can now track more than one ionization state of the same element in a single simulation — e.g. C I *and* C II simultaneously, linked as a chain — which was not possible in v1.x (each metal could only have one ionization state modeled at a time). This is implemented at the C-solver level (`chain_parent` bookkeeping in `io.c`, `soe.c`, `difeq.c`), not just in the Python wrapper, so recombination, advection, and column-density calculations correctly cascade through the chain (e.g. C II's population is drawn from C I's ionized fraction).

To use it, just include both states in your species list — chain linkage (same element, one fewer electron) is detected automatically:

```python
sim.add_metals(['CI', 'CII'])
```

Chain parent/child species must share the same mass fraction (`HX`) in `phys_params.inp`; `add_metals()` handles this for you automatically.

---

## Performance

Adding metals and ramping metallicity used to require rewriting `src/rate_coeffs.h` and recompiling the C source on every step. Rate coefficients are now written to `inputs/rate_coeffs.inp` and read at runtime, so metals steps no longer trigger a recompilation — this is the main driver of the ~1-2 hour → ~3 minute speed-up for metals + ramping to 10Z/100Z.

---

## New heating & cooling flags

`sim.windsoln.flags_tuple` has grown from 4 flags to 8. In order:

`integrate_outward, tidalforce, linecool, bolo_heat_cool, conduction, recombo_cool, free_free_cool, molec_layer`

| Flag | Status | Notes |
|---|---|---|
| `linecool` | Renamed from `lyacool` | Governs Lyman-α **and** metal line cooling |
| `conduction` | **New**, off by default | Self-consistent conductive heat flux. See warning below. |
| `recombo_cool` | **New**, on by default | Recombination cooling (`cool_rec`) |
| `free_free_cool` | **New**, on by default | Free-free / bremsstrahlung cooling (`cool_free`) — previously not implemented at all |
| `molec_layer` | **New** | Separate erfc multiplier for the mean-molecular-weight (mu) molecular-to-atomic transition, independent of `bolo_heat_cool`. Defaults to `bolo_heat_cool`'s value when loading older solution files. |

New convenience methods: `sim.turn_on_conduction()` / `sim.turn_off_conduction()`, `sim.turn_off_line_cool()`, and `sim.turn_on_molecular_layer(kind=)` / `sim.turn_off_molecular_layer(kind=)` (`kind='bolo_only'`, `'molec_only'`, or `'both'`) — replacing the old `turn_on_bolo()`/`turn_off_bolo()`.

> **⚠️ Warning:** Conductive cooling makes the ODEs very stiff and should only be turned on after all other ramping is complete. The temperature dependence of free-free and recombination cooling makes Wind-AE numerically unstable at high metallicity in some cases — if a ramp fails, try turning these flags off during ramping and back on afterward.

**Also new:** metal line cooling now covers Fe, Mg, Ca, and Ne (via Fe II, Mg II, Ca II, and Ne III transitions), in addition to the existing O and C lines — previously only planned/not implemented.

**Also new:** the optical (`kappa_opt`) and IR (`kappa_IR`) opacities used in the molecular/bolometric layer, and the adiabatic index (`gamma`), are now first-class per-planet parameters stored on the solution and tunable via `write_physics_params()` / the `physics` object, rather than fixed constants.

---

## Renamed / changed APIs

| Old (v1.x) | New (v2.0) | Notes |
|---|---|---|
| `sim.base_bcs()` | `sim.find_base_bcs()` | Same behavior; renamed to distinguish "find" from "ramp"/"converge" functions |
| `sim.self_consistent_Ncol()` | `sim.find_self_consistent_Ncol()` | Same behavior |
| `sim.turn_on_bolo()` / `sim.turn_off_bolo()` | `sim.turn_on_molecular_layer()` / `sim.turn_off_molecular_layer()` | Generalized with a `kind=` argument to independently control bolometric heating/cooling and/or the mu transition |
| `sim.erf_velocity()` | `sim._erf_velocity()` | Now a private/internal helper. Use `sim.converge_mol_atomic_transition(polish=True, width_factor=...)` instead, or `sim.ramp_molecular_erfc(v_drop, rate)` for full manual control |
| `sim.ramp_spectrum(..., kind='full')` | `sim.ramp_spectrum(...)` | `kind` argument removed — mono vs. full is now inferred automatically from `goal_spec_range` |
| `sim.ramp_to_user_spectrum(..., updated_F=...)` | `sim.ramp_to_user_spectrum(..., updated_Ftot=...)` | Renamed for consistency |
| `wind_sim(csv_file=...)` | `wind_sim(csv_file=None, workdir=None)` | New `workdir` parameter (see above); `csv_file` default changed to `None` |
| `input_handler()` | `input_handler(workdir=None)` | New parameter, ties into the workdir system |
| `metal_class(windsoln_object)` | `metal_class(windsoln_object, workdir=None)` | New parameter |
| `write_flags(lyacool, tidalforce, bolo_heat_cool, integrate_outward, ...)` | `write_flags(integrate_outward, tidalforce, linecool, bolo_heat_cool, conduction, recombo_cool=1, free_free_cool=1, molec_layer=None, ...)` | Reordered and extended — see flags table above |
| `write_physics_params(..., phys_file=<hardcoded path>)` | `write_physics_params(..., kappa_opt=4e-3, kappa_IR=1e-2, gamma=5/3, phys_file=None)` | New physics parameters; `phys_file` now honors `workdir` when `None` |
| `energy_plot(..., all_terms=, CII_line_cool=, CIII_line_cool=, OII_line_cool=, OIII_line_cool=, sub_sonic=)` | `energy_plot(..., plot_dom_lines=True, N_top_lines=2)` | Simplified: automatically plots the dominant cooling term(s) at each of several sample radii instead of requiring each line species to be toggled manually |
| `six_panel_plot(..., first_plotted=True, ax=0)` / `quick_plot(..., first_plotted=True, ax=0)` | `six_panel_plot(..., ax=None)` / `quick_plot(..., ax=None)` | `first_plotted` removed — whether this is the first plot on a shared axis is now inferred from `ax is None` |
| `regrid(q_arr=None, simple=True)` | `regrid(dq=8e-3, write_to_guess=False, q_arr=None, simple=True)` | New `dq` (resolution) and `write_to_guess` (write result directly to `saves/windsoln.csv` as the next guess) parameters |
| `ramp_base_bcs(Kappa_opt=0.004, Kappa_IR=0.01, ...)` | `ramp_base_bcs(user_Rmin=None, user_rho_Rmin=None, user_T_rmin=None, Kappa_opt=None, Kappa_IR=None, ...)` | New `user_Rmin`/`user_rho_Rmin`/`user_T_rmin` let you directly override the computed base boundary conditions; `Kappa_opt`/`Kappa_IR` now default to the current solution's values |
| `polish_bcs(...)` | `polish_bcs(..., add_conduction=False)` | New flag to automatically add conduction if it's found to be significant (>70% of photoionization heating) |

---

## New functions

- `sim.direct_solve(system)` — writes inputs directly from a target `system` and re-solves without any intermediate ramping. Useful (though not always successful, given the relaxation method's sensitivity) for re-solving after a small parameter change.
- `sim.ramp_molec_adjust(goal)` — ramps the mean molecular weight adjustment factor (`molec_adjust`) to a goal value independently.
- `input_handler.write_solver_params(M, itmax, rhoscale)` — writes solver parameters (previously baked into other write calls).
- `spectrum.write_mono_spectrum(mono_nm, species, savefile=...)` — writes a monofrequency spectrum file directly, bypassing the full smoothing/binning pipeline (faster for monofrequency-only workflows).

## Newly documented (pre-existing but previously undocumented)

The following public methods already existed but were missing from the docs; they've now been added in their respective sections: `wind_solution.add_user_vars()`, `wind_solution.current_metallicity()`, `wind_solution.alpha_rec()`, `wind_solution.calc_roche_lobe()`, and `wind_simulation.ramp_molecular_erfc()`.

---

## Other under-the-hood fixes

- **Chain-aware metal addition.** `add_metals()` now correctly seeds the mass fraction of a newly added ionization state (e.g. adding C II when C I is already present) from its parent species' mass fraction, rather than always starting from a near-zero seed. When adding multiple species at once, each new species' column density is now seeded independently from the same pre-existing column, rather than chaining off each other's already-tiny seed values.
- Recombination cooling, and the associated `recomb_*`/`advec_*` rate columns, now use a chain-aware effective-density calculation consistent with the C solver's `soe.c`, which is the fix underlying the electron-counting corrections mentioned above.

---

## Known issues
- Ramping to high metallicity with recombination and/or free-free cooling on can be numerically unstable in some cases (see warning above).
- Knudsen number calculations still only include H-H collisions (multispecies planned).
- Converting spectrum ``kind`` from ``'mono'`` (monofrequency) to ``'full'`` (multifrequency) occasionally has issues. There is no trouble going from ``'full'`` to ``'mono'``.