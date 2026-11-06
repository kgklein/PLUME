> Historical implementation report, archived before the Stage 5B.1 manual.
> Use [Normalization](../normalization.md) for the current user reference.
>
> The source below originally lived at repository-root normalization.md; its
> relative links retain that historical location. See the
> [original at the archived revision](https://github.com/kgklein/PLUME/blob/a853ca0ed4882e9738e4fcc29b8efb1cf4f64013/normalization.md)
> or the current manual's developer/history index for navigable current references.

# Dual normalization: source audit and implementation report

Historical report. [Stage 3B.2](tests/COLD_REFERENCE_MIXED_REPORT.md) now supports cold-reference
mixtures containing warm species; historical measurements below are retained.

This report records the original dual-normalization stage. Stage 2 now permits
non-reference cold species in collisionless `d_ref` dispersion calculations;
the adapter skips their zero temperatures. See the current [input guide](input.md)
and [mixed-dispersion report](tests/MIXED_DISPERSION_REPORT.md) for supported
functionality. The original validation results below remain historical records.

Stage 3A adds a fully cold path that bypasses the adapter altogether, using
`d_ref` with `betap=0`. Its scope and tests are recorded in the
[fully cold report](tests/FULLY_COLD_DISPERSION_REPORT.md); cold-reference
mixtures containing warm species remain deferred to Stage 3B.

PLUME accepts `rho_ref` and `d_ref` inputs while retaining the finite-temperature
kinetic calculation in its existing variables. The parameter definitions and a
matched input pair are documented in [README.md](README.md) and [input.md](input.md).

## Source audit

| Subsystem | Findings |
| --- | --- |
| Parameter storage and input | `src/vars.f90` declares global parameters, the `specie` type, and scan state. `read_in_params` and `spec_read` in `src/functions.f90` read `&params` and indexed species namelists. Validation precedes kinetic evaluation. |
| Reference species | Species 1 supplies the reference anisotropy and reference frame checks. The Python wrapper likewise assumes species 1 is the reference. No new reference selector was introduced. |
| Susceptibility | `disp0` and `disp1` in `src/disprels.f90` implement collisionless and weakly collisional kinetic responses. Both use `betap/vtp**2` for the plasma-frequency scale, signed `muS/Qs` for gyrofrequency, thermal Bessel arguments, and thermal resonance arguments. |
| Dispersion and roots | `disp` dispatches the susceptibility calculation. `map_search`, `refine_guess`, and `rtsec` use dimensionless frequency `omega/Omega_R`. The determinant and susceptibility are identical for physically equivalent inputs. |
| Eigenfunctions and heating | `calc_eigen` uses the wave tensor, susceptibility, and derived thermal state for fields, density, velocity, total heating, and seven resonance heating channels. Its finite-frequency derivative calculation is unchanged. |
| Diagnostics | `src/fpc.f90` uses thermal velocities for gyrotropic and Cartesian distributions and correlations. Its file headers describe derived legacy parameters even for a `d_ref` input. The public single-frequency diagnostic uses the same derived state. |
| Scans | `om_scan`, `om_double_scan`, and `map_scan` mutate input parameters through scan pointers. Their setup helpers cover both wavevectors, magnitude/angle, beta, global velocity, and all six species parameters. `read_radial_input` reads `.rad` tables; `radial_scan` interpolates their values and continues roots through wavevector scans. |
| Python | `linfpclib/linfpc.py` constructs and reads namelists, computes collision frequencies, and loads numeric outputs. `contours/plot_plume_resonances.py` interprets temperature widths, drift/resonance velocities, and wavevector labels. |

The existing gyroradius is **perpendicular**:
`rho_R = v_tperp,R / abs(Omega_R)`. For an anisotropic reference,

\[
k d_R=\frac{k\rho_R}{\sqrt{\beta_{\parallel,R}\alpha_R}},\qquad
\alpha_R=T_{\perp,R}/T_{\parallel,R}.
\]

Omitting the anisotropy factor changes the physical system when `alphS` for
species 1 differs from unity. The sign of `Qs` is retained in cyclotron terms;
plasma-frequency magnitudes use `Qs**2`.

Temperature-sensitive expressions include the Bessel argument proportional to
`1/tauS`, the thermal resonance argument, the drift factor
`kpar*vvS/sqrt(betap*alpha_R)`, eigenfunction factors `vtp**2/betap`, and FPC
factors containing `sqrt(tauS/(muS*betap))`. None is a cold susceptibility.
Zero and negative kinetic temperatures are rejected before evaluating them.

## Architecture and scan semantics

Approach A was selected for this stage. It requires fewer changes than
rewriting both susceptibility kernels and their diagnostics in inertial units,
and retains the legacy floating-point expressions. Approach B would avoid a
thermal adapter at very small beta and offer a more direct cold-reference path,
but would change substantially more numerical code.

The public globals remain the authoritative input state in the chosen
convention. `sync_solver_state` derives `solver_*` values before every `disp`
evaluation. Susceptibility, eigenfunction, and FPC routines read that derived
state. In `rho_ref`, the adapter copies inputs without changing their values.
In `d_ref`, it uses

\[
k\rho_R=(kd_R)\sqrt{\beta_R\alpha_R},\qquad
v_{t\parallel,R}/c=\sqrt{\beta_R}\,(v_{A,R}/c),\qquad
\tau_s^{(\rho)}=\beta_R/\tau_s^{(d)}.
\]

The input state can represent zero temperature independently of a thermal
reference length. This release's kinetic adapter deliberately rejects it. A
future cold-reference calculation must dispatch before this adapter rather than
create a fictitious reference temperature. Extremely small positive beta still
inherits the conditioning limits of the legacy thermal calculation.

`d_ref` beta scans and species-1 temperature scans update both redundant
parameters. Other input values stay fixed. Scanning both independently in a
two-dimensional run is rejected. A legacy beta scan instead holds its legacy
wavevectors, thermal velocity, and temperature ratios fixed. To compare those
physical paths, convert **every point**, including all affected inputs.
Radial tables use the selected input convention; linear interpolation of
reciprocal temperatures in one convention generally differs from linear
interpolation in the other. Matching endpoints alone does not establish an
identical interpolated physical path.

`Kn` uses the selected length. With
`Kn_d = Kn_rho * sqrt(beta_R*alpha_R)`, the equivalent collision frequency is
`sqrt(muS*tauS_d)/(sqrt(2)*Kn_d)`. Numeric scan output remains in input units.
FPC output retains its documented derived thermal units. Reusing an output name
for a legacy run removes a stale `d_ref` metadata file.

## Verification

Run the available automated suite with:

```sh
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The suite compares raw complex determinants, every susceptibility element,
electric/magnetic fields, polarization, density/velocity fluctuations, heating
and all seven heating channels, and converged complex roots. Cases include
isotropic and anisotropic proton/electron plasmas, electron drift, drifting
alpha particles, oblique propagation, and reference beta 1, 0.3, and 0.02.
Root checks also require convergence before the secant iteration limit and a
scaled determinant below `1e-8`.

Integration checks cover parallel/perpendicular wavevector scans,
two-dimensional scans, frequency maps, beta scans converted point by point,
velocity and species-temperature scans, radial interpolation, both FPC
coordinate systems, weak collisions, example inputs, Python round trips, input
validation, and metadata when an output name is reused.

High-precision comparisons use `1e-10 + 5e-8*max(abs(a),abs(b))`. Rounded mode
files use `1e-4*max(1,abs(a),abs(b))`; FPC files use their printed precision.
The existing map fixtures contain singular grid points, so the map comparison
requires finite values in its nonsingular test grid.

The completed automated run passes **15 tests**. Across the high-precision
equivalence cases, the maximum `abs(a-b)/max(1,abs(a),abs(b))` is `4.668e-13`;
the maximum absolute difference in a frequency component is `6.939e-18`.
The legacy probe comparison is exact. Python compilation and `git diff --check`
also pass. FORD is not installed in the local environment, so its rendered
documentation build was not run locally; the report is included in the existing
documentation workflow.

Legacy regression builds commit
`c6c7616e1fc58eb6118e1f9423a3f0df6d52ef20`, the last pre-normalization revision,
and requires exact probe records for two existing input files. The reference
commit is pinned so merging this feature cannot turn the test into a comparison
against itself. Existing expected data files are not modified.

The full-resolution `example_map_par.in` was also run in isolation during the
initial implementation: its dispersion map and all four mode files matched the
pre-change build byte for byte. The first six `.roots` columns matched exactly.
Some remaining legacy `.roots` species columns differ because the old
non-map-scan path writes an uninitialized parameter array. The new `d_ref` path
initializes those columns. The tracked golden files additionally correspond to
older input ranges and do not reproduce the current example; they were retained.

## Modified files and future cold-fluid requirements

The initial implementation changed `src/vars.f90`, `src/functions.f90`,
`src/disprels.f90`, `src/fpc.f90`, `linfpclib/linfpc.py`, and
`contours/plot_plume_resonances.py`. It added the Fortran probe, Python test suite,
and matched example pair, and updated README, input/output documentation, FORD
source comments, and `.gitignore`. This follow-up pins the regression baseline,
extends path coverage, fixes the single-frequency heating array allocation,
aligns Python's beta/temperature tolerance with Fortran, and completes this
report.

The future cold-fluid stage must:

- Add cold susceptibility and dispatch each species to its kinetic or fluid response.
- Permit `tauS=0` and `betap=0` through a path independent of the thermal adapter,
  including a fully cold proton/electron plasma.
- Define cold-species eigenfunction, heating, collision, and FPC behavior;
  thermal distributions and temperature anisotropy are undefined for a cold species.
- Extend scan, diagnostic, and limiting-case tests across mixed kinetic/fluid
  systems and cold reference species without adding an artificial temperature.

No cold susceptibility has been added in this stage.
