# Stage 2: mixed kinetic/cold dispersion relations

## Audit and implementation

Input reading (`functions::read_in_params`, `spec_read`) leads through
`vars::sync_solver_state` to `disprels::disp`. Supplied guesses use
`refine_guess -> rtsec -> disp`; map samples and refined map candidates use the
same evaluator. Previously the adapter rejected zero temperatures, and `disp0`
would otherwise divide by them in thermal/Bessel/resonance arguments.
`disp1` has the same thermal assumptions. `test_disp` also called eigenfunction
and heating calculations unconditionally. FPC builds thermal distributions.

The existing tensor conventions are compatible: both the validated cold kernel
and the kinetic routines return ordinary dimensionless `chi_s`. Thus the species
sum needs neither a sign/transpose change nor extra frequency scaling. The
reference remains warm, so the established perpendicular-gyroradius conversion
`k_rho=k_d*sqrt(beta_R*alpha_R)` remains valid for every warm component.

Classification is exactly `tau_s==0`, derived from the authoritative input
state; there is no extra species-model flag or small-temperature threshold.
Negative temperatures are rejected. Cold components require collisionless
`d_ref`; species 1 has `tau_s=betap>0`. The adapter copies cold parameters without
converting their temperature. It converts only positive-temperature components.
Cold anisotropy is ignored; no artificial thermal speed is supplied.

The small `disp` dispatcher retains the legacy all-kinetic path and delegates
mixed evaluations to `disp_mixed`. The latter evaluates cold responses using
native `k*d_R`, `v_A/c`, frequency and signed species inputs before any thermal
expression. In `disp0`, cold components bypass the thermal calculations and
contribute their stored `susc` entries directly to the six independent tensor
elements. Every warm kinetic expression and the final electromagnetic/vacuum
assembly remain unchanged. Cold entries also bypass thermal cache checks.

Keeping the mixed checks in a separate routine is numerically significant under
the existing GNU fast-math build: enlarging `disp` changed compiler optimization
of unrelated diagnostic arithmetic at roundoff level. The final separation
preserves the pinned legacy regression **exactly**, without changing expected
values or loosening the regression assertion.

## Exact parallel limit

The legacy general kinetic formula contains divisions by the perpendicular
Bessel argument at `kperp=0`. Mixed runs use its analytical parallel limit in
`disp_mixed_parallel`, with no wavevector floor. Only the `n=+/-1` harmonics
contribute to transverse response and `n=0` contributes longitudinally. For a
warm species in native inertial units, define

\[
v_t/v_A=\sqrt{\mu_s\tau_s},\quad
\zeta_n=(w-k_\parallel u_s-ng_s)/(k_\parallel v_t),\quad
F_n=\alpha_s-1+
\frac{\alpha_s(w-k_\parallel u_s)+(1-\alpha_s)ng_s}{k_\parallel v_t}Z(\zeta_n).
\]

Then

\[
\chi_{xx}=\chi_{yy}=\frac{A_s}{2w^2}(F_{-1}+F_1),\quad
\chi_{xy}=\frac{iA_s}{2w^2}(F_1-F_{-1}),\quad
\chi_{zz}=\frac{2A_s}{(k_\parallel v_t)^2}[1+\zeta_0Z(\zeta_0)].
\]

The other drift-coupling components vanish at `kperp=0`. This uses PLUME's
existing causal `zet_in`, including negative parallel wavevectors. It is only
selected for mixed runs; the all-kinetic path is untouched. Exact perpendicular
propagation (`kpar=0`) still needs the warm perpendicular limit and is explicitly
unsupported for mixed runs. Parallel equilibrium drifts and arbitrary oblique
directions with nonzero `kpar` are supported.

## Failure handling and diagnostic boundaries

`vars::dispersion_status` is reset on every evaluation: 0 means valid, 1 means
invalid cold input/frequency, 2–4 are the cold kernel's physical poles, and 5
means nonfinite numerical arithmetic. Failed evaluations invalidate `susc`,
`lam` and resonance arrays with NaNs, then return a NaN determinant. The next
valid evaluation recomputes the response. IEEE helper checks live in the
kernel's existing precision-preserving compilation unit so production
`-ffast-math` cannot remove them. No damping, collisions or numerical floors
are introduced.

Mixed `rtsec` retains its original secant update and maximum iteration count.
It rejects nonfinite evaluations with `iflag=-1`, requires a scaled determinant
residual below `1e-8` at convergence, and reports exhaustion with `iflag=-2`.
Supplied-guess failures terminate with an informative diagnostic. All-kinetic
failure behavior is unchanged. Map searches mask isolated invalid points only
for candidate selection; the actual map determinant and magnitude remain NaN.
Failed candidates are omitted, valid candidates continue, and one summary
warning reports invalid samples. Mask values never replace a physical tensor.

Mixed runs support `option=-1`, `0`, or `1`, with `nscan=0`. The direct evaluator
prints its status and explicitly omits supplementary diagnostics. Eigenfunction
and FPC entry points are guarded before thermal evaluation or velocity-grid
allocation. `susc_low` for cold components is NaN, never a fabricated resonant
heating result. Parameter/radial scans and cold-containing collisions are
rejected. Temperature scans that reach or leave zero are rejected even when
the initial state is all kinetic. Comprehensive mixed scans are deferred.

Charge and current balance retain warning-only semantics. Input reading reports
`sum(D_s/Q_s)` and `sum(D_s*vv_s/Q_s)` once, with neither automatic density/drift
adjustments nor frame changes. Nonstandard homogeneous equilibria can require
external forces or currents to sustain the background.

## Independent verification

`mixed_probe.f90` directly calls the production evaluator and existing root
finder, emitting every species tensor, the assembled wave matrix, determinant,
status, root and convergence flag. Tests build production and checked variants;
the latter uses bounds checks and traps for invalid/divide-by-zero/overflow.

The test-only Python reference adapts the NHDS/Stix bi-Maxwellian harmonic-sum
formulas using independent SciPy `ive` and `wofz` implementations. It retains
finite electron mass and uses the Stage 1 independent cold-fluid linear solve
for cold components. The NHDS BSD notice and attribution are retained in
`NHDS_LICENSE.txt`. It neither calls PLUME nor copies its expanded kinetic
expressions. An independent SciPy root solve uses this tensor; checks at 60 and
90 harmonics verify truncation adequacy for the selected roots. The tests compare
the same branches, since broad guesses can converge to different legitimate
branches with different root solvers.

Tests include warm protons/cold electrons; warm protons/electrons with drifting
cold alphas; both drift and parallel-wavevector signs; multiple cold species;
nonunit masses/densities; anisotropic warm references; growing/damped complex
evaluation frequencies; nonneutral/current-carrying states; exact parallel and
oblique propagation; singular samples followed by successful evaluations;
unsupported configurations; the example; and Python input construction.
Warm species susceptibilities are also required to remain bitwise identical
when a different component is changed from warm to cold.

Componentwise independent tensor comparisons use `rtol=2e-7, atol=1e-5` for the
large, ordinary susceptibilities at `v_A/c=1e-3`. Root comparisons separately
check real and imaginary parts with `rtol=1e-8, atol=1e-10`, convergence before
75 iterations, and normalized residual below `1e-10`. These tolerate the
unchanged kinetic special-function approximations while resolving every
species contribution. The largest measured tensor error in the oblique matrix
is `9.640e-8` using `abs(actual-reference)/max(1,abs(reference))`; the independent
root comparisons are much closer (results below).

The alpha-species kinetic-to-cold sequence holds all physical parameters fixed
and uses `tau_d=0.1, 0.01, 0.001, 0.0001, 0.00001`. Errors are
`max(abs(chi_kinetic-chi_cold))/max(abs(chi_cold))`:

| Alpha drift | kperp*d_R | First temperature error | Last temperature error |
| --- | --- | --- | --- |
| 0 | 0 | `8.940e-3` | `8.822e-7` |
| 0 | 0.001 | `8.940e-3` | `8.826e-7` |
| 0 | 0.05 | `8.825e-3` | `8.711e-7` |
| 0.25 | 0 | `1.151e-2` | `1.132e-6` |
| 0.25 | 0.001 | `1.151e-2` | `1.132e-6` |
| 0.25 | 0.05 | `1.140e-2` | `1.122e-6` |

Each sequence converges monotonically by approximately a factor of ten per
temperature step. No machine-precision cold limit is asserted: very small warm
temperatures inherit cancellation and special-function conditioning limits.

## Final verification

`make -j2` succeeds. The full requested command

```sh
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

passes **38 tests**: all 27 previous tests and 11 mixed-dispersion tests. The
pinned pre-normalization regression matches exactly. The largest independent
root difference is `2.113e-11` in `omega/Omega_R`, across warm-proton/cold-electron
parallel and oblique modes, drifting cold alphas with warm protons/electrons,
and an oblique cold-electron drift case. Both real and imaginary frequency
comparisons and all residual/convergence assertions pass.

The checked build passes the mixed tensor and pole-recovery tests with bounds
checks and floating-point traps enabled. The production example and Python
wrapper inputs also run successfully. Python compilation and `git diff --check`
pass. Tests used GNU Fortran 11.4.0; Intel Fortran and rendered FORD documentation
were not tested locally because those tools are unavailable.

## Changed files and remaining work

Production changes: `src/vars.f90` (classification, adapter, restrictions/status),
`src/functions.f90` (input/scan validation and warning labels),
`src/disprels.f90` (assembly, parallel limit, mixed-only failure/diagnostic guards),
`src/cold_susceptibility.f90` (IEEE helpers; cold equations unchanged),
`src/fpc.f90` (entry guards only), and `Makefile` (module dependency).
`linfpclib/linfpc.py` accepts supported cold inputs and rejects unsupported
combinations. The existing normalization test's link command includes the new
dependency; its assertions and pinned expected data are unchanged.

New tests: `mixed_probe.f90`, `mixed_reference.py`, `test_mixed_dispersion.py`,
and `NHDS_LICENSE.txt`. The example is
`inputs/example/mixed_warm_protons_cold_electrons.in`, with an ignore exception.
README, input/output/install documentation and historical normalization/cold
reports point to current functionality. This report records the Stage 2 design
and results. `ford_project.md` and `.github/workflows/doc.yml` publish this report
with the static documentation. Unrelated working-tree changes were preserved.

Future stages must support a cold reference before the thermal adapter, add
cold eigenfunctions/heating/FPC and collision models as separately specified,
and define mixed/model-transition scans. The warm perpendicular limit is also
still needed. No Stage 3 implementation is included here.
