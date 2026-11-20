# Stage 3B.2: cold-reference mixed kinetic/fluid dispersion

## Milestone A checkpoint — completed before production dispatch

The audit follows the normalization, Stage 0/1, Stage 2, Stage 3A and auxiliary
normalization reports. `disp0`'s species response is now held once in
`src/kinetic_species.inc`, included by the unchanged legacy context and the
new, separately compiled `kinetic_core` module. The shared harmonic cutoff is
in `kinetic_harmonics.inc`. No expression in the extracted block was changed.

`calc_chi_kinetic(species,omega,state,chi,ierr,bessel,zet_in)` takes a converted
warm species and a `kinetic_normalization` containing beta, velocity, kpar,
kperp and normalization anisotropy. Physical `species%alph_s` remains distinct.
Callbacks use PLUME's existing Bessel/Z functions without a module dependency
cycle. The core returns ordinary chi and never looks up physical species 1.
It has local work arrays and disables cache reuse. The legacy cache, dispatch,
root arithmetic and all existing numerical paths remain unchanged at this
checkpoint. IEEE exception state is saved/restored around the new calculation;
unrepresentable arithmetic returns status 5, including in trapping builds.

Production does not depend on Python extraction. The independent Stage 3B.1
prototype remains a test-only extraction; its loader now resolves the shared
Fortran includes, with its original transformations and numerical assertions
intact. The new core and prototype are compared in production and checked builds
at all four auxiliary scales, both signs of parallel wavevector, both frequency
half-planes, isotropic/anisotropic/drifting species and cold-reference states.
Physical-state bit snapshots and deliberately poisoned cold anisotropy are
retained. Tests also cover invalid geometry and numerical failure/recovery.

**Checkpoint results before Milestone B:** all 59 pre-existing tests passed,
including the exact pinned legacy comparison. Four new core-interface tests
passed. Maximum componentwise `abs(core-reference)/max(1,abs(reference))` across
prototype comparisons and scale invariance was `2.477e-9`, within the unchanged
Stage 3B.1 `atol=1e-10, rtol=5e-9` tolerance. No existing expected values or
numerical tolerances were changed. The first milestone's build compiled the
new objects separately without linking them into production; cold-reference
mixed input still failed at this checkpoint.

## Milestone B: native cold-reference dispatch

Input reading and `disp` recognize a cold physical reference using exact
`spec(1)%tau_s==0`. `validate_cold_reference` requires collisionless `d_ref`,
`betap=0`, positive Alfvén speed, valid species parameters and supported
root/determinant options without scans. It allows warm additional species,
but rejects exact `kpar=0` when any are warm. The original thermal adapter is
never called for a cold reference. Warm-reference all-kinetic/mixed cases keep
their original adapter and evaluation paths; fully cold evaluation retains
its original cold responses and native wave assembly.

The shared cold preflight still evaluates each cold response in native inertial
units using `calc_chi_cold`. Its equations and statuses are unchanged. For an
oblique cold-reference mixture, `cold_reference::calc_warm_auxiliary` finds the
first positive-temperature species and chooses its native temperature as
`b_star`. It uses the validated auxiliary adapter to prepare scales and converts
one warm species at a time. Cold components never undergo temperature conversion.
Each warm response is returned by the explicit kinetic core and stored in the
original species slot. No species ordering, mass/charge/density ratio, physical
anisotropy or drift changes. No `solver_spec` allocation is needed in this path.

The existing native assembler (`disp_all_cold`, whose historical name is
retained) also assembles the new mixture after all responses have been stored:

\[
 \epsilon=I+\sum_s\chi_s,\qquad
 \Lambda=\epsilon+nn^T-n^2I,\qquad
 n=(\kappa_\perp,0,\kappa_\parallel)/(aw).
\]

It adds identity and electromagnetic terms once and uses the existing determinant
expression. Its input is ordinary susceptibility, with neither extra frequency
scaling nor a transpose/conjugation conversion.

For exact `kperp=0`, `disp_mixed_parallel` already evaluates warm responses using
physical `sqrt(mu_s*tau_s_d)` and skips cold entries. The new path reuses that
routine without changing its equations. It needs no auxiliary thermal scale.
A separate IEEE guard protects this new caller and returns numerical status 5
if the native thermal-speed products are unrepresentable. The old warm-reference
caller is unchanged. Exact perpendicular warm propagation remains unsupported;
no small wavevector is substituted.

The existing `rtsec`, `refine_guess` and `map_search` methods are unchanged.
Cold poles retain statuses 2–4; invalid input is 1 and numerical failure is 5.
Auxiliary range errors map to status 5. A failed warm response invalidates all
stored tensors, including any species computed earlier in that evaluation.
The new core and parallel guard save and restore the caller's IEEE halting
modes and exception flags, so checked builds return a failure instead of aborting
on new-path arithmetic overflow. No floors, damping or trial-frequency offsets
are added. Map search retains NaNs for invalid samples and skips them during
candidate selection. Input charge/current balance warnings remain warning-only,
with no automatic density, drift or reference-frame adjustment.

The direct evaluator prints raw `det(Lambda)` for either cold-reference case;
thermal output scaling remains unchanged for warm references. The normalization
sidecar identifies the first warm species and its provisional `b_star`. The
same policy is recorded for exact-parallel input, although that analytical
limit uses physical species speeds directly. Supplementary eigenfunctions,
heating/FPC and cold-containing collisions and scans remain unavailable.

## Normalization mapping and numerical conditioning

For `a=v_A,ref/c`, `theta_s=tauS_s_d`, and the first warm species `j`,

\[
 b_*=\theta_j,\quad B=b_*,\quad V=a\sqrt{b_*},\quad
 K_i=\kappa_i\sqrt{b_*},\quad t_s=b_*/\theta_s,\quad\alpha_*=1.
\]

`alpha_*` is independent of every physical `alpha_s`. The physical reference
stays cold (`betap=tauS_1=0`); `B` is an internal coordinate scale, not a
replacement physical beta. Cancellation gives

\[
 A_s=\frac{B\mu_sD_s}{V^2Q_s^2}=\frac{\mu_sD_s}{a^2Q_s^2},\quad
 g_s=\mu_s/Q_s,\quad \delta_s=w-K_\parallel u_s/\sqrt B
 =w-\kappa_\parallel u_s,
\]
\[
 \frac{v_{t\parallel,s}}{v_A}=\sqrt{\mu_s\theta_s},\quad
 \zeta_{ns}=\frac{\delta_s-ng_s}{\kappa_\parallel\sqrt{\mu_s\theta_s}},\quad
 \lambda_s=\frac{\kappa_\perp^2Q_s^2\alpha_s\theta_s}{2\mu_s}.
\]

All nine kinetic components remain those derived in the Stage 3B.1 report.
The unchanged legacy normalization still uses physical `alpha_R` and
`k_rho=k_d*sqrt(beta_R*alpha_R)`; the new normalization does not silently replace
that physical anisotropy. There is no fictitious reference population or thermal
reference speed. The new core can receive an explicit normalization state and
is independently tested at `b_*=0.01, 0.1, 1, 10`, while production chooses the
first warm species automatically.

The new core is compiled separately with IEEE-preserving flags. The legacy
compilation flags and arithmetic remain intact. Prototype/core discrepancies
are consistent with the previously measured cancellation in the unchanged
anisotropic electron harmonic sums; a mathematically arbitrary scale is not
bitwise arbitrary in floating-point arithmetic. Selecting the first warm
species is provisional and does not guarantee representable intermediate
products for arbitrarily extreme temperature spreads. Such detected failures
return status 5 instead of changing the physical state. No optimization of
scale selection or general kinetic conditioning was introduced.

The legacy Bessel cache keeps its original keys and reuse behavior. New core
calls use local arrays and recompute Bessel functions, so changing the first
warm species or auxiliary scale cannot reuse stale entries. Harmonic cutoff 128,
signed-wavevector continuation and the existing special functions are shared.
The production Makefile compiles Fortran modules/includes directly; Python
extraction is used only by the independent historical prototype tests.

## Independent validation

Ten new Milestone B tests compare each species susceptibility, all nine elements
of `epsilon` and `Lambda`, and the determinant against the existing independent
NHDS/Stix SciPy harmonic sum and cold-fluid linear solve. Configurations include:

- Cold protons with warm electrons; drifting warm electrons of both drift signs;
  anisotropic warm electrons.
- Cold protons with warm electrons/warm alphas; cold electrons/warm drifting
  alphas; warm electrons/cold drifting alphas.
- Four-component plasmas with different positive temperatures, masses and
  signed charges, including nonunit densities.
- Nonneutral/current-carrying equilibria, including a moving reference species;
  warnings appear once and input parameters are not adjusted.
- Exact parallel and oblique wavevectors, both signs of `kpar`, complex frequencies
  in both half-planes, and `v_A/c=0.001, 0.05`.
- Native evaluation with poisoned thermal solver scalars and unallocated
  `solver_spec`; varied cold anisotropy placeholders give bitwise equal results.
- Poles followed by valid evaluations, and checked numerical failures in both
  oblique auxiliary arithmetic and native parallel thermal-speed products.
- Supplied roots, a resolved parallel frequency map, the example and Python
  wrapper, and explicit rejection of unsupported diagnostic/scan/collision or
  exact-perpendicular requests.

Core/prototype tests are independent of production dispatch. The prototype
continues to generate its separate test module and uses its original numerical
comparisons. The production core receives its normalization through Fortran,
not through that test extraction.

Root references are calculated independently with SciPy using the physical
NHDS/Stix/cold tensor. Broad guesses sometimes find different legitimate
branches with SciPy and PLUME. Tests therefore identify a branch with the
independent solver, then initialize PLUME at `0.98*w_reference + 1e-4*i`.
This is an unsolved nearby guess, not the final root; convergence, real/imaginary
frequency agreement, and both independent/production residuals are checked.
No root method or stopping condition is changed to force a match.

Reordering the non-reference warm electron/alpha species changes the chosen
auxiliary beta from 0.2 to 0.03. Tests remap species susceptibilities and compare
the assembled tensor and the same independently identified root, in both builds
and at exact-parallel and oblique propagation. This directly tests the numerical
selection policy while keeping the physical plasma fixed.

### Tolerances and measured differences

| Check | Tolerance | Maximum measured discrepancy |
| --- | --- | --- |
| Milestone A core/prototype and scale invariance | componentwise `atol=1e-10, rtol=5e-9` | `2.477e-9` scaled element error |
| Independent mixed species/epsilon/wave tensor | componentwise `atol=1e-5, rtol=2e-7` | `6.223e-9` scaled element error |
| Species-order invariance | componentwise `atol=1e-5, rtol=5e-9` | `1.891e-9` scaled element error |
| Root real and imaginary components | `atol=1e-10, rtol=1e-8` | `2.258e-12` absolute complex-frequency difference |
| Order-invariant roots | absolute difference `<5e-10` | all comparisons pass |
| Root scaled determinant | `<1e-8`, PLUME's existing convergence criterion | `3.239e-10` production residual |
| Rounded map roots | absolute difference `<5e-6` | resolved parallel branch passes |

Scaled element errors are `abs(actual-reference)/max(1,abs(reference))`.
The determinant residual is `abs(det(Lambda))/max(abs(Lambda))**3`.
The largest new residual occurs for a drifting parallel branch. Its frequency
still agrees tightly with the independent solution; the residual tests use
PLUME's existing `1e-8` criterion, not a stricter new stopping rule. Previous
test tolerances and numerical expectations are unchanged.

Representative dimensionless roots (finite electron mass `mu_e=1836`):

| Physical configuration | Wavevector / parameters | Root `omega/Omega_R` |
| --- | --- | --- |
| Cold protons, stationary warm electrons | `kpar=0.06`, `kperp=0`, `a=0.001`, `theta_e=0.2`, `alpha_e=0.8` | `0.05764490456558`, negligible imaginary part |
| Same, oblique damped mode | `kperp=0.05` | `0.05917021837673 - 1.5890035767e-6 i` |
| Drifting warm electrons | `u_e=-0.15`, `kpar=0.06`, `kperp=0.05` | `0.01903007290150 - 4.0714240679e-4 i` |
| Cold protons, warm electrons/warm alphas | `theta_alpha=0.03`, `alpha_alpha=1.2`, `D_alpha=0.05`, `u_alpha=0.3`, same wavevector | `0.04318326542752 - 1.9419515790e-5 i` |
| Cold protons/electrons, growing warm alpha beam | `theta_alpha=0.03`, `alpha_alpha=1`, `u_alpha=2`, `D_alpha=0.05`, `D_e=1.1`, `u_e=0.2/1.1`, `kpar=0.3`, `kperp=0.1`, `a=0.01` | `0.30574584309028 + 0.07964933433436 i` |

## Full regression and scope

With GNU Fortran 11.4.0, the requested build and test commands succeed:

```sh
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

**73 tests pass:** the 59 pre-existing tests, four Milestone A tests and ten
Milestone B tests. The pinned legacy comparison remains **exact**. Checked
builds use bounds checking and floating-point traps. Existing warm-reference
mixed and fully cold results pass unchanged, and the Stage 3B.1 prototype
numerical tests remain intact. Python compilation and `git diff --check` pass.
Intel Fortran and rendered FORD documentation were not tested locally because
those tools are unavailable.

Obsolete input-rejection assertions were updated only after the Milestone A
checkpoint: previously rejected cold-reference/warm cases now either test
successful production use or the still-unsupported warm perpendicular limit.
The old prototype linkage assertion now requires the production adapter/core
while continuing to forbid the test-extracted module in `plume.e`. Test build
commands include the new module dependencies, and the prototype loader expands
the shared Fortran includes. No expected numerical values or existing numerical
tolerances were altered.

The example `inputs/example/mixed_cold_protons_warm_electrons.in` runs from the
repository root and prints `0.059170 - 1.5890e-6 i`. It uses `d_ref`, `betap=0`,
a cold proton reference, warm electrons and a supplied complex-frequency guess,
without scans or supplementary diagnostics.

## Modified files and remaining work

| File | Purpose |
| --- | --- |
| `src/kinetic_species.inc` | Single shared source of unchanged collisionless species expressions. |
| `src/kinetic_harmonics.inc` | Shared original harmonic cutoff. |
| `src/kinetic_core.f90` | Explicit normalization/species interface, local work arrays, numerical-failure guards. |
| `src/cold_reference.f90` | First-warm-species scale selection, one-species conversion/evaluation, ordered warm responses. |
| `src/disprels.f90` | Include shared expressions; new cold-reference dispatch, reuse native/parallel assembly, raw direct output. |
| `src/vars.f90`, `src/functions.f90` | Native cold-reference classification/validation and metadata; bypass thermal adapter. |
| `src/auxiliary_normalization.f90` | Header updated for production integration; conversion implementation unchanged. |
| `Makefile` | Module and include dependencies, IEEE-preserving core/adapter flags, include files in source packaging. |
| `linfpclib/linfpc.py` | Accept supported cold-reference warm mixtures, retain restrictions and zero-beta consistency. |
| `tests/kinetic_core_probe.f90`, `tests/test_kinetic_core.py` | Milestone A comparisons and numerical failure/state tests. |
| `tests/test_cold_reference_mixed.py` | Independent production assembly/roots, ordering, geometry, failure, map, wrapper and example tests. |
| `tests/auxiliary_kernel.py` | Resolve includes for the historical prototype extraction. |
| `tests/test_auxiliary_normalization.py`, `tests/test_mixed_dispersion.py`, `tests/test_fully_cold_dispersion.py`, `tests/test_normalization.py` | Link/build dependencies and obsolete support-boundary expectations only. |
| `inputs/example/mixed_cold_protons_warm_electrons.in`, `.gitignore` | Working example and tracking exception. |
| `README.md`, `input.md`, `output.md` | Current capability, normalization policy, output/status and limitations. |
| `normalization.md`, `tests/AUXILIARY_THERMAL_NORMALIZATION.md`, `tests/FULLY_COLD_DISPERSION_REPORT.md` | Historical reports point to current support. |
| `tests/COLD_REFERENCE_MIXED_REPORT.md`, `ford_project.md`, `.github/workflows/doc.yml` | This report and its documentation publication. |

The cold susceptibility equations, kinetic physics, root-finding method,
physical reference definitions and species ordering are unchanged. Unrelated
workspace modifications were preserved. Further work remains separately scoped:
warm exact-perpendicular limits; cold-fluid eigenfunctions, heating and FPC;
collisions; cold-containing/model-transition scans; optimized scale selection
and wider conditioning guarantees. None is implemented by Stage 3B.2.
