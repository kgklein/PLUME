# Stage 4A: mixed cold/kinetic wavevector scans

## Audit and implementation plan (before source changes)

`plume` reads plasma inputs, then options 1/2 initialize roots through
`map_search` or `refine_guess` before reading `scan_input_N`. Option 1 calls
`om_scan` for serial scans; option 2 calls `om_double_scan`; option 4 calls
`map_scan -> map_search` independently at each wavevector. Option 5 interpolates
physical plasma parameters as well as wavevectors and remains outside this stage.

`vars::sync_solver_state` and `validate_cold_reference` reject every cold-containing
run with `nscan>0` and all options except -1/0/1. These are implementation-stage
restrictions, not susceptibility limitations for nonzero parallel wavevectors.
The Python wrapper repeats the blanket rejection. Species classification is
already exact tauS=0, including arbitrary cold/warm species ordering.

`get_out_name`, `get_double_out_name` and `set_map_pointers` configure pointer-based
scan coordinates, names, spacings and established numeric formats. Global style 0
supports kperp/type 0 and kpar/type 1. Style -1 supports vector endpoints/type 0,
angle/type 1 and magnitude at fixed angle/type 2. Other style/type combinations
change the plasma and are not enabled for cold-containing runs here.

Both continuation routines call `sync_solver_state` directly at scan points;
that is invalid for a cold reference. Their secant calls already use the correct
production `disp`, but they store returned roots without checking failure flags.
The two-dimensional routine saves row-start roots and restores them after each
inner traversal. Existing scans can request heating/eigenfunctions, which remain
unavailable for cold-containing configurations. Existing susceptibility tensor
output needs no new physics and may be retained.

The compound geometry code contains shared mutable coordinate state in map scans
and an uninitialized theta_q in an inner fixed-angle two-dimensional path.
A small cold-scan coordinate adapter will avoid those assumptions while retaining
the existing continuation loops and IO helpers. It will support Cartesian kperp/
kpar grids and the existing ordered (theta, magnitude) polar grid; overlapping
coordinate definitions will be rejected. Legacy all-kinetic geometry is preserved.

`map_scan` also hardwires search bounds using reference anisotropy in d_ref.
That is undefined for cold species. Cold-containing map scans will retain the
user's supplied &maps bounds. All-kinetic span behavior will remain unchanged.

The dispersion evaluator already dispatches warm-reference mixed, fully cold and
cold-reference mixed correctly. Auxiliary b_star is the first warm species'
physical tauS, independent of wavevector. Bessel cache keys use physical lambda;
kpar changes reuse coefficients and kperp changes invalidate/resize them. Neither
the evaluator, auxiliary adapter, Z backend, nor harmonic/cache policy needs changes.

Validation will precede root initialization, reject kpar=0 (including explicit
90-degree angle endpoints), and prohibit cold diagnostics/physical-parameter
scans. Failures will produce NaN output rows and warnings, retaining the last
valid continuation seed. No new root method, damping or coordinate floors will
be added. The original numeric column order and precision will be retained.

## Completed implementation

A new, separately compiled `scan_support` module handles preflight geometry,
initial-point setup, and numerical acceptance/failure isolation. It contains
**no susceptibility equations, root-finding iteration or branch tracker**.
The existing `om_scan`, `om_double_scan`, `map_scan`, `disp`, and `rtsec` remain
the production scan/dispersion/root framework.

Input validation now permits cold-containing `option=1/2/4` wavevector scans,
while retaining the collisionless d_ref requirement and all diagnostic and
plasma-parameter/radial restrictions. `read_scan_input` validates the complete
scan before root initialization. Option 1/2 read scans before initializing
supplied/map guesses. Newly enabled compositions initialize at the actual
first point; historical all-kinetic initialization is retained. Optional
`heating`, `eigen` and `tensor` flags have explicit false input defaults.

Supported scan coordinates and existing input meanings:

| Style/type | Coordinate | Stage 4A workflows |
| --- | --- | --- |
| 0/0 | kperp, swi -> swf | Serial, Cartesian grid, independent maps |
| 0/1 | kpar, swi -> swf | Serial, Cartesian grid, independent maps |
| -1/0 | Input vector -> (swi,swf) | Serial, independent maps |
| -1/1 | Input angle -> swi degrees, fixed magnitude | Serial, first axis of polar grid, independent maps |
| -1/2 | Input signed magnitude -> swf, fixed angle | Serial, second axis of polar grid, independent maps |

Cartesian grids require distinct kperp/kpar axes; polar grids retain ordered
(theta,magnitude) axes. Ambiguous overlapping axes are rejected. Linear and
logarithmic spacing retain their meaning; logs require positive scanned
coordinates. Negative kpar is supported with linear spacing. The historical
signed magnitude is `kpar/cos(atan(kperp/kpar))`, so negative-propagation
magnitude scans use negative endpoints. Immutable native coordinate frames
avoid shared-state aliasing in newly enabled compound/grid/map traversals.
All-kinetic geometry calculations are kept on their existing pathway.

Preflight checks every wavevector sample (and both independent grid axes),
including hidden nres substeps. kpar=0 fails with a Stage 4A diagnostic.
Degree-space checks catch exact 90-degree endpoints despite cos(pi/2) rounding.
Checking a linear numerator before division also detects, without an epsilon,
the zero sample in a -.1 -> .2 scan with three steps, which ordinary arithmetic
could misrepresent as ~1e-17. No coordinate is floored or replaced.
Fully cold single-wavevector perpendicular functionality is unchanged.

## Normalization and physical-state integrity

Outputs remain in the input convention. For warm-reference d_ref,
`k_rho=k_d*sqrt(betap*alphS_ref)` remains inside the unchanged adapter. Wavevector
scan points use that adapter only when the physical reference is warm.
Fully cold and cold-reference mixed scans bypass it: no reference thermal
speed, gyroradius or solver_spec allocation is constructed.

Cold-reference oblique responses still use the first warm species' physical
`b_star=tauS(j)`. Because all species parameters and v_A/c remain fixed in a
wavevector scan, b_star is constant; only the derived K components change.
The exact-parallel path continues to use physical species thermal velocities.
Mass/charge/density, temperatures, anisotropies, drifts, beta and v_A/c retain
their established definitions and values. Test probes compare their binary
storage before/after scans and reject any cold-reference solver_spec allocation.
Cold anisotropy placeholders, including negative values, do not affect roots.
Charge/current warnings remain once-per-input, with no balancing constraints.

Cold-containing `map_scan` uses the supplied &maps frequency window at each
point; it does not apply the legacy reference-anisotropy span heuristic.
All-kinetic map-span behavior is unchanged. Existing normalization/Bessel/Z
metadata remains in use; no auxiliary-scale input is introduced.

## Continuation and numerical failures

The secant algorithm, bracketing multipliers, tolerance, maximum iteration
count and species dispatch are unchanged. At each point a thin guard calls the
existing solver, evaluates its returned frequency, and checks finite tensor/
determinant, status and `abs(det(Lambda))/max(abs(Lambda))**3 < 1e-8`.
Successive divisions avoid overflowing the normalization cube. The wrapper
saves/restores IEEE flags and traps, including in legacy mode, so failed scan
arithmetic does not abort a checked run or poison the next evaluation.
The legacy all-kinetic flag 75 cannot distinguish exhaustion from last-iteration
success; it is conservatively rejected in scans. Accurate/cold failure flags
retain their existing meanings. This is acceptance checking, not a new solver.

Only accepted roots replace each branch's continuation seed. Failed frequencies
and susceptibility slots are NaN, warnings identify coordinates and failure
flags, and later points retry from the last valid seed. No finite stand-ins or
artificial damping are introduced. Two-dimensional scans retain row-start root
save/restore. Independent map searches retain their existing failure masking
and do not associate branches between maps.

A real coarse-scan benchmark demonstrates recovery: warm protons/cold electrons,
kpar=.01, kperp=.05 -> 20 in four output steps with nres=2. Trials at kperp
2.54375 and 5.0375 fail; the latter writes a NaN row. Later points return valid
roots, including ~4.151441 at kperp=10.025. Their matrices and roots agree with
independent evaluation, and checked builds survive. **Recovery selects a
higher-frequency branch** after this deliberately large gap. This is an
identified continuation limitation, not proof that the original branch was
retained. Reduce steps near crossings, poles, strong damping or failed points.
No new branch-tracking algorithm is supplied.

## Output and diagnostic boundaries

The main root-only `.modeN` layout is unchanged: kperp, kpar, betap, vtp,
omega_r, omega_i, followed by six input parameters per species. Species order,
column meanings, filenames, scientific precision and grid blank separators
are retained. Existing optional susceptibility tensor files retain 4+18*nspec
real columns. They require no cold perturbation physics.

Failed output rows have NaN frequency/tensor slots and intact physical parameter
columns. Failures at hidden nres substeps warn even when no row is scheduled.
Existing all-kinetic eigenfunction/heating requests remain available; failures
skip their calculation and mark their slots invalid. Cold-containing requests
for eigenfunctions or heating are rejected before root initialization. FPC,
cold collisions and radial/plasma-parameter scans remain unavailable. Stage 4B's
future kinetic-species heating/FPC policy is **not implemented** here; cold
quantities remain unavailable, not physical zeros.

## Independent validation

`test_mixed_scans.py` runs normal and checked production executables and the
actual continuation loops. Checked flags are `-O0 -g -fdefault-real-8
-fcheck=all -ffpe-trap=invalid,zero,overflow`; production flags retain GNU's
existing optimization/fast-math settings. The new validation module alone is
compiled with precise IEEE semantics. No kinetic/cold/Z/Bessel routines were
changed. Test-only probes assert physical-state invariance and expose final
roots at full precision without changing production formats.

The reference is the independent NHDS/Stix bi-Maxwellian/cold implementation
with SciPy wofz and a converged 512-harmonic reference. Each saved scan row is
compared with a new standalone PLUME root calculation in a fresh process at the
same physical point, plus independent SciPy roots. Both frequency components
and scaled matrix determinants are checked; a small determinant alone is not
used as proof of root/branch agreement.

The matrix covers all-kinetic, fully cold proton/electron, warm protons/cold
electrons, cold protons/warm electrons, and a four-species drifting system
with two warm and two cold populations. It includes both signs of kpar,
serial forward/reverse kperp traversals, Cartesian and polar grids, vector/
angle/magnitude scans, log spacing, exact mixed/cold parallel propagation,
supplied and map-initialized roots, independent map scans, both normalizations
where allowed, both Z backends and all three Bessel spellings. Perpendicular
samples, diagnostic requests, and cold plasma-parameter scans fail explicitly.

**16 new tests pass**, comparing **291 mode-file rows**.

| Comparison (absolute frequency-component differences) | Maximum |
| --- | ---: |
| Printed omega_r vs standalone root | 4.262464e-07 |
| Printed omega_i vs standalone root | 3.752771e-10 |
| Printed omega_r vs independent root | 4.262461e-07 |
| Printed omega_i vs independent root | 3.754970e-10 |
| Full-precision final omega_r vs independent root | 7.290279e-13 |
| Full-precision final omega_i vs independent root | 5.479777e-14 |
| Scaled independent determinant at rounded saved frequency | 4.026568e-11 |
| Scaled determinant at standalone refined frequency | 1.041095e-17 |

The mode file has seven significant digits. Thus the ~4.27e-7 maximum real
frequency difference is output rounding (the coarse recovery/large-argument
cases have order-unity frequencies), not a 4e-7 solver discrepancy. Mode-row
comparisons allow `1e-11 + 6e-7*abs(reference)` against standalone values, and
`3e-10 + 1e-6*abs(reference)` against independent roots. Full-precision final
root comparisons use 2e-9 absolute acceptance and achieve the errors above.
Map `.roots` files retain their lower existing precision; comparisons allow
6e-5 relative complex-frequency error. All normalized residual checks require
less than 1e-8. Forward/reverse traversal agrees within retained output precision
on the selected non-crossing branches.

Representative final roots at kpar=.09, kperp=.05, v_A/c=.01:

| Plasma | PLUME root (full precision) | Independent root | Real error | Imaginary error |
| --- | --- | --- | ---: | ---: |
| Kinetic p/e | 0.0907182667350623-8.8663986e-05i | 0.0907182667357914-8.86639859e-05i | 7.290e-13 | 5.479e-14 |
| Cold p/e | 0.0885984956050465-3.12694318e-63i | 0.0885984956050908-1.36820955e-70i | 4.431e-14 | 3.127e-63 |
| Warm p, cold e | 0.0914398151231068-5.15990039e-05i | 0.0914398151231447-5.15990039e-05i | 3.790e-14 | 7.942e-16 |
| Cold p, warm e | 0.0878844079634319-5.47392357e-06i | 0.0878844079638679-5.47392357e-06i | 4.360e-13 | 2.708e-15 |
| Four-species drifting mixture | 0.137773295801836-6.81522007e-05i | 0.137773295801836-6.81522007e-05i | 0.000e+00 | 0.000e+00 |

The complete reproducible measurements, including both build modes and signs,
are in `MIXED_SCAN_RESULTS.json`. No historical expected numerical results or
existing tolerances were changed.

## Cache checks

The harmonic selector, Bessel cache, kinetic/cold susceptibility includes,
auxiliary adapter and Faddeeva sources are byte-for-byte unchanged. Cache keys
remain physical lambda, independent of complex frequency, kpar and auxiliary
normalization. Cold species acquire no valid Bessel entries.

A warm-proton/cold-electron scan kperp=10 -> 50 -> 10 at kpar=.01 increases
lambda_p from 50 to 1250, selecting up to N=263, then reduces the active cutoff
while retaining capacity. Roots at every saved point agree with standalone and
independent calculations. A cold-reference/warm-electron parallel-only scan
retains one Bessel miss across the whole traversal, confirming auxiliary-state
and physical-key consistency.

| Test / build | Final active N | Capacity | Hits | Misses |
| --- | ---: | ---: | ---: | ---: |
| growth/False | 54 | 271 | 693 | 17 |
| parallel/False | — | — | 219 | 1 |
| growth/True | 54 | 271 | 693 | 17 |
| parallel/True | — | — | 219 | 1 |

## Regressions and reproducibility

Baseline at task start: **125 tests passed**. The final suite passes **141
tests**, including all prior tests plus the 16 new scan tests. The clean build,
normal builds, checked scan/dispersion builds and examples pass. The pinned
pre-normalization legacy numerical regression remains exact. A second exact
regression pins pre-Stage-4A commit
`8331925c7953dcc68f29f14c1d7b78dd9198bd5a`: serial all-kinetic kpar/kperp mode
files are **byte-for-byte identical** in rho_ref and d_ref with legacy Z and
fixed128. No expected data was regenerated to make tests pass.

Reproduce from the repository root:

```bash
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/measure_mixed_scans.py \
  --output /tmp/MIXED_SCAN_RESULTS.json

./plume.e inputs/example/mixed_scan_warm_reference.in
./plume.e inputs/example/mixed_scan_cold_reference.in
./plume.e inputs/example/fully_cold_scan.in
```

The three supplied examples continue Alfvén/ion-cyclotron branches from
kpar*d_ref=.06 to .1 at fixed kperp*d_ref=.05, use finite electron mass,
write nine root rows, and request no unsupported supplementary diagnostics.
Python wrapper write/run tests cover the same three compositions.
GNU Fortran 11.4.0 on Linux x86_64 was used; the full toolchain string and
platform are recorded in the JSON. No commit or push was made.

## Remaining limits

- No exact perpendicular scans, including fully cold scans; single-point fully
  cold perpendicular calculations remain intact.
- No cold-containing plasma-parameter/radial scans, model transitions, collisions,
  heating, FPC or eigenfunctions. Existing all-kinetic capabilities are retained.
- Geometry does not extend any existing susceptibility limit. Exact parallel
  kinetic limits used here are the already supported mixed/cold paths.
- No branch identity guarantee across mode crossings, poles, strong damping or
  coarse failed intervals; reduce continuation steps and compare branches.
- Independent maps do not associate branches across wavevectors. Their supplied
  frequency window must contain and resolve the desired modes.
- Legacy Z retains its historical numerical limitations; accurate/default Z is
  recommended. Bessel tolerances retain their established interpretation and
  do not certify total susceptibility or root accuracy.
- Standard mode/map files retain their existing precision; use test probes or
  direct evaluations when examining smaller numerical differences.

## File inventory

The inventory is relative to the task-start working tree. Pre-existing data
outputs, plotting-script modifications, and unrelated untracked files are
preserved. Existing numerical reports/artifacts and protected numerical sources
are unchanged.

Modified:

- `.github/workflows/doc.yml`
- `.gitignore`
- `Makefile`
- `README.md`
- `ford_project.md`
- `input.md`
- `linfpclib/linfpc.py`
- `output.md`
- `src/disprels.f90`
- `src/functions.f90`
- `src/plume.f90`
- `src/vars.f90`
- `tests/plasma_dispersion_audit.py`
- `tests/test_adaptive_dispersion.py`
- `tests/test_auxiliary_normalization.py`
- `tests/test_cold_reference_mixed.py`
- `tests/test_fully_cold_dispersion.py`
- `tests/test_kinetic_core.py`
- `tests/test_mixed_dispersion.py`
- `tests/test_normalization.py`

Added:

- `inputs/example/fully_cold_scan.in`
- `inputs/example/mixed_scan_cold_reference.in`
- `inputs/example/mixed_scan_warm_reference.in`
- `src/scan_support.f90`
- `tests/MIXED_SCAN_REPORT.md`
- `tests/MIXED_SCAN_RESULTS.json`
- `tests/measure_mixed_scans.py`
- `tests/mixed_scan_probe.f90`
- `tests/scan_failure_probe.f90`
- `tests/test_mixed_scans.py`
