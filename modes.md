<a id="calculation-modes-maps-and-continuation"></a>

# Calculation modes, maps and continuation

[Manual home](./) · [Input reference](input.md) · [Output schemas](output.md)

The executable reads one input file, builds species responses, and follows
`option`. The same production dispersion evaluator and secant root method serve
kinetic, cold and mixed systems. Their geometry restrictions differ; consult the
[capability matrix](models.md#capability-matrix) before requesting diagnostics.
Create the input's `data/<dataName>` directory and use distinct output names.

<a id="workflow-selector"></a>

## Workflow selector

| option | Operation | Root initialization | Supplementary output |
|---|---|---|---|
| -1 | Evaluate at supplied complex frequency, without refining it | guess_N | Text diagnostics and corrected .energy where valid; off-mode fields can be unavailable |
| 0 | Single-wavevector frequency map and root search | maps | .map and .roots; no scan eigen/heating/tensor blocks |
| 1 | Supplied roots or map initialization, then serial scans | use_map chooses maps or guess_N | .modeN; eigen/heating/tensor according to each scan |
| 2 | Two-coordinate grid, nscan=2 | use_map chooses maps or guess_N | .modeN grid blocks; special flag rules below |
| 3 | Deprecated historical hard-wired calculation | none | Stops; not a supported workflow |
| 4 | Repeat frequency maps along one scan, nscan=1 | maps at each point | Tagged .map/.roots files only |
| 5 | Experimental radial-model continuation | maps or guess_N | Historical radial files; all kinetic only |
| 6 | Historical gyrotropic FPC | maps or guess_N | Gyro velocity-grid files; all kinetic only |
| 7 | Historical Cartesian FPC | maps or guess_N | Projection/distribution/debug files; all kinetic only |

`use_map` does not change option=-1, 0 or 4 into a different workflow.
`option=1, nscan=0, use_map=.false.` refines supplied guesses and prints roots
when writeOut is true; it does **not** itself write a standalone .roots file or
scan eigenfunctions. Use a diagnostic scan or direct option=-1 for those outputs.

<a id="supplied-guesses-and-direct-evaluation"></a>

## Supplied guesses and direct evaluation

Provide `guess_1` through `guess_nroot_max` with `g_om=omega_r/Omega_R` and
`g_gam=gamma/Omega_R`. For root refinement, two nearby points are constructed
from the guess and passed to the existing complex secant method. These are not
guaranteed brackets. Several guesses can converge to the same branch or fail.
A negative gamma is a damped mode under exp(-i omega t).

Existing single-wavevector examples:

```bash
mkdir -p data/example
./plume.e inputs/example/cold_proton_electron.in
./plume.e inputs/example/mixed_cold_protons_warm_electrons.in
```

These are supplied-root workflows; inspect stdout and warnings. For option=-1,
`test_disp` evaluates the first supplied frequency as-is and prints the determinant
and requested/default direct diagnostics. Supplying nroot_max>1 does not make it
an arbitrary-frequency output scan. A finite determinant is not evidence of a
root; an off-mode matrix may correctly have unavailable eigenfunctions.

The raw API determinant is det(Lambda). Direct stdout for a warm reference
historically displays det(Lambda)*solver_vtp^6, while a cold reference displays
the raw determinant. [Numerical methods](numerics.md#tensors-and-the-determinant)
defines the matrices and scaling. Never compare these differently scaled numbers
as if they were the same residual.

<a id="frequency-map-search"></a>

## Frequency-map search

A map samples complex frequency grid points at fixed plasma and k. Local
minima of determinant magnitude seed secant refinement; finite successful
candidates are deduplicated using the historical absolute frequency gap 1e-5.
The internal candidate storage is 500, not a user-selectable input.
A minimum is not necessarily a physical zero, and a map cannot guarantee that
all branches were found. Map resolution and bounds require a convergence study.

The workflow distinguishes:

1. **Grid points:** every sampled determinant, written to .map, including invalid
   points as NaN on guarded paths.
2. **Refined candidates:** roots found and deduplicated before final selection,
   written to .roots. These are not grid frequencies.
3. **Selected roots:** subsequent filters include gamma/Omega_R<1 and, when
   positive_roots=true, omega_r/Omega_R>-1e-10.
4. **Continued branch slots:** scans follow nroot_max slots using the existing
   continuation; they do not continue every .roots row.

`nroot_max` is the supplied-guess/continuation count, **not** a cap on lines in
.roots. The file is written before the later map filters. The optional internal
sorting block is disabled by default: do not assume branches are automatically
ranked by least damping. A too-large requested count is not a way to manufacture
missing roots. For incomplete searches, inspect candidates and initialize a
bounded supplied-guess run deliberately. The historical all-kinetic legacy
backend also retains weaker failure handling than the guarded accurate path.

The existing [example_map_par.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_map_par.in)
is a secondary map lesson, not the introductory run:

```bash
./plume.e inputs/example/example_map_par.in
```

It writes `data/example/dispersion_map_par.map` and `.roots` (and any scans
requested by the file). Use the matching input to interpret units. The audit
found more .roots rows than nroot_max here; this is expected from writer order.

<a id="linear-versus-historical-logarithmic-maps"></a>

### Linear versus historical logarithmic maps

For linear flags, i=0..nr and j=0..ni give

\[
 w_{r,i}=\mathrm{omi}+i\frac{\mathrm{omf}-\mathrm{omi}}{\mathrm{nr}},\qquad
 w_{i,j}=\mathrm{gami}+j\frac{\mathrm{gamf}-\mathrm{gami}}{\mathrm{ni}}.
\]

All coordinates are omega/Omega_R. nr and ni count intervals; there are
(nr+1)(ni+1) grid points. Signed imaginary endpoints are ordinary linear bounds.

**Use linear map grids for physical endpoint control.** Historical logarithmic
maps have an unresolved implementation inconsistency, intentionally not changed
by the documentation or reader corrections. If the respective log flag is true,
the actual code instead evaluates

\[
 w_{r,i}=10^{\mathrm{omi}+i[\log_{10}(\mathrm{omf})-\log_{10}(\mathrm{omi})]/\mathrm{nr}},
\]
\[
 w_{i,j}=-10^{|\mathrm{gami}|+j[\log_{10}|\mathrm{gamf}|-\log_{10}|\mathrm{gami}|]/\mathrm{ni}}.
\]

The origins are exponent-like, while increments use physical-bound ratios.
Neither an exponent-bound convention nor a physical-endpoint convention
explains both. Real logarithmic endpoints must be positive and imaginary
magnitudes nonzero for these logarithms to exist; all imaginary input signs are
discarded and the result is negative. For omi=.1, omf=1 the range is approximately
1.258925..12.58925; for gami=-1e-4, gamf=-1e-2 it is -1.000230..-100.0230.
Warnings identify this behavior. No current numerical-mode switch repairs or
selects a different map convention. A future correction needs a separate
compatibility decision; [historical evidence](https://github.com/kgklein/PLUME/blob/main/tests/PLUME_5B0_CORRECTIONS.md)
records the originating source change.

<a id="serial-wavevector-scans"></a>

## Serial wavevector scans

Set option=1 and give nscan indexed scan_input_N groups. After initialization,
they run **serially**, passing the last state into the next scan; they are not
independent resets. Match the initial input wavevector and supplied root to the
first scan point, particularly for the unchanged kinetic-only initialization.
Cold-containing runs preflight and prepare their first native scan coordinate.

For an ordinary scan, M=ns*nres continuation intervals are traversed. The code
writes the endpoints and every nres-th point: normally ns+1 rows per branch.
Hidden substeps help continuation but do not add rows. Positive logarithmic
wavevector scans interpolate log10 of coordinates consistently; this is separate
from the defective logarithmic **frequency-map** convention.

| scan_style | scan_type | Meaning of swi/swf |
|---|---|---|
| 0 | 0 or 1 | Initial/final native kperp or kpar |
| -1 | 0 | Target vector (swi=kperp_final, swf=kpar_final); start is the current vector |
| -1 | 1 | Target angle in degrees in swi; swf unused; magnitude fixed |
| -1 | 2 | Target signed magnitude in swf; swi unused; angle fixed |
| 0 | 2 or 3 | Historical all-kinetic beta or vtp scan |
| s=1..nspec | 0..5 | Historical population-s tau, mu, alpha, Q, D or drift scan |

The compound angle is atan(kperp/kpar); its historical magnitude is signed
kpar/cos(angle), not a new quadrant-independent atan2 convention. Negative
parallel propagation is supported for nonzero points. Logarithmic scanned
coordinates must be positive, so use linear traversal for negative kpar.

Supported mixed/cold scans change wavevectors only. They keep temperatures,
reference normalization and drifts fixed. kpar=0 is rejected even at hidden
points; there is no small-kpar substitution. Exact-parallel kperp=0 scans work
for cold/mixed populations but remain unsupported for all-kinetic populations.

```bash
./plume.e inputs/example/mixed_scan_warm_reference.in
./plume.e inputs/example/fully_cold_scan.in
```

These write branch `.modeN` files. The historical cold-reference scan example
`mixed_scan_cold_reference.in` contains a failed written point in the current
local replay; it is useful for studying explicit NaNs/recovery, not as an
all-success introductory benchmark. The smaller
`mixed_fields_cold_reference.in` provides a successful cold-reference scan.

<a id="two-dimensional-wavevector-grids"></a>

## Two-dimensional wavevector grids

Use option=2, nscan=2. Supported coordinate pairs are distinct Cartesian
kperp/kpar axes, or ordered polar axes (scan_input_1 angle, scan_input_2 magnitude).
The first scan is outer and the second is inner. Each coordinate uses its own
ns/nres; output has blank separators between outer blocks. The existing large
[example_guess_double.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/example_guess_double.in)
is historical, not the missing small curated 2D lesson. It is not recommended
for a first run; see the [examples index](examples.md).

Eigenfunctions and heating in the grid use **scan_input_2** switches. Tensor
output requires `tensor=.true.` in **both** groups. The code preserves original
species order and native coordinates. Independent tests exercise small Cartesian
and polar grids; a dedicated small public input remains a curriculum follow-up.

<a id="root-failure-and-continuation"></a>

## Root failure and continuation

Successful points reuse the existing last-root continuation strategy. Guarded
scan refinement checks finite roots/tensors, status and a determinant scaled
by the maximum wave-matrix element. A failed written point has NaN root/diagnostic
slots, with parameters retained; a failed hidden point warns without adding a
row. The continuation seed remains the last valid one. A failed response must
not become a plausible converged root or contaminate subsequent species caches.

This is not a new branch-tracking algorithm. Near crossings, poles, strong
damping or genuine branch termination, forward and reverse traversal can choose
different branches. Check residuals, polarization and isolated single-wavevector
refinements; reduce step size deliberately rather than assuming file continuity
proves branch identity.

<a id="repeated-maps-and-historical-parameterradial-workflows"></a>

## Repeated maps and historical parameter/radial workflows

Option=4 remaps at each scanned point rather than continuing a root. It accepts
only one scan. Eigen/tensor flags and legacy heating are ignored with warnings;
corrected map heating is rejected. Files remain determinant maps and candidate
roots, without supplementary blocks. Use a separately initialized supported
root/scan workflow for diagnostics. No heating is added to maps.

All-kinetic plasma-parameter scans retain legacy meanings and input-dependent
units; corrected heating is restricted to wavevector scans. d_ref beta and
reference tau scans synchronize each other; independent 2D variation is rejected.
Temperature scans cannot reach or leave tauS=0. These are not supported cold/mixed
parameter scans. [Normalization](normalization.md#scan-paths) explains why a beta
scan follows different physics in the two conventions.

Option=5 is experimental radial interpolation. It reads `modelName0.rad` for
radius/beta/vtp and `modelNameN.rad` for each species' six parameters, and selects
one of seven k_range layouts using k_scan=0..6. It is all-kinetic, with legacy
radial heating only; corrected heating and cold species are not supported.
No small public standalone radial model is validated as a primary example.
The synthetic interpolation test is a developer fixture, not a promise of a
complete user workflow. Input fields and exact old record layouts are retained
in [Input](input.md#experimental-radial-controls) and [Output](output.md#experimental-radial-output).
