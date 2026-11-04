# Stage 3A: fully cold plasma dispersion

Historical report. [Stage 3B.2](COLD_REFERENCE_MIXED_REPORT.md) now supports cold-reference
mixtures containing warm species; historical measurements below are retained.

## Audit and implementation choice

The Stage 2 tensor and cold-kernel conventions are compatible without any
transpose, sign or frequency scaling changes. `susc(species,row,column)` is
ordinary dimensionless susceptibility; `lam` is the electromagnetic wave
matrix. Cartesian axes are `B0` along z and `k=(kperp,0,kpar)`.
Frequencies use `omega/Omega_R` and `exp(i*k.r-i*omega*t)`; positive imaginary
frequency means growth. The existing finite-electron-mass cold kernel reads
only density, mass, signed charge and parallel drift plus explicit frequency,
wavevectors and Alfvén speed. Its equations need no change.

The preimplementation audit identified these reference-temperature assumptions:

| Location | Thermal dependency and Stage 3A treatment |
| --- | --- |
| `vars::sync_solver_state` | Requires positive beta and reference temperature, constructs `sqrt(beta*alpha_R)`, `vtp*sqrt(beta)` and `beta/tau_s` for warm species. Both input reading and `disp` now dispatch all-cold states **before** this adapter. |
| `functions::spec_read` | Collision-frequency expressions contain thermal speeds; its existing exact-zero-temperature branch returns zero collision frequency before those expressions. Retained unchanged. |
| `disprels::disp0`, `disp1` | Thermal Bessel/resonance arguments use beta, reference alpha, species temperature and derived `solver_*` state. Neither is called for an all-cold state. Kinetic equations are unchanged. |
| `disprels::disp_mixed_parallel` | Uses positive warm-species thermal speeds; the all-cold branch precedes this selection too. Warm-reference mixed behavior is unchanged. |
| `test_disp` | Previously scaled printed determinants by `solver_vtp**6`. All-cold output now explicitly prints the raw determinant, without accessing that variable. Existing kinetic/mixed scaling is retained. |
| `calc_eigen`, `fpc` | Eigenfunctions/heating and velocity-space distributions use temperature ratios, beta and thermal velocity. Existing cold-containing entry guards remain in place; unsupported input options are rejected before entry. |
| `om_scan`, `om_double_scan`, `map_scan`, `radial_scan` | Thermal conversions, scan-bound heuristics and supplementary diagnostics remain inappropriate for cold reference species. Cold-containing parameter/radial scans stay rejected. |
| `map_search`, `refine_guess`, `rtsec` | Direct frequency evaluation and refinement use `disp`; no thermal normalization is intrinsic to these methods. Native input columns and existing singular-point handling can be reused. |
| Python input wrapper | Its warm-reference restriction needed a fully cold exception; reference beta consistency, collision/diagnostic restrictions and zero-temperature collision-frequency behavior remain intact. |

There is no unresolved mathematical convention and no need to rewrite the
kinetic susceptibility. The selected path adds only native cold wave-matrix
assembly and bypasses the thermal state rather than supplying a fictitious
reference temperature.

## Dispatch, validation and tensor assembly

`all_species_cold()` derives classification from the authoritative species
array, requiring every `tau_s==0` exactly. It is not a stored species-model flag.
`has_cold_species()` and all positive-temperature classification retain their
Stage 2 meanings; no small-temperature threshold is introduced.

At input reading and each dispersion evaluation, all-cold states take
`validate_all_cold` instead of `sync_solver_state`. Validation requires `d_ref`,
`betap=0`, positive `vtp=v_A,ref/c`, valid mass/charge/density, collisionless
operation and supported options `-1`, `0`, `1` with `nscan=0`. Used species and
wave/frequency inputs also pass the cold kernel's established finite-input
checks. Cold anisotropy is never used in a scale or susceptibility. The thermal
`solver_spec` array is not allocated by this path and no `solver_*` scalar is
read. No gyroradius or reference thermal speed is constructed.

A cold reference with any warm species receives an explicit Stage 3B error.
`rho_ref` remains all-kinetic only. Warm-reference mixed and legacy kinetic
states continue through the established adapter; its anisotropic conversion
remains `k_rho=k_d*sqrt(beta_R*alpha_R)`.

The small `disp` dispatcher delegates cold-containing evaluations to the
existing `disp_mixed` preflight/status wrapper. It calls `calc_chi_cold` for each
cold species using native `kpar`, `kperp`, `vtp`, and stores each response.
For all-cold states it then selects `disp_all_cold` instead of either thermal
assembly path. This routine constructs

\[
\epsilon=I+\sum_s\chi_s,\qquad
n_x=\frac{k_\perp d_R}{(v_A/c)w},\quad
n_z=\frac{k_\parallel d_R}{(v_A/c)w},
\]
\[
\Lambda=\epsilon+\mathbf n\mathbf n^T-(n_x^2+n_z^2)I.
\]

Identity is added once. The wave terms subtract `nz**2`, `nx**2+nz**2`,
`nx**2` from the three diagonals, and add `nx*nz` to xz and zx. The determinant
uses PLUME's existing six-independent-component expression. No cold-fluid
susceptibility equations are duplicated. Exact parallel, perpendicular,
oblique and zero-wavevector evaluations use the actual input components.

The secant update, convergence criteria, root refinement and map algorithms
are unchanged. Cold status codes remain 0 valid, 1 invalid input, 2 zero
frequency, 3 zero Doppler shift, 4 cyclotron pole, 5 nonfinite arithmetic.
Failure invalidates the tensor arrays and determinant; subsequent valid
calls recompute them. Map candidate selection skips isolated invalid samples
without replacing their physical outputs. No damping, frequency floors or
pole offsets are added. Existing shared messages containing “mixed” also apply
to fully cold evaluations.

Charge neutrality and equilibrium current checks remain warning-only, emitted
once during input reading. No densities, drifts or reference frames are changed.
Nonstandard homogeneous backgrounds may require external currents or forces.

## Independent verification

`fully_cold_probe.f90` calls the production `disp` and existing `rtsec`, emitting
all species tensors, all wave-matrix entries, determinant, status and roots.
Before input reading it poisons thermal solver scalars with NaNs and checks that
`solver_spec` remains unallocated after input and each evaluated frequency.
It is compiled against both production objects and a separate checked build:
`-O0 -g -fdefault-real-8 -fcheck=all -ffpe-trap=invalid,zero,overflow`.
This catches unintended visits to the thermal adapter, kinetic/Bessel routines
or thermal reference-scale arithmetic.

Independent references use the Stage 1 NumPy momentum/continuity linear solve,
including the perturbed magnetic force and drift density current, and literal
paper circular-response equations. They never call the Fortran cold kernel.
For stationary species the scalar checks retain finite electron mass:

\[
R=1-\sum_s\frac{A_s}{w(w+g_s)},\quad
L=1-\sum_s\frac{A_s}{w(w-g_s)},\quad
P=1-\sum_s\frac{A_s}{w^2},\quad S=(R+L)/2.
\]

Independent bracketed scalar roots verify `n^2=R,L` and `n^2=RL/S`. The ordinary
root follows `w^2=sum(A_s)+(kperp/a)^2`. Oblique/drifting roots use SciPy's root
solver on the independent full fluid determinant. Multiplication by `(a*w)^2`
conditions that test-only determinant away from poles; neither root solver is
added to production. Each final root is also checked using the unscaled wave
matrix normalized by its largest element.

The twelve new tests cover:

- Every species susceptibility, dielectric and wave-matrix element, plus the
  determinant, for proton/electron and multispecies systems, nonunit densities
  and masses, signed charges/drifts/wavevectors, `v_A/c=0.001, 0.03, 0.1`, and
  both frequency half-planes. No contribution is omitted or added twice.
- R/L roots at `|kpar*d_R|=0.2, 0.6`, in both propagation directions.
- O/X roots at exactly `kpar=0`, `kperp*d_R=0.2`, `v_A/c=0.1`.
- Oblique roots at 15, 45 and 75 degrees, magnitudes 0.12 and 0.4, and both
  signs of parallel wavevector.
- A charge-neutral, current-balanced alpha beam and its reversed-drift
  counterpart: `D_alpha=0.05`, `mu_alpha=0.25`, `Q_alpha=0.5`, `u_alpha=2`,
  `D_e=1.1`, `u_e=0.2/1.1`, `mu_e=1836.15267343`, `Q_e=-1`;
  protons have unit mass/charge/density and zero drift. At
  `(kpar*d_R,kperp*d_R)=(0.3,0.1)`, `v_A/c=0.01`, the unstable and damped
  conjugate roots are both recovered.
- Checked-build invariance with cold anisotropy placeholders 0, 1, -10, and
  `1e100`; zero beta and poisoned/unallocated thermal state throughout.
- Frequency, Doppler and signed cyclotron poles followed by valid evaluations;
  a failed zero-frequency root request followed by successful refinement.
- Warning-only nonneutral/current-carrying states, unsupported configurations,
  Python input construction and both executable examples.
- Maps containing poles (invalid samples retained as NaN), and a separate
  resolved perpendicular map that finds the X-mode root.

Well-conditioned tensor comparisons use componentwise `atol=rtol=5e-13`.
Direct determinants agree with NumPy at `rtol=5e-12`. Root real/imaginary parts
use `atol=2e-10, rtol=5e-10`, require successful refinement before 75 iterations,
and determinant residual `abs(det(Lambda))/max(abs(Lambda))**3 < 1e-10` in
both implementations. Rounded map-root records use absolute tolerance `5e-6`.
The existing cold-kernel near-pole conditioning tests remain unchanged; no
machine-precision claim is made at singular points.

## Results

GNU Fortran 11.4.0: `make -j2` succeeds. The full requested command

```sh
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

passes **50 tests**: all 38 existing tests and 12 Stage 3A tests. The pinned
pre-normalization regression remains **exact**. No expected numerical values
were altered. One existing Stage 2 invalid-input case now uses a cold reference
with a *warm* electron, preserving its unsupported-mixture assertion now that
an entirely cold plasma is supported.

| Measurement | Result |
| --- | --- |
| Maximum `abs(actual-reference)/max(1,abs(reference))`, every cold tensor element and both builds/oracles | `2.560e-13` |
| Maximum independent root difference, `omega/Omega_R` | `5.684e-14` |
| R, L at `kpar*d_R=0.2`, `v_A/c=0.01` | `0.220901945332439`, `0.180954505513887` |
| O, X at `kpar=0`, `kperp*d_R=0.2`, `v_A/c=0.1` | `428.624856188952`, `0.198913156806301` |
| Oblique alpha beam growing/damped roots | `0.304781932722338 +/- 0.0791683896443822 i` |

Stationary-mode imaginary parts are at roundoff level (approximately `1e-14`
or smaller). The production examples run successfully: the parallel L example
prints `0.058212` and the perpendicular X example prints `0.19994`, with
negligible imaginary parts. Python compilation and `git diff --check` pass.
Intel Fortran and rendered FORD documentation were not tested because those
tools are unavailable locally.

## Modified files and remaining restrictions

| File | Change |
| --- | --- |
| `src/vars.f90` | Derived all-cold classification, native validation, Stage 3B boundary, metadata. |
| `src/functions.f90` | Dispatch all-cold input validation before the thermal adapter. |
| `src/disprels.f90` | All-cold dispatch and native wave-matrix assembly; direct output without thermal scaling. |
| `linfpclib/linfpc.py` | Accept fully cold inputs, retain unsupported-mixture/diagnostic restrictions. |
| `tests/fully_cold_probe.f90` | Production evaluator/root probe with poisoned thermal state. |
| `tests/test_fully_cold_dispersion.py` | Independent tensors, roots, checked builds, maps, boundaries and examples. |
| `tests/test_mixed_dispersion.py` | Update only the cold-reference unsupported-mixture input, not expected numerical data. |
| `inputs/example/cold_proton_electron.in` | Minimal parallel L example. |
| `inputs/example/cold_proton_electron_perpendicular.in` | Exact perpendicular X example. |
| `.gitignore` | Include both new input examples. |
| `README.md`, `input.md`, `output.md` | Current support matrix, units, workflows and diagnostic boundaries. |
| `normalization.md`, `tests/MIXED_DISPERSION_REPORT.md` | Point historical audits to Stage 3A. |
| `tests/FULLY_COLD_DISPERSION_REPORT.md` | This audit, design and verification report. |
| `ford_project.md`, `.github/workflows/doc.yml` | Publish the Stage 3A developer report. |

The cold susceptibility kernel, kinetic equations, root method, FPC and heating
calculations are unchanged. No Makefile change is needed. Unrelated workspace
modifications and pre-existing untracked files were preserved.

Stage 3B must represent warm-species thermal speeds without dividing by a cold
reference temperature, and dispatch such mixtures before the current thermal
adapter. It needs independent mixed tensors and kinetic limits with cold
references while preserving the two existing warm-reference paths. Cold
fluid density/velocity eigenfunctions, heating/FPC, collisions, model-transition
scans and the warm exact-perpendicular limit remain separately scoped work.
None is enabled by Stage 3A.
