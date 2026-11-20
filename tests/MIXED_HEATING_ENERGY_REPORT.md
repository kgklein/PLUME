# Stage 4B.3A: corrected collisionless heating and complete wave energy

## Result and scope

The default is now `heating_convention='corrected'`. Collisionless all-kinetic,
fully cold, warm-reference mixed, and cold-reference mixed configurations can
request diagnostics at a single eigenfrequency or in the existing serial/2D
wavevector scans. Every kinetic species retains its total absorption/emission
diagnostic and all seven historical selected channels. Cold species have NaN
kinetic heating, but their susceptibility derivatives contribute to energy.

The previous [audit](MIXED_HEATING_ENERGY_AUDIT.md) and its
[measurements](MIXED_HEATING_ENERGY_AUDIT_RESULTS.json) are retained unchanged
as historical evidence. The new reproducible numerical observations are in
[MIXED_HEATING_ENERGY_RESULTS.json](MIXED_HEATING_ENERGY_RESULTS.json).
This implementation does not change susceptibility physics, normalizations,
the root method, Z evaluation, Bessel policy/cache, species moment equations,
or FPC. In particular, the separate collisional zz correction is preserved;
finite-rate xz/yz questions remain outside this task.

`heating_convention='legacy'` preserves the original `calc_eigen` implementation,
including its old normalizer, derivative, scratch behavior and output arithmetic.
Exact historical reproduction additionally requires `z_backend='legacy'` and
`bessel_mode='fixed128'`. None of these three controls selects the other two.
Previously supported all-kinetic collisional diagnostics require explicit
legacy mode. New corrected heating rejects collisional configurations.

## Mathematical definitions and normalization

PLUME uses exp(i k.r - i omega t), B0 along z, k=(kperp,0,kpar), and ordinary
species susceptibilities chi_s. Let w=omega/Omega_ref, kappa=k*d_ref,
a=vA_ref/c, and E=deltaE/deltaEx with Ex=1. The accepted complex mode supplies
the field polarization E. Corrected response snapshots and energy differentiation
use its **real laboratory frequency** wr, without replacing wr by a species
Doppler frequency.

Define chi_H=(chi+chi†)/2 and chi_A=(chi-chi†)/(2i). Then

```
H_E = E† E
B_r = (kappa cross E)/(a*wr)
H_B = B_r† B_r
D_s = d[wr*chi_H,s(wr,k)]/dwr
H_species(s) = E† D_s E
H_total = H_E + H_B + sum_s H_species(s)
W = abs(deltaEx)^2/(16*pi) * H_total
Q_s = E† chi_A,s(wr,k) E = Im(E† chi_s(wr,k) E)
Ps_s = Q_s/H_total
```

Here the subscript on Q_s denotes a diagnostic numerator, not the input charge
ratio also historically named Qs. The response remains in Gaussian field units.
Differentiating the dimensionless frequency product introduces no extra
Omega_ref factor: d[omega*chi]/domega = d[w*chi]/dw. No extra omega, 2*pi,
period, or factor of two is inserted in Ps.

For a weak isolated mode, -wr*sum(Ps_kinetic) approximately predicts gamma.
This is checked independently, never imposed. Ps is signed normalized
absorption/emission, not dimensional power or irreversible thermodynamic
heating. A cold collisionless real-frequency susceptibility is Hermitian away
from poles; its reversible energy contribution is not kinetic heating.

For rho_ref input, the diagnostic adapter constructs

```
kappa = K_rho/sqrt(beta_ref*alpha_ref)
a = vtp/sqrt(beta_ref)
tau_d,s = beta_ref/tau_rho,s
```

The perpendicular-reference anisotropy is essential. d_ref inputs pass through
unchanged. Physical reference mass/charge/density, drifts, species anisotropies,
frequency and k are fixed. Warm diagnostic snapshots reuse the validated
auxiliary thermal adapter; cold references never enter a thermal reference
conversion. An optional explicit auxiliary beta exists only in the reusable
Fortran API/tests, not as a new user input. It changes no physical energy.

### Seven selected channels

All numerators use the **same** wr and field E as Q_s:

| Channel | Numerator |
|---|---|
| 1 | abs(Ey)^2 Im chi_yy,n=0 (TTD) |
| 2 | Im[conj(Ey) Ez chi_yz,n=0] (TTD cross) |
| 3 | Im[conj(Ez) Ey chi_zy,n=0] (LD cross) |
| 4 | abs(Ez)^2 Im chi_zz,n=0 (LD) |
| 5 | Im[E† chi_n=0 E] |
| 6 | Im[E_perp† chi_n=+1 E_perp] |
| 7 | Im[E_perp† chi_n=-1 E_perp] |

Every channel is divided by the same H_total. Channels 1--4 sum to 5 when the
n=0 tensor structure applies; 5--7 are not an exhaustive partition of total
absorption. Summing all seven double-counts n=0. Higher harmonics and +/-1
longitudinal/cross terms are deliberately omitted. One independent finite-
kperp stress case has a 28.4487% remainder beyond channels 5--7; no artificial
closure is forced. Fixed128 corrected snapshots initialize **all** harmonic
entries, including the five unused n=0 entries. The legacy path is untouched.

## Implementation architecture and state isolation

```
read_in_params / scan preflight
  -> test_disp / om_scan / om_double_scan (only when heating requested)
     -> heating_diagnostics::calc_corrected_eigen
        -> existing calc_em_eigen at accepted complex mode
        -> existing calc_species_moments at accepted complex mode
        -> native_heating_state (explicit physical d_ref state)
        -> heating_energy::calculate_heating
           -> local response_snapshot at wr, including selected harmonics
           -> analytic cold / centered kinetic energy derivatives
           -> explicit component/fraction/availability calculation
     -> unchanged low_n=true mode columns + versioned companion writer
```

The new modules keep complexity outside the compiler-sensitive legacy kernel.
`heating_response` uses the shared kinetic includes and existing special functions,
not a second susceptibility formulation. Two scalar exact-parallel expressions
are extracted verbatim into small includes shared by the existing dispatcher and
the snapshot evaluator. Their mathematical operations are unchanged. The
optional harmonic output from `calc_warm_auxiliary` initializes scratch and
reuses the existing adaptive/shared kinetic implementation, including when
fixed128 is requested. Calls without optional outputs retain their established
production route. The Bessel cache policy and physical keys are unchanged.

All wr and derivative susceptibility/harmonic tensors are local arrays.
Derivative evaluations never call production `disp`; they cannot overwrite
`susc`, `susc_low`, `lam`, or `dispersion_status`. The probe verifies these global
arrays **bit for bit** around every calculation, including failed calculations,
checked builds and auxiliary-scale changes. `spec` is read-only and no solver
thermal state is changed. A valid root and E/B/moments survive unavailable
heating. Only legacy/FPC callers retain their historical global-state behavior.

### Frequency derivatives and numerical resolution

Warm derivatives use centered differences of the full product wr*chi_H:

```
C(h) = [(wr+h)*chi_H(wr+h) - (wr-h)*chi_H(wr-h)]/(2*h)
R(h) = [4*C(h)-C(2*h)]/3
```

Thirteen candidate relative steps range from 0.064 to 0.000015625, halving
at each step. Two consecutive adjacent Richardson differences must agree in
both matrix Frobenius norm and the E-contracted scalar at 1e-6, using max(1,
norm/value) comparison scales. Each species selects the best stable three-
estimate window; it does **not** blindly choose the smallest h. Full response
failure at any sampled point returns failure rather than silently retaining
an earlier estimate. No fallback to the historical normalizer is used.

The tolerance was selected against the independent high-precision audit, whose
ordinary energy target was approximately 1e-6 relative. The electron derivative
can have matrix norm approximately 2.73e7 while its E contraction is about
-179.1464554. Successive projected estimates become noisy at smaller steps;
a matrix-only test can therefore be misleading. At the ordinary kinetic
benchmark the chosen relative electron step is 0.004, not 1e-6, and its actual
species energy error is 8.95e-8 relative. Cold-reference warm electrons select
0.008, with actual species error 1.10e-7 relative. Selected warm-species steps
in the ordinary benchmark set span 0.000125--0.008.

The exported matrix/projected error estimates are conservative neighboring-
window changes; they may include the previous combined convergence bound.
The uncertainty in H combines projected adjacent-window differences and
floating-point error estimates. It is an empirical resolution estimate, not
a rigorous global error bound. Close poles, under-resolved kinetic sums or
Z-backend errors can defeat finite-difference validation; energy then remains
unavailable or needs a separate convergence study. In particular, selecting
legacy Z is allowed but does not guarantee modern absolute accuracy.

Cold derivatives differentiate the independently validated fluid linear solve
analytically, with no thermal parameters. For normalized velocity v and relative
density eta,

```
delta = wr-kpar*u,  g=mu/Q
M = [[-i*delta,-g,0],[g,-i*delta,0],[0,0,-i*delta]]
M v = g*(E+a*u*zhat cross B_r)
M v' = -g*a*u*(zhat cross B_r)/wr + i*v
eta' = kappa.v'/(a*delta) - eta/delta
D_s E = Hermitian_part{ i*Ds/(Q*a^2) * (v'+a*u*eta'*zhat) }
```

In code each matrix column is obtained before Hermitian symmetrization; the last
line specifies that matrix operation, not symmetrization of a vector. The
moment solver is reused with an effective forcing; cold chi equations are not
duplicated or changed. Off-center warm derivative snapshots skip cold response
entirely, avoiding spurious sampling of an unrelated cold pole. The central
real-frequency cold response must itself be regular.

## Independent validation and benchmarks

References are independently structured NHDS/Stix harmonic tensors, SciPy Z,
65-digit mpmath derivatives (including Z'=-2(1+z Z)), differentiated cold
analytical tensors, and a Gaussian fluid momentum/continuity solve. New tests
also use fresh independently differentiated drifting cases rather than only
replaying cached audit values. Arbitrary precision is test-only and optional;
without it the reference uses independent fourth-order differences. The
recorded run had mpmath available.

Representative H values, using the accepted mode field at kpar=.06,
kperp=.05, a=.01 except the noted growth/high-a cases:

| Plasma | Corrected H | Relative error against independent high precision |
|---|---:|---:|
| Warm protons/electrons | 20342.6800956026 | 7.95e-10 |
| Warm protons/cold electrons | 20384.8000697559 | 3.32e-11 |
| Cold protons/warm electrons | 20516.8791817050 | 1.02e-9 |
| Fully cold protons/electrons | 20607.6296161399 | below displayed precision |
| Interleaved cold protons, warm electrons/alphas, cold alphas | 115849.972219727 | 1.89e-11 |
| Weak anisotropy-driven growth (kpar=.3) | 59263.7354430223 | 1.01e-12 |

Both optimized and checked results are recorded separately in the JSON.
Maximum species energy scaled error in this benchmark set is 1.11e-7.
Additional drifting off-root real-frequency cases, both kpar signs, agree with
fresh high-precision species energies to 2.69e-9 scaled and raw kinetic Q to
3.97e-15 scaled. These off-root calculations explicitly qualify Ps when weak-
mode consistency is not applicable, without hiding raw response diagnostics.

| Check | Maximum observed discrepancy/result |
|---|---:|
| Independent seven raw channel numerators, scaled by max(1,abs(reference)) | 5.83e-11 |
| Stationary cold field/fluid energy identity, scaled | 3.99e-16 |
| rho_ref/d_ref energy comparison, both Z backends and fixed128/adaptive | bit-identical in tested cases |
| Auxiliary beta .01,.1,1,10, energy/fractions scaled | 4.35e-9 |
| Scan versus independently initialized point diagnostic | bit-identical for tested energy values |
| Weak-mode gamma consistency, relative | 6.97e-5 |
| Cold-containing scan successful diagnostic records | 58 |
| Explicit failed-root records / subsequent recoveries | 2 / 2 |

Single-point, serial and 2D tests include all supported compositions and
interleaved ordering, drift, signed kpar, exact-parallel and oblique response.
The pure parallel channels have their existing analytical zeros; no kperp
floor is introduced. High-harmonic omission is independently recognized.
All five example inputs run to completion, without root or energy warnings.

### Vacuum electric energy and signed fractions

For a stationary finite-mass cold test with a=.3, mass ratio 4, kpar=.2,
kperp=.1 and wr=7.46110067295017,

```
H_total = 53675.3078351476
H_E = 26837.6690117475
f_E = f_E_abs = 0.500000281212615
```

Thus omitting vacuum electric energy would produce approximately 50% error.
At the ordinary kinetic low-frequency benchmark, f_E is only 4.92195e-5.
No universal small-a/low-frequency assumption removes this term.

Fractions are f_E=H_E/H, f_B=H_B/H, f_species=sum(Hs)/H and f_s=Hs/H;
they are signed and unclamped. Synthetic cancellation/negative-energy tests
check values outside [0,1], exact zero H and finite-but-unresolved H. Signed
fractions are NaN when abs(H)<=10*uncertainty. The noncancelling measure

```
f_E_abs = H_E/(abs(H_E)+abs(H_B)+sum_s abs(Hs))
```

is available when its finite denominator is nonzero, even if signed fractions
are unresolved. It is **not a physical energy fraction**. Ordinary resolved
fractions sum to one; cold species retain their individual signed Hs/f_s.

## Failure, qualification and compatibility

A numerical status is independent of applicability and root/field validity:

| Status | Meaning |
|---|---|
| 0 | Resolved formal energy |
| 1 | Invalid diagnostic inputs |
| 2 | Unavailable/failed Ex=1 field representation |
| 3 | wr=0 |
| 4 | Response failure or real-frequency cold pole |
| 5 | Derivative not resolved |
| 6 | Vanishing/unresolved H |
| 7 | Unsupported corrected collisions |
| 8 | Failed scan root |
| 9 | Nonfinite/unrepresentable arithmetic |

Warm normalized Ps and all seven channels are NaN on numerical failure.
Cold/absent and numerical/applicability species statuses are distinct. Raw
Q/channels and partial formal components may still be available in the API
when a later normalization fails; consumers must check statuses. Nonfinite
arithmetic results are explicitly assigned NaN. Negative H is not a failure.
IEEE exception flags/halting settings are restored to the caller.

Applicability bit 1 flags abs(gamma/wr)>0.05. Bit 2 flags disagreement between
the independently computed weak-mode gamma and accepted gamma, beyond rounding/
derivative uncertainty and a 10%-of-combined-magnitudes criterion. These are
conservative warnings, **not** a proof of applicability or a comprehensive
near-degeneracy classifier. Formal resolved energy remains exposed, while
normalized kinetic Ps/channels are unavailable for flagged modes. Raw Q remains
available. Reactive benchmarks have H approximately -9993.59 (mixed) and
-9651.50 (fully cold), strong growth and both bits set; no energy closure is
forced. A finite real-frequency formula is not presented as an exact conserved
energy of such modes.

Tests exercise zero wr, exact real-frequency cold poles, invalid Ex forcing,
zero/unresolved energy, negative energy, strong/reactive modes, recovery after
failure, and checked floating-point traps. A test-only oscillatory perturbation
to Z deliberately makes the derivative nonconvergent: status 5 is required,
valid E/B/moments remain finite, globals remain bit-identical, and the next
unperturbed calculation in the same process succeeds. This perturbation is
never used as a physical reference or production numerical method.

One cold-reference 2D continuation point at (kperp,kpar)=(.052,.06) fails the
existing root residual validation (flag -2) in both normal/checked tests. It
is not relabeled converged: the companion records status 8 and NaNs, then the
next point recovers. Root algorithms/tolerances are unchanged. Serial corrected-heating records
also retain the accepted real frequency rather than applying the historical
output-only abs(wr)<1e-15 rounding, keeping wr and snapshots consistent; legacy
and non-heating output rounding are unchanged. This is an
explicit scan limitation, not a diagnostic accuracy claim for a failed root.

The historical `calc_eigen`, `disp0`, `disp1` and `rtsec` bodies remain unchanged.
Pinned pre-normalization regressions and existing full-file scan comparisons
remain exact with explicit legacy heating/Z/fixed128. No golden numeric values
or prior accuracy tolerances were changed. Tests that intentionally characterize
legacy defects explicitly select legacy; the old Python reader-offset defect
assertion is now a correctness assertion. The historical audit JSON is not
regenerated or rewritten. Legacy numerical issues remain available for replay,
not silently corrected in that convention.

FPC source/behavior, collisional tensor includes, kinetic/cold equations, field
and moment kernels, Bessel selector/cache, normalization equations and root
method are unchanged. Protected source hashes are recorded in the results.
Corrected collisional diagnostics, new FPC availability, maps, radial/plasma-
parameter corrected heating, and exact perpendicular warm-species propagation
remain deferred.

## Output and public interfaces

`heating_convention` is an optional `&params` character input, default corrected.
For scans, existing `heating` switches request it; in 2D, scan_input_2 selects
the diagnostics. `option=-1` uses corrected fields/moments/heating at a valid
supplied mode. Its historical four-entry stdout Ps_split block is retained;
all seven are available in the companion. No heating is added to maps.

The numeric low_n=true mode files retain column order and output precision:
root/plasma prefix; optional fields and moments; **all species totals**, then
**seven channels per species** (+1 before -1); species parameter suffix.
Cold heating is NaN in its existing slots. Corrected low_n=false explicitly
marks historical unassigned Ew as NaN when that old slot exists; heating-only
no longer emits the orphan Ew format-reversion record. These unassigned
historical values are never repurposed as corrected energy.

The companion `data/<dataName>/<outputName>.energy` has schema
`plume-energy-v1`. Eleven integer context/status fields precede sixteen scalar
real diagnostics, then fifteen reals per species, then a quoted source-file
string. Exact field order and status definitions are documented in
[output.md](../output.md#corrected-energy-companion-plume-energy-v1) and its
`# columns:` header. It contains native-unit wavevectors, accepted frequency,
H components, fractions, raw Q, normalized Ps/channels, derivative resolution,
and numerical/applicability/model availability flags. Scan context includes
workflow, scan, outer index, internal point index, mode and root_valid, so
missing/failed modes cannot shift later associations. `point` is not an output
row counter when nres>1. Single points use source `stdout`; scan records name
the actual mode file. Single/scan runs reset their prior companion; legacy
writes no new companion. Existing data directories must exist as for ordinary
PLUME output.

Python set_params/write/load understands heating_convention. The new
`load_plume_energy` parses versioned metadata, numeric columns, failed rows and
quoted paths. The existing low_n=true sweep loader now uses correct heating-only
offsets and seven-channel strides; marker tests check exact column identity.
Its ambiguous autodetection remains unchanged; supply nspec/eigen/heating when
needed. low_n=false uses the documented explicit schema or companion, not that
seven-channel parser.

## Build, test and reproduction

Recorded toolchain: GNU Fortran 11.4.0 (Ubuntu), NumPy 2.2.6, SciPy 1.15.3,
mpmath 1.4.1, Linux x86_64. New IEEE-sensitive modules use production precision
`-fdefault-real-8` and `-fno-fast-math`. Legacy compilation flags are retained. A repeated clean parallel build exposed
the existing implicit functions.mod ordering in disprels; the Makefile now
explicitly depends on functions.o, preventing that build race.
The final clean rebuild produces byte-identical objects and executable to
those used for the passing full suite; this check is recorded in the JSON.
Checked probes/executables use `-O0 -g -fdefault-real-8 -fcheck=all
-ffpe-trap=invalid,zero,overflow`. No runtime Python/high-precision dependency is
introduced.

```bash
make clean
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
# With optional mpmath importable for independent 65-digit references:
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/run_corrected_heating.py \
  --require-high-precision --output tests/MIXED_HEATING_ENERGY_RESULTS.json
./plume.e inputs/example/heating_corrected_cold_reference.in
./plume.e inputs/example/heating_legacy.in
```

The pre-change baseline passed 183 tests. After implementation a clean normal
build, all original tests, exact legacy regressions, and the new corrected tests
pass. Final full discovery: **201 tests passed, zero failures/skips, in 137.858 s**
(183 prior tests plus 18 new corrected-diagnostic tests). The exact legacy
regressions passed. The measurements runner independently reran the 18 new
tests with optimized and checked probes. Count/time are also in the results JSON. Both optimized and checked validations are part of discovery. The
corrected path intentionally makes up to 26 extra warm snapshot evaluations
per energy calculation; this stage does not optimize performance.

## Complete file inventory for this task

Pre-existing staged/unstaged changes were preserved. Relative to the task-start
workspace, this task adds or changes only the following files:

- `src/heating_response.f90`, `src/heating_energy.f90`, `src/heating_diagnostics.f90`
  (new explicit snapshots, derivative/energy/heating API, production adapter/IO).
- `src/kinetic_parallel_response.inc`, `src/kinetic_parallel_longitudinal.inc`
  (new shared, verbatim existing scalar equations).
- `src/cold_reference.f90` (optional selected harmonic/auxiliary-scale output).
- `src/disprels.f90` (requested diagnostic dispatch/writers, shared parallel
  includes; historical heat/root/thermal susceptibility bodies unchanged).
- `src/functions.f90`, `src/scan_support.f90`, `src/vars.f90` (input, supported
  configurations, default control and metadata).
- `Makefile` (objects/dependencies/IEEE flags).
- `linfpclib/linfpc.py` (control validation/roundtrip, companion reader,
  heating-block offsets/stride).
- `tests/corrected_heating_probe.f90`, `tests/heating_fraction_probe.f90`,
  `tests/test_corrected_heating.py`, `tests/run_corrected_heating.py` (new probes,
  independent validation, state/failure/IO tests and reproducible measurements).
- `tests/plasma_dispersion_audit.py`, `tests/test_adaptive_dispersion.py`,
  `tests/test_auxiliary_normalization.py`, `tests/test_cold_reference_mixed.py`,
  `tests/test_fully_cold_dispersion.py`, `tests/test_kinetic_core.py`,
  `tests/test_mixed_dispersion.py`, `tests/test_mixed_eigenfunctions.py`,
  `tests/test_mixed_scans.py`, `tests/test_normalization.py` (test link
  dependencies and explicit legacy input for historical tests; old executables
  receive their original namelist fields).
- `tests/test_heating_energy_audit.py` (reader correctness assertion; retained
  legacy/audit expectations).
- `tests/MIXED_HEATING_ENERGY_REPORT.md`, `tests/MIXED_HEATING_ENERGY_RESULTS.json`
  (new report/measurements; original audit evidence untouched).
- `inputs/example/heating_corrected_kinetic.in`,
  `inputs/example/heating_corrected_warm_reference.in`,
  `inputs/example/heating_corrected_cold_reference.in`,
  `inputs/example/heating_corrected_cold.in`, `inputs/example/heating_legacy.in`
  (working examples).
- `.gitignore` (include new examples), `README.md`, `input.md`, `output.md`,
  `ford_project.md` (focused current-convention and schema documentation).

## Remaining limitations

Real-frequency energy and Ps remain weak-mode diagnostics. Flags cannot resolve
every mode degeneracy; legacy Z accuracy or under-resolved harmonics can pass a
local derivative check, so independent Z/Bessel convergence still matters.
Ill-conditioned derivatives are unavailable instead of silently repaired.
The empirical error estimator is not a proof of accuracy for arbitrary plasma
parameters. Collisional corrected heating/operator validation and FPC require
separate tasks; cold kinetic heating will remain unavailable rather than zero.
No later stage is implemented, and no changes are committed or pushed.
