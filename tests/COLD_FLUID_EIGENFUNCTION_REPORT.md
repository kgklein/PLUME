# Stage 4B.2: cold-fluid and mixed-species moments

## Milestone 1 audit (before production changes)

`disprels::calc_eigen` already calculates kinetic moments from each species'
ordinary susceptibility, not from explicit thermal velocity-space quadrature.
After the existing Ex=1 field elimination, its `eigen_L` block computes current
from `susc(species,row,column)`, removes the equilibrium-drift density current,
and uses continuity to obtain density. The separate `heat_L` block and FPC
callers are outside this task. Legacy all-kinetic single-point and scan callers
must continue to use that original block, without arithmetic changes.

Stage 4B.1 introduced `electromagnetic_fields::calc_em_eigen` (refresh `disp`,
then solve the complete `lam`) and `calc_mixed_eigen` (existing scan argument
adapter). The latter currently assigns NaN to every species moment. Serial
`om_scan`, 2D `om_double_scan`, and `test_disp` already dispatch cold-containing
fields to this module. The smallest integration is a species-resolved moment
adapter after a successful field calculation; no dispersion/root/cache changes
or new eigenvector algorithm are needed.

### Output units and kinetic formulas

Let e=deltaEx/B0, a=vA/c, w=omega/Omega_R, kappa=k*d_R, u=U_s/vA.
The existing columns are

- E=deltaE/deltaEx and B=deltaB/deltaEx (Gaussian units), with Ex=1;
- v=deltaU_s/(c*e), for all three velocity components;
- eta=(delta n_s/n0_s)/e, a relative species density perturbation.

The kinetic current conversion in the existing block is

    F = -i*w*(Q_s/D_s)*(vtp_solver**2/betap_solver)*chi_s E
      = -i*w*a**2*(Q_s/D_s)*chi_s E.

It gives F=v+a*u*eta*z_hat, not just v when u is nonzero. Existing transverse
velocities are Fx,Fy. The longitudinal expression in `calc_eigen` is
`(Fz-vx*kappa_perp*u/delta)/(1+kappa_parallel*u/delta)`;
its density expression is `(kappa dot v)/(a*delta)`, with
`delta=w-kappa_parallel*u`. Algebraically the same current removal is

    eta = (kappa dot F)/(a*w)
    v_z = Fz-a*u*eta.

This compact form reuses the existing susceptibility-to-current moment method
and continuity; it does not introduce a different kinetic closure or duplicate
lengthy thermal formulas. It avoids removable intermediate divisions by delta.
The historical block will remain literally unchanged for exact regressions.
The new mixed-only implementation will use the compact physical form and be
validated against independent velocity-space moments, including finite drift.

For rho_ref, kappa=K/sqrt(beta_R*alpha_R) and a=vtp/sqrt(beta_R).
The old density factor `K/(vtp_solver*sqrt(alpha_R))` becomes kappa/a.
For auxiliary normalization, V**2/B=a**2 and K/sqrt(B)=kappa with alpha_*=1.
No physical reference temperature, anisotropy, species ordering, or auxiliary
scale remains in physical moments. Signed Q enters the current conversion;
D must be positive to interpret density relative to n0_s. Zero-density species
have no defined relative-density eigenfunction and will be marked unavailable.

### Independently normalized cold-fluid equations

Use exp(i*k.r-i*omega*t), B0 along z, k=(kperp,0,kpar), and signed g=mu_s/Q_s.
With the same field/output units as above,

    -i*delta*v = g*(E + v cross z_hat + a*u*z_hat cross B)
    delta*eta = (kappa dot v)/a.

Set H=E+a*u*z_hat cross B. The normalized momentum matrix is

    [-i*delta, -g, 0; g, -i*delta, 0; 0, 0, -i*delta] v = g*H.

The direct solution uses the factored denominator (delta-g)*(delta+g):

    vx = g*(i*delta*Hx-g*Hy)/((delta-g)*(delta+g))
    vy = g*(g*Hx+i*delta*Hy)/((delta-g)*(delta+g))
    vz = i*g*Hz/delta
    eta = (kappa_perp*vx+kappa_parallel*vz)/(a*delta).

This uses neither temperature nor susceptibility. The independent current test
will check `D_s/Q_s*(v+a*u*eta*z_hat) = -i*w*a**2*chi_s E` using the unchanged
production cold susceptibility. It must retain both the magnetic drift force
and the drift density current. Exact Doppler/cyclotron poles remain failures;
nearby values are evaluated without floors, damping or offsets.

Planned implementation: a separate precise, standalone species-moment module;
per-species status/NaN isolation; unchanged E/B method and historical kinetic
block; the current mixed scan adapter and single-point output will populate
existing moment slots. Heating, energy, FPC, pressure perturbations and new
propagation limits remain excluded. Tests will cover direct off-root response,
independent Gaussian velocity integration, physical and auxiliary normalization,
root/scan diagnostics, failure recovery, checked builds and exact legacy files.

## Implemented design and call paths

The audit design is implemented in the separately compiled, explicit-state
`src/species_moments.f90` module. It exposes `calc_cold_moments`,
`calc_kinetic_moments`, and array dispatch `calc_species_moments`. The cold
routine solves the momentum matrix above analytically and then applies
continuity. It never obtains velocity by multiplying susceptibility. The warm
routine uses the same current/continuity method already in `calc_eigen`, in its
compact physical coordinates derived above. The historical `calc_eigen` body
is unchanged; no second thermal/harmonic response implementation is introduced.

```text
single-point test_disp / serial om_scan / 2D om_double_scan
    -> electromagnetic_fields::calc_mixed_eigen
       -> unchanged calc_em_eigen -> existing disp -> complete E/B and species chi
       -> species_moments::calc_species_moments (original species order)
          tauS == 0 -> direct cold momentum and continuity
          tauS >  0 -> existing kinetic current removal and continuity
```

The scan loops need no changes: they already use the mixed adapter when any
species is cold. The adapter now fills their existing U/n arguments after
successful fields. Single-point `test_disp` uses the same adapter and prints
`ns`, `Us`, and a per-species status list, alongside the existing E/B blocks.
A valid field solution is not invalidated by a failed species calculation.
The species arrays `susc` and total `lam` are only read by the moment kernels;
probes verify the species susceptibilities are bitwise identical before and
after diagnostics. No additional dispersion or derivative evaluations are
introduced by moment calculation itself.

Only physical native d_ref inputs enter cold-containing production moments.
For kinetic currents, the old coefficient `vtp_solver**2/betap_solver` is a^2;
using that identity avoids both a cold reference division and unnecessary
auxiliary state. Converting rho_ref wavevectors/velocity as described above
reproduces the same warm physical moments. Physical species anisotropies and
temperatures are not overwritten. No reference species is reordered, changed,
or replaced. All-kinetic rho_ref/d_ref callers still execute the original
`calc_eigen` block, retaining its exact arithmetic and output.

### Current normalization and the essential drift term

In Gaussian units, with epsilon=deltaEx/B0,

    delta J_s = n_R*q_R*c*epsilon * (D_s/Q_s)
                * [v + a*u*eta*z_hat]
              = -i*Omega_R*w*epsilon*B0/(4*pi) * chi_s E.

Since `Omega_R*B0/(4*pi*n_R*q_R*c)=vA**2/c**2=a**2`, the test identity is

    (D_s/Q_s) * [v + a*u*eta*z_hat] = -i*w*a**2*chi_s E.

This fixes normalization and signed charge independently of any root/eigenvector
algorithm. Both `a*u*z_hat cross B` in momentum and `a*u*eta*z_hat` in current
are required. Omitting either loses the drifting oblique response. D is not
silently absorbed into eta: eta is relative to each species' n0, not n_ref.

## Numerical contract and failure isolation

All used inputs must be finite: w, kappa, a, E, B (cold), chi (warm), mu, Q, D,
and drift. Require a>0, mu>0, Q!=0 and D>=0. The cold kernel does not inspect
temperature, anisotropy, collision rate, or any global state. Array dispatch
validates the temperature selector, with exact zero and no small-temperature
threshold. The explicit warm current method does not need a thermal argument.

| Status | Meaning and response |
| --- | --- |
| 0 | Finite valid moments. |
| 1 | Invalid/nonfinite input or invalid species classification. |
| 2 | Exact zero frequency (the existing dispersion convention excludes it). |
| 3 | Exact cold Doppler pole, delta=0. |
| 4 | Exact cold cyclotron pole, delta=+/-mu/Q. |
| 5 | Nonrepresentable/nonfinite intermediate arithmetic or final result. |
| 6 | Zero-density species: relative-density and bulk moments of an absent population are undefined. |

Failure assigns quiet NaNs to that species' velocity/density block and returns
its status, not a fabricated zero. Other species and valid E/B remain available.
The mixed adapter's optional `moment_status(:)` is separate from its unchanged
field status. It warns with species index, status, and wavevector on failed
moment blocks. If fields fail, all moments remain unavailable. Heating arrays
remain NaN and heating requests remain rejected before evaluation.

Cold exact poles are reported without floors or damping. Nearby points are
calculated directly using `(delta-g)*(delta+g)`. In contrast, warm delta=0
need not be a pole: current and continuity determine finite moments when the
kinetic susceptibility is finite. The compact elimination cancels the old
removable intermediate delta divisions algebraically. No physical resonance
is removed or changed. IEEE exception flags/halting modes are saved and restored
around arithmetic; checked builds can therefore observe a numerical failure
without crashing or contaminating a subsequent successful evaluation.

Zero-density moments are deliberately unavailable, even though a test-particle
limit could be assigned by an additional convention. No such convention is
introduced here. Cold anisotropy placeholders may be zero, negative or NaN.
These are never used by the cold kernel.

## Independent physics validation

`tests/cold_moment_reference.py` supplies independently structured oracles:

1. A NumPy complex linear solve of the normalized 3x3 momentum matrix, followed
   by continuity. The Fortran cold routine uses a factored analytic solution,
   so the oracle does not reproduce its implementation algebra.
2. An independent 4x4 solve of susceptibility current plus continuity for
   kinetic moments. Its chi comes from the existing independent NHDS/Stix
   bi-Maxwellian implementation with SciPy Bessel/Faddeeva functions and 512
   harmonics, rather than from the Fortran production response.
3. Direct linearized Vlasov velocity-space quadrature, without deriving moments
   from a susceptibility. It solves

       [-i*(w-kappa dot v)-g*d/dphi]*(f1/(epsilon*f0))
           = -(g/a)*(E+a*v cross B) dot grad_v(log f0)

   for a drifting bi-Maxwellian, with velocities expressed in vA units. A
   Fourier gyroangle derivative, Gauss-Hermite parallel quadrature and
   Gauss-Laguerre squared-perpendicular-speed quadrature integrate f1 to density
   and velocity. Bulk velocity subtracts `u*eta*z_hat` from the perturbed flux.
   The fine grid uses 24x24 velocity nodes and 32 gyroangles; comparison uses
   16x16 nodes and 24 angles. This real-velocity contour is valid for Im(w)>0.
   Lower-half-plane tests instead use the independently Landau-continued
   NHDS/Stix oracle, not an incorrect real-axis damped-mode integral.

The direct cold test contains 192 ordinary off-root comparisons across both
builds: proton, electron, alpha and nonunit signed charge/mass/density species;
zero and both signed drifts; signed kpar; parallel/oblique kperp; and growing/
damped complex frequencies. Every velocity component and relative density is
compared. Momentum, continuity and the current/susceptibility identity are
checked independently. The production `calc_chi_cold` is used only for the
current check, never to generate cold moments. Additional tests approach each
Doppler/cyclotron resonance from positive/negative imaginary offsets of order
1e-3 through 1e-7, without regularization.

Six direct Vlasov reference configurations cover proton/electron/alpha warm
species, anisotropy, drift, both kpar signs and temperatures differing by orders
of magnitude. Off-root mixed susceptibility tests cover both imaginary
half-planes. Production eigenfrequency tests cover fully cold, warm-reference
mixed, cold-reference mixed and multiple interleaved cold/warm species, with
parallel/oblique and signed-kpar cases. Reordering non-reference species changes
which warm species supplies the auxiliary scale; moments remain in the original
physical species order and agree after the explicit test permutation.

## Measured results and tolerances

The deterministic results are in
[COLD_FLUID_EIGENFUNCTION_RESULTS.json](COLD_FLUID_EIGENFUNCTION_RESULTS.json).
Unless stated otherwise, scaled error means componentwise
`max(abs(actual-reference)/max(1,abs(reference)))`, not pure relative error
near zeros. Residuals use sums/norms of both sides plus a unit scale.

| Comparison | Maximum measured error | Acceptance |
| --- | ---: | --- |
| Cold fluid moments versus independent linear solve, ordinary points | 1.66e-15 | rtol=5e-12, atol=5e-12 |
| Cold momentum residual | 3.48e-15 | 5e-12 normalized |
| Production mixed moment continuity residual | 1.87e-16 | 3e-11 normalized |
| Cold moment current versus production cold susceptibility | 6.41e-16 | 2e-11 scaled |
| Warm moments versus independent NHDS/Stix + 4x4 solve | 9.72e-12 | rtol=3e-7, atol=2e-7 |
| Warm moments versus direct Vlasov quadrature | 2.19e-12 | rtol=1e-9, atol=1e-10 |
| Coarse/fine Vlasov quadrature difference | 3.17e-12 | 1e-9 scaled |
| New warm current removal versus original calc_eigen moments | 7.59e-16 | rtol=3e-10, atol=2e-10 |
| rho_ref/d_ref equivalent warm moments | 1.96e-16 | 3e-10 scaled |
| Auxiliary b_star=0.01,0.1,1,10 | 2.08e-11 | 2e-7 scaled |
| Non-reference species permutation | 6.43e-12 | rtol=2e-7, atol=2e-7 |
| Near-pole cold linear-solve comparison | 4.00e-11 | rtol=1e-8, atol=1e-7 |
| Printed scan moments versus freshly initialized single point | 4.78e-7 | 5e-6 scaled |
| Existing Stage 4B.1 E/B at identical frequency | 0 | 2e-14 scaled/relative |

Near-pole comparisons have a separately documented conditioning allowance:
large velocities and a near-singular reference linear solve amplify roundoff.
They are not evidence of finite values at exact poles. The independent kinetic
susceptibility and printed scan precision also limit comparisons; the tighter
cold/current and full-precision moment checks isolate the new implementation
from those effects. No existing susceptibility or root tolerances were changed.

## Single-point, scans and output compatibility

Actual executable serial and 2D scans cover four cold-containing composition
classes, both propagation directions and both normal/checked builds (32 scan
runs). All existing field/root/plasma columns and scientific precision are
preserved. Each printed moment is compared with a separately initialized
single-wavevector calculation at the printed omega and k, whose moments are
also checked against the independent oracles. Formatting explains the largest
printed comparison error. No scan coordinates, root-continuation state, Bessel
cache, or numerical backend selection are modified.

The first six columns remain native kperp/kpar/betap/vtp and complex frequency.
Columns 7--18 are B then E pairs. The next 6*nspec values are complex Ux/Uy/Uz
per species; the next 2*nspec are complex densities. Existing species parameter
columns follow. Eigen-only rows retain 18+14*nspec real values. The same loader
and writer conventions remain valid. No heating/FPC columns are added.

Single-point cold-containing `option=-1` prints the existing field blocks, then
`ns`, `Us`, and species statuses when the electromagnetic gauge is available.
Normalization sidecars now identify moment units, each species model, and
failure/status semantics. All-kinetic output layouts and calculations remain
unchanged. Modes with Ex=0 retain the Stage 4B.1 unavailable status, with no
alternative normalization. Exact perpendicular scan support remains deferred;
existing single-wavevector geometry limits are preserved.

Failure tests include zero frequency, exact Doppler/cyclotron poles, nonfinite
inputs, invalid physical parameters, absent species, unrepresentable arithmetic,
successful evaluations after failures in the same probe process, failed fields,
and failed scan points followed by successful continuation. A zero-density third
species produces a status-6/NaN moment block while the first two species and
E/B remain valid. A subsequent populated case succeeds. A real coarse scan
failure produces NaNs rather than stale moments and subsequently recovers.

The existing four field examples now demonstrate moments as well:

- `mixed_fields_warm_reference.in`: warm protons, cold electrons;
- `mixed_fields_cold_reference.in`: cold protons, warm electrons;
- `cold_fields_parallel.in`: fully cold parallel eigenfunction scan;
- `cold_fields_single.in`: fully cold oblique single-point eigenfunction.

The one new example `mixed_moments_multispecies.in` adds cold protons, drifting
warm electrons, warm alphas and cold alphas. It is charge/current balanced but
no such constraint is enforced. All five run successfully: the four scans each
produce three finite E/B/moment rows; the single-point case prints its valid
moments. Heating/energy/FPC remain explicitly unavailable.

## Regression, toolchain and reproducible commands

Baseline before changes: **154/154 passed** (102.472 s). Full suite after
integration: **169/169 passed** (113.622 s), including 15 new moment tests.
Both the pinned pre-normalization fixed128/legacy regression and complete
byte-for-byte historical eigen-only scan comparisons remain exact. The original
all-kinetic `calc_eigen` procedure body is unchanged. Stage 4B.1 assertions that
expected all mixed moment slots to be NaN were updated to require finite moments;
no established numerical expected values or regression tolerances were changed.
All their independent electromagnetic/nullspace/Faraday checks remain intact.

Toolchain: GNU Fortran 11.4.0, Linux x86_64, promoted default real/complex
precision; NumPy 2.2.6 and SciPy 1.15.3 for independent tests. The small new
module is compiled with the existing precise-module override
`-fno-fast-math`, so IEEE finite checks and failure handling survive the normal
production O3/fast-math flags. Checked builds use O0, `-fcheck=all`, and
`-ffpe-trap=invalid,zero,overflow`. Probes link the standalone cold moment kernel
without dispersion/root machinery, set global beta and vtp to zero, and verify
no physical/solver species arrays are allocated. Intel/other compilers were not
available for validation; no numerical kind or third-party backend is changed.

```bash
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
MPLCONFIGDIR=/tmp/plume-mpl OPENBLAS_NUM_THREADS=1 \
  python3 tests/measure_species_moments.py \
  --output /tmp/COLD_FLUID_EIGENFUNCTION_RESULTS.json
./plume.e inputs/example/cold_fields_single.in
./plume.e inputs/example/mixed_moments_multispecies.in
```

The measurement command separately passes all 15 new tests and writes its
observations without rewriting any regression expectations. `git diff --check`
passes. All work is local and uncommitted.

## Changed-file inventory

Changes are relative to the task-start snapshot; unrelated user edits and
outputs, including `plotter/example/cycle.sh`, are preserved.

| Files | Purpose |
| --- | --- |
| `src/species_moments.f90` (new) | Independent cold momentum/continuity, existing kinetic current method, status/isolation and species dispatch. |
| `src/electromagnetic_fields.f90` | Populate existing mixed moment slots; optional per-species statuses; E/B methods unchanged. |
| `src/disprels.f90` | Cold-containing single-point moment output only. |
| `src/vars.f90` | Moment units/model/status metadata only. |
| `Makefile` | New separate precise object and dependency. |
| `tests/species_moments_probe.f90`, `tests/moments_eigen_probe.f90` (new) | Standalone and production-call probes, including susceptibility preservation. |
| `tests/cold_moment_reference.py`, `tests/test_species_moments.py` (new) | Independent fluid/current/Vlasov references and 15 tests. |
| `tests/measure_species_moments.py`, `tests/COLD_FLUID_EIGENFUNCTION_RESULTS.json` (new) | Reproducible numerical measurements. |
| `tests/em_auxiliary_probe.f90` | Add moment measurements at explicit test-only auxiliary scales. |
| `tests/test_mixed_eigenfunctions.py` | Update obsolete mixed-moment availability assertions; retain independent field and exact regression checks. |
| `tests/plasma_dispersion_audit.py`, `tests/test_adaptive_dispersion.py`, `tests/test_auxiliary_normalization.py`, `tests/test_cold_reference_mixed.py`, `tests/test_fully_cold_dispersion.py`, `tests/test_kinetic_core.py`, `tests/test_mixed_dispersion.py`, `tests/test_normalization.py` | Explicit build/link lists include the new dependency; no numerical expectation changes. |
| `README.md`, `input.md`, `output.md`, this new report | User/developer conventions, availability, equations and validation. |
| `ford_project.md`, `.github/workflows/doc.yml`, `.gitignore` | Publish the report and expose the new example. |
| Four existing `inputs/example/*fields*.in` examples listed above | Update comments to describe enabled moments. |
| `inputs/example/mixed_moments_multispecies.in` (new) | The missing multi-cold/warm example. |

Task-start hashes confirm unchanged kinetic/cold susceptibility sources,
normalization adapters, harmonic selection/cache, Z/Faddeeva sources and vendor
code, root/scan support, main program, FPC, Python wrapper, and all earlier report
and measurement evidence. Exact procedure comparisons additionally confirm
unchanged `disp`, `disp0`, `disp1`, `disp_mixed`, `disp_all_cold`,
`disp_mixed_parallel`, `rtsec`, `zet_in`, `om_scan`, and `om_double_scan`.

## Remaining limitations and stopping boundary

This stage adds only bulk velocity and relative-density eigenfunctions, in the
existing Ex=1 gauge. It does not add pressure perturbations, heating, wave-energy
normalization, FPC, collisional cold responses, parameter/radial scans, or
perpendicular scans. Existing kinetic diagnostic behavior is preserved, including
its previously audited energy limitations. Real-velocity Vlasov quadrature does
not validate Landau continuation in the lower half-plane by itself; that part
uses the independent analytical susceptibility reference. Finite precision near
poles, extreme parameter ranges and ill-conditioned field gauges remains a
physical/numerical limitation with explicit unavailable statuses, not floors.

Stage 4B.3A should separately audit/restore kinetic-species heating availability
in mixtures; FPC follows in 4B.3B. No work on either stage is included here.
