# Stage J.1: adaptive Bessel harmonics and caching

## Result and checkpoint record

Adaptive harmonic selection is the default. Each warm species has its own
cutoff, reusable coefficient arrays and work storage. This applies to both
input normalizations, warm-reference mixtures, auxiliary-normalized
cold-reference mixtures, and the existing weakly collisional all-kinetic path.
Cold susceptibilities and analytical mixed parallel limits remain unchanged.

The checkpoints were completed in order:

1. **A:** the unchanged baseline passed 73 tests. The standalone selector then
   passed three tests against independently summed SciPy Bessel tails, including
   a checked GNU build. Production susceptibility evaluation was still unchanged.
2. **B:** four isolated cache tests passed before production integration. These
   covered coefficients, adjacent orders, reuse, resizing, tolerance changes,
   species-count changes, fixed counts, capped warnings and recovery.
3. **C:** production dispatch, controls, metadata and wrapper support were added.
   Twelve production tests now cover tensors, roots, cache keys, auxiliary
   invariance, ordering, collisions, diagnostics/scans, input validation,
   arithmetic-failure recovery, and wrapper round trips. The complete suite
   passes **92 tests**, including all 73 pre-existing tests and the exact pinned
   pre-normalization regression. Existing numerical expectations and tolerances
   were not relaxed.

The new example `inputs/example/adaptive_bessel.in` runs successfully and selects
263 proton harmonics and 7 electron harmonics. Benchmark measurements are saved
in [ADAPTIVE_BESSEL_BENCHMARK.json](https://github.com/kgklein/PLUME/blob/main/tests/ADAPTIVE_BESSEL_BENCHMARK.json).

## Source audit and integration architecture

The audit used the normalization and Stage 0/1, 2, 3A, 3B.1 and 3B.2 reports.
Before this change, `disprels::disp0` and `disp1` summed orders `-128..128`, with
Bessel orders through `+/-129` for adjacent-order combinations. The collisionless
species expressions were already shared in `kinetic_species.inc`.
`kinetic_core` used the same expressions but recalculated all Bessel values for
each auxiliary evaluation. Legacy caches used the reference normalization and
individual thermal ratios as keys, not just the physical Bessel argument.

The new `adaptive_kinetic` module receives the existing converted species and
`kinetic_normalization`. It calls the cache with an independently calculated
physical argument, then evaluates the shared susceptibility expressions at the
selected range. `kinetic_collisional_species.inc` extracts the existing `disp1`
species block verbatim apart from trailing whitespace; it is included in both
its original context and the new runtime-range evaluator. Neither species
block changes the Vlasov-Maxwell equations, drift terms, harmonic ordering,
normalization factors, low-order heating contributions, or Z function.

`disp` and the root finder are unchanged. The new branch is at the beginning of
`disp0`/`disp1`, after their declarations. Warm-reference states still pass
through `sync_solver_state`; cold-reference warm states still use the validated
auxiliary adapter. `cold_reference::calc_warm_auxiliary` selects the new evaluator
when appropriate, passing the same native physical argument regardless of
auxiliary scale. It never constructs a thermal scale for a cold population.
Existing cold preflight, mixed failure handling, native wave assembly and
exact-parallel limits are retained. Identity and electromagnetic terms are
still added once. No new propagation limit or parameter-scan capability is added.

### Fixed compatibility and compiler sensitivity

`bessel_mode='fixed', bessel_nfixed=128` retains the historical compiled
susceptibility context, old Bessel evaluator and old caches. Other fixed counts
use the runtime-range evaluator and vector coefficients. All fixed counts mean
maximum absolute order N, i.e. `2*N+1` susceptibility harmonics, not N total terms.

An initial attempt to dispatch adaptive calculations inside the small `disp`
routine changed GNU optimizer decisions in unrelated eigenfunction arithmetic.
The pinned regression detected last-bit density/velocity differences. Moving the
branch into the large, non-inlined susceptibility routines restored exact results.
The legacy assertion remains exact against commit
`c6c7616e1fc58eb6118e1f9423a3f0df6d52ef20`.

New numerical modules use default real/complex kinds, matching PLUME's
compiler-promoted 64-bit arithmetic, and compile with trailing
`-fno-fast-math` (Intel rule: `-fp-model precise`). This preserves IEEE validation,
compensated summation, and exception handling. Legacy compiler flags are unchanged.
The independent Stage 3B.1 extraction remains test-only: it strips the new dispatch
and continues to test the original fixed-128 equations. Production has no Python
source-extraction dependency.

## Physical argument and convergence criterion

With native `d_ref` temperature theta_s and wavevector kappa_perp,

\[
 \lambda_s=\frac{\kappa_\perp^2Q_s^2\alpha_s\theta_s}{2\mu_s}.
\]

For historical `rho_ref`, this is
`kperp_rho**2 * Qs**2 * alphS / (2*muS*tauS_rho*alphS_ref)`.
The reference perpendicular gyroradius factor is retained. A change in auxiliary
beta or normalization anisotropy does not change the physical argument.

Let `p_n=exp(-lambda)*I_|n|(lambda)`. These positive symmetric weights have
sum one, second moment lambda, and generating function

\[
 M(t)=\sum_{n=-\infty}^{\infty}p_n e^{nt}
     =\exp[\lambda(\cosh t-1)].
\]

For `m=N+1` and any `t>0`, bounding both positive and negative tails gives

\[
 E_0(N)\le B_0=2\exp[\lambda(\cosh t-1)-mt],
\]
\[
 E_2(N)\le B_2=B_0(\cosh t+\lambda\sinh^2t).
\]

The second expression follows by differentiating M twice and dividing by lambda.
Choose `t=asinh(m/lambda)`, which minimizes B0. This gives conservative bounds,
not a claim to choose the smallest N satisfying the exact tails. The implementation
finds the first crossing of both bounds in their convergent decreasing region.
No test of a single small low-order coefficient can stop the selection prematurely.

Bounds are evaluated logarithmically. With `h=hypot(lambda,m)`, use
`lambda*(cosh(t)-1)=m*m/(h+lambda)` and
`log(B2)=log(B0)+log(h+m*m)-log(lambda)`.
The implementation avoids overflow in `hypot`, `h+lambda`, and `m/lambda`;
`asinh` uses a logarithmic expression for very large ratios. Probabilistic bounds
are clipped above at one. There is no subtraction of a nearly complete sum from
unity, no tolerance floor, and no Bessel approximation inside selection.
At exactly lambda=0, the selector returns N=0 and zero tails analytically.
Existing production parallel limits handle their finite n=+/-1 response directly.

Representative selected cutoffs (a cap of 2000 was used for this table):

| lambda | tol=1e-6 | tol=1e-10 | tol=1e-14 |
| ---: | ---: | ---: | ---: |
| 0 | 0 | 0 | 0 |
| 1e-8 | 1 | 2 | 2 |
| 0.01 | 3 | 5 | 6 |
| 1 | 9 | 12 | 15 |
| 10 | 21 | 27 | 32 |
| 100 | 61 | 76 | 88 |
| 1000 | 190 | 236 | 273 |
| 10000 | 602 | 745 | 863 |

The default 512 cap deliberately warns for some broad distributions in this table.
A user must increase it to resolve such calculations.

## Coefficients and cache

The historical `bessim0/1` use polynomial approximations with limited accuracy;
`bessim(n,x)` starts its downward recurrence from an order depending on n alone.
At large x this does not reliably resolve low orders of the broad distribution.
Simply increasing the susceptibility count therefore does not fix its coefficients.
These routines remain unchanged for fixed-128 reproduction.

The dynamic cache uses the same modified-Bessel recurrence, evaluated for the
whole array. Downward ratios satisfy

\[
 r_n=I_n/I_{n-1}=\lambda/(2n+\lambda r_{n+1}).
\]

The starting order resolves both positive tail bounds below
`0.01*epsilon(real)` and adds 32 orders of recurrence depth. Products of ratios
produce `I_n/I_0`; a compensated positive sum normalizes all coefficients with
`Gamma_0+2*sum(Gamma_n)=1`. No unscaled exponentially large I_n is formed.
Coefficient tests cover lambda=0, 1e-300, 1e-12 through 1e6, including adjacent
orders. Maximum relative discrepancy from SciPy `ive`, for reference values
above 1e-200, was **4.80e-13**; tests use `atol=2e-16, rtol=2e-12`.

Each cache entry contains the exact physical argument, policy/tolerance/limit,
active cutoff, allocated capacity, valid computed range, coefficients through
N+1, adjacent averages, Z work arrays, convergence bounds/status, and counters.
Capacity grows geometrically and is retained when N shrinks. Only the needed
coefficient range is calculated; unused retained capacity is not mistaken for
valid coefficients. Growth allocates replacement storage before discarding old
storage, so failure cannot leave a partly allocated reusable entry.

Keys use exact equality, without approximate matching. Frequency, parallel
wavevector, density, drift and auxiliary scale do not enter them. Changes in
physical perpendicular temperature, wavevector, anisotropy, mass or charge
cause recalculation when they change lambda. Tolerance changes reselect the
range; coefficients are reused when their computed range suffices. Species-count
changes preserve surviving slots, and reordered slots cannot reuse a different
argument. Cold species allocate no coefficient/work arrays.

Counters distinguish a complete request hit, a request miss, and a coefficient
array evaluation. A miss need not evaluate new coefficients. Fixed-128 cache
instrumentation is deliberately not retrofitted into its compatibility path;
its optional metadata labels detailed statistics unavailable, while evaluating
the final physical tail bounds for its unresolved flag.

## Caps, failures and metadata

Default controls are adaptive mode, tolerance 1e-10, maximum 512, and fixed count
128. Both counts must be positive and at most 46338 with the current int32
`n*n` expressions. Tolerance must be finite and strictly between zero and one.

At an unresolved adaptive cap, evaluation continues with **exactly the cap**.
The warning includes species, lambda, tolerance, maximum and both bounds.
Repeated unchanged requests hit the cache and do not repeat it. A later changed
request may warn again. Unresolved history survives later converged points and
is written to the `.bessel` metadata file. With `bessel_diagnostics=.true.`, the
file also includes per-species active counts, capacity, bounds and counters.
Normal numeric output layouts are unchanged. Standalone callers can request the
same table through `write_bessel_stats` and the sidecar through
`finish_bessel_diagnostics`.

Arithmetic/allocation/range failure is distinct from a capped valid array.
The dynamic evaluator returns status 5 and invalid tensors on arithmetic failure,
saving/restoring IEEE flags and trapping modes. A checked test deliberately
requests an unrepresentable recurrence order, then verifies an identical valid
response from the same cache. No artificial damping, frequency offset, temperature
floor, or Bessel-argument floor is introduced.

## Susceptibility and root validation

Every complex susceptibility and wave-matrix element is compared across proton,
electron and alpha populations; isotropic/anisotropic and drifting states;
positive/negative parallel wavevectors; both frequency half-planes; wavevectors
with widely separated species arguments; both normalization conventions; and
warm/cold reference mixtures. Additional tests exercise weak collisions,
all-kinetic eigenfunctions/heating and existing wavevector scans. Independent
physics uses the existing SciPy NHDS/Stix implementation with 1024 harmonics.

Maximum scaled discrepancy `abs(a-b)/max(1,abs(b))` for adaptive versus the
1024-harmonic evaluator in the broad tensor grid was **2.54e-11**. Including
auxiliary-scale, ordering and collisional comparisons, the aggregated maximum
was **7.90e-10**. Auxiliary scales 0.01, 0.1, 1 and 10 reuse the same physical
coefficient entry. Cached and freshly computed responses at identical states
compare exactly in the direct cache-key test.

Well-conditioned independent comparisons retain `rtol=2e-7, atol=1e-5`.
Cancellation-sensitive stress points expose the unchanged Z implementation's
`error=1e-7` stopping rule: the largest scaled independent discrepancy was
**7.63e-4**. The high-count fixed evaluator reproduces it; at representative
small/moderate arguments fixed 128 also reproduces it. Stress tests therefore
separately require tight harmonic convergence and verify that adaptive evaluation
does not enlarge the high-count evaluator's independent discrepancy. Their
independent envelope is 1e-3. This is an inherited special-function/conditioning
limit, not an achieved Bessel truncation error or a tolerance change to old tests.
The Z function is not modified in this task.

Root tests compare fixed 128, high-count fixed, and adaptive tolerances 1e-6,
1e-10 and 1e-13, requiring actual secant convergence before the iteration cap
and normalized determinant residual below 1e-8. Cases include warm and cold
references, stationary and drifting species, and growing/damped modes. The
aggregate complex-root difference in the common comparison matrix, including
fixed 128, was **3.72e-9**. The benchmark adaptive roots have real-component
errors at most **2.94e-11**, imaginary-component errors at most **2.96e-12**, and
normalized residuals below **2e-18**.

At lambda_p=1250 (`kperp=50`, `kpar=.01`, `v_A/c=.02`, theta_p=1,
theta_e=.2, mu_e=1836), the first ion-Bernstein branch gives:

| Evaluator | Complex frequency |
| --- | --- |
| Independent 1024-harmonic Stix calculation | 1.010537151312035 - 0.005693196262751 i |
| Adaptive at tol=1e-13 | 1.010537151340713 - 0.005693196260767 i |
| Historical fixed 128 | 1.010612074375771 - 0.005619207676129 i |

Both production runs converge, but the fixed-128 answer is under-resolved.
The default tolerance uses N_p=263 and N_e=7. Tightening the tolerance and the
1024-harmonic fixed reference confirm the adaptive branch.

## Performance measurements

`python3 tests/benchmark_bessel.py` reproduces the benchmark. It runs 100 supplied-
guess root searches per sample and reports the median of three samples. CPU times
cover the solver; wall times also include process launch/input/output. GNU Fortran
11.4.0 production flags were used. No timing threshold is asserted in tests.

| Case | Adaptive species N | Fixed128 CPU s | Adaptive CPU s | Adaptive wall s | Speedup vs fixed128 | Cache-only speedup |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Small arguments | 4, 2 | 0.321 | 0.0387 | 0.0416 | 8.28x | 1.04x |
| Moderate arguments | 20, 3 | 0.351 | 0.0828 | 0.0854 | 4.24x | 1.03x |
| Large arguments | 263, 7 | 0.360 | 0.685 | 0.687 | 0.526x | 1.03x |
| Three warm species | 4, 2, 3 | 0.557 | 0.0590 | 0.0673 | 9.43x | 1.27x |

The large-argument adaptive result is about 1.9x slower than the inaccurate
fixed-128 result, but about 7.4x faster than the converged fixed-1024 calculation.
Fixed-1024 CPU times for the other cases were 4.56, 5.04 and 9.61 seconds.
Compiler flags and root iteration counts also influence timings: the three-
species fixed128 search took 36 iterations versus 45 for the accurate evaluator.

Each adaptive species evaluated its Bessel array once per process. Per-species
hits/misses were 3399/1, 3799/1, 3799/1 and 4799/1 respectively: hit rates exceed
99.97%. The cache-only comparison invalidates entries before every trial
frequency while retaining allocation capacity; it uses the same adaptive N and
physics. Thus its speedup isolates coefficient reuse, while the larger overall
improvement mainly comes from reducing Z evaluations and harmonic sums. Complete
wall/CPU timings, root errors, residuals, counts and counters are in the JSON
artifact. Timing variation should be expected on other machines.

## Verification commands and files

Final verification:

```sh
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

**92 tests passed** (34.0 seconds in the recorded final full run). Checked probes
use `-O0 -g -fdefault-real-8 -fcheck=all -ffpe-trap=invalid,zero,overflow`.
The pinned legacy records remain exact. The example exits successfully with the
expected branch. Python compilation and whitespace checks are also performed.
FORD is not installed locally; the report is registered in the existing FORD/CI
workflow, but rendered documentation was not built locally.

Files added in this stage:

- `src/harmonic_selection.f90`: independent positive-tail selector.
- `src/bessel_cache.f90`: dynamic vector coefficients, cache, counters and cap history.
- `src/adaptive_kinetic.f90`: shared-equation evaluator, physical argument, controls and metadata.
- `src/kinetic_collisional_species.inc`: shared existing collisional species equations.
- `tests/harmonic_selection_probe.f90`, `bessel_cache_probe.f90`,
  `adaptive_core_probe.f90`, `harmonic_probe.f90`: isolated and production probes.
- `tests/test_harmonic_selection.py`, `test_bessel_cache.py`,
  `test_adaptive_dispersion.py`: 19 additional tests.
- `tests/benchmark_bessel.py`, `ADAPTIVE_BESSEL_BENCHMARK.json`, and this report.
- `inputs/example/adaptive_bessel.in`: working large-argument convergence example.

Existing files changed:

- `src/vars.f90`, `functions.f90`, `plume.f90`: input controls, validation and metadata finalization.
- `src/disprels.f90`, `cold_reference.f90`: narrow dispatch and shared-source integration.
- `Makefile`: module ordering, include dependencies and IEEE compilation flags.
- `linfpclib/linfpc.py`: optional controls, validation and input round trips.
- `tests/auxiliary_kernel.py`, `auxiliary_probe.f90`, `kinetic_core_probe.f90`:
  explicitly preserve the original fixed mathematical prototype.
- `tests/test_normalization.py`, `test_mixed_dispersion.py`,
  `test_auxiliary_normalization.py`, `test_kinetic_core.py`,
  `test_cold_reference_mixed.py`, `test_fully_cold_dispersion.py`: new build
  dependencies, explicit historical modes, and exclusion of nonnumeric metadata
  from numeric FPC parsing. Original numerical assertions are unchanged.
- `README.md`, `input.md`, `output.md`, `ford_project.md`,
  `.github/workflows/doc.yml`, `.gitignore`: documentation and example registration.

The pre-existing untracked `kinetic_species.inc`, `kinetic_harmonics.inc`, old
examples and license file were preserved. Unrelated user files, including
`plotter/example/cycle.sh`, were not changed by this work.

## Remaining numerical limitations

Bessel-tail convergence does not guarantee convergence near a resonance or in a
strongly cancelling tensor component. The Z accuracy limit above remains.
Tolerances below floating-point precision can select more harmonics but cannot
create more coefficient accuracy. The conservative bounds may retain more orders
than exact tails would require. Coefficient-generation temporary storage grows
with the resolved distribution width, which can exceed the requested susceptibility
cap; extreme unrepresentable recurrence orders return a numerical failure.

The explicit compatibility path keeps legacy cache behavior and lacks detailed
new hit/miss instrumentation. Intel rules are provided but only GNU was available
for validation. Existing exact-perpendicular warm restrictions, cold diagnostic
and collision restrictions, and cold-containing scan restrictions remain intact.
No root algorithm, cold susceptibility, plasma dispersion function, normalization,
or additional plasma physics was changed. No later stage is implemented.
