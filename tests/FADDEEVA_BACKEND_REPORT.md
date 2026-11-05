# Faddeeva backend integration report

## Milestone 1 checkpoint (completed before production integration)

The bundled source and header are byte-for-byte OpenSpecFun mirror files at
commit e98d4f427a20af69cb6f592ff65634f00e97b638. Author Steven G. Johnson,
Massachusetts Institute of Technology; upstream MIT license is retained and
compatible with PLUME's BSD-2-Clause license. The original MIT host could not
be reached. See third_party/faddeeva/README.md for retrieval date, checksums,
source URLs and the inactive mirror OpenLibm include hooks.

Normal make compiled the C99 source and ISO_C_BINDING module without new
runtime dependencies. No production backend choice/default had changed at this
checkpoint. All 110 pre-existing tests passed before changes.

Five standalone tests passed: all 32,783 audit arguments, signed-k convention,
frozen arbitrary-precision values, invalid/range/failure recovery, checked GNU
builds, complex ABI/kinds, and a fast-linked subnormal-result check. Maximum
scale-aware error versus SciPy was 2.132e-16; versus selected rounded 90-digit
references it was 5.001e-16. The native (no Python callback) susceptibility stress
intervention had maximum scaled error 6.113e-14.

A real build-specific issue was detected and fixed before this checkpoint:
GNU -ffast-math link startup enables FTZ and DAZ. Fortran gradual-underflow
control alone left DAZ active and incorrectly zeroed a representable result at
zeta=1e308. The small C ABI bridge now saves/restores fenv state and evaluates
normalization/sign operations in the default IEEE environment with the caller's
rounding mode. This also safely records failures with caller traps enabled.
No upstream numerical source was modified.

Native-code performance comparisons were completed before integration and saved
in FADDEEVA_STANDALONE_BENCHMARK.json. Both scalar loops and susceptibility
loops execute entirely in compiled code; Python callback timings are not used.
The production default may now be integrated (Milestone 2).

## Production integration and precision

`z_backend='faddeeva'` is now the default. `z_backend='legacy'` retains the
historical algorithm. Both operate independently of `bessel_mode='adaptive'`,
`'fixed'` and the added `'fixed128'` alias. Fixed128 requires `bessel_nfixed=128`;
other explicit fixed counts remain supported. No physical input or susceptibility
expression changed. The public interface remains `disprels::zet_in(kpar,zin)`:

\[
Z(\zeta)=i\sqrt{\pi}\,w(\zeta),\qquad
\operatorname{zet\_in}(k,\zeta)=
\begin{cases}Z(\zeta),&k>0,\\-Z(-\zeta),&k<0.\end{cases}
\]

The continuation is the entire Landau/Faddeeva continuation, rather than
complex conjugation of an upper-half-plane result. There is no extra branch
prescription or tolerance input in the new backend. Upstream `relerr=0` requests
machine accuracy. Both charge signs, propagation directions and imaginary
frequency half-planes retain their established argument definitions. The
upstream algorithm combines continued fractions at sufficiently large arguments
with Zaghloul–Ali Algorithm 916 methods elsewhere and local Taylor expansions
to avoid cancellation. These numerical branches are upstream code, not new
PLUME approximations.

The C99 implementation is compiled statically into PLUME with `cc -std=c99
-x c -O3 -fno-fast-math`; the C bridge and Fortran wrappers use strict IEEE
semantics. No CMake, C++ runtime, libcerf, SciPy/Python runtime dependency or
package installation is introduced. `make clean` also removes the new root-level
objects, and a clean parallel build was validated. `ldd plume.e` lists only
ordinary Fortran/C runtime libraries and libm (including the existing GNU
quadmath runtime), with no new third-party shared library.

The interface uses `real(c_double)` and `complex(c_double_complex)` with the
C99 `double complex` return and by-value argument ABI. Both this complex kind
and PLUME's compiler-promoted default complex have **128-bit storage** under
`-fdefault-real-8`. The actual ABI is exercised, rather than assumed.
`nrtype::DP` is deliberately not used: under these flags it is wider than the
production real kind. Unsupported kind choices produce compile-time type
mismatches rather than reinterpreted complex storage.

## Legacy compiler compatibility

The historical body is moved to `src/plasma_z_legacy.f90` as `zetout_legacy`.
Its numerical statements, constants, stopping tests, `error=1.e-7`, and early
unassigned-result returns are unchanged. Only its name and pure/intent
interface declarations change. This option intentionally retains historical
accuracy and undefined-result limitations; it is not an alternate robust mode.

A dispatcher that mutated a global Z status/counter made the GNU optimizer
assume plasma globals could change across every Z call. That changed kinetic
arithmetic by a few ulps and failed the exact regression. The final value path
is a **pure**, separately compiled value interface, with assigned quiet NaNs as
its failure result. It does not mutate counters or plasma globals. The explicit
standalone status interface remains available. The C bridge restores its
caller's floating-point environment, so this pure value contract also preserves
exception flags and modes. Dispersion guards classify invalid results after the
evaluation. Test-only linker instrumentation counts calls without adding a
production counter. The small public `disp` dispatch and the kinetic equations
remain unchanged; endpoint guards surround the calculation.

The pinned pre-normalization comparison at commit
`c6c7616e1fc58eb6118e1f9423a3f0df6d52ef20` passes **exact numeric equality** for
both original input fixtures with explicit `z_backend='legacy'` and
`bessel_mode='fixed128'`. No golden output or original assertion was relaxed.
Existing characterization/prototype fixtures explicitly select legacy Z, and
the independent J1.B.1 defects remain reproducible. Probe build lists add the
new dependencies. The previous collisional zz correction is byte-for-byte
preserved; choosing legacy Z does not undo that separate physics correction.

## Numerical validation

The final replay uses the original **32,783-input** audit grid, including
32,573 recovered signed stress arguments. Deduplication retains the original
`(kpar,zeta)` keys before applying the sign convention, so physically identical
effective arguments from different signs remain present as in J1.B.1.
The upper/lower imaginary-axis, weak growth/damping, real and nearly real
axes, old radius-3 and imaginary-1 transitions, old continued-fraction failure
region, large arguments and near-Z-zero cases are included.

The metrics are absolute error, relative error where the reference is nonzero,
and `abs(actual-reference)/max(1,abs(reference))`. Near a zero the last metric
is meaningful even when relative error is large. SciPy wofz is a related
algorithm, so the 11 independently generated frozen 90-digit mpmath values in
the unchanged audit artifact are also checked. Maximum scaled discrepancy is
**2.132e-16** versus SciPy and **5.001e-16** versus those rounded high-precision
values. Acceptance is 1e-13 for ordinary representable Z values.

| Argument regime | Inputs | Maximum scale-aware Z error | Genuine range failures |
| --- | ---: | ---: | ---: |
| asymptotic | 44 | 1.4302e-16 | 0 |
| nearly real | 33 | 1.2671e-16 | 0 |
| continued fraction | 12 | 1.2665e-16 | 0 |
| small series | 10 | 1.4021e-16 | 0 |
| radius 3 boundary | 50 | 2.1311e-16 | 0 |
| nearly real asymptotic | 12 | 5.5511e-17 | 0 |
| imaginary 1 boundary | 12 | 2.0008e-16 | 0 |
| large arguments | 24 | 1.3878e-17 | 0 |
| domain and range | 9 | 1.4711e-16 | 2 |
| near zero | 4 | 3.3881e-21 | 0 |
| recovered stress | 32573 | 1.3536e-16 | 0 |

Range cases `-27i` and `-30i` genuinely exceed binary64 range. `-14i` and
`-26i` are accepted. NaN/Inf inputs and nonfinite k are rejected, and a
representable subnormal `Z(1e308)` is retained even in a fast-math-linked
executable. The fast-link test reads returned hexadecimal storage bits to
avoid formatted-I/O DAZ effects. These are binary64 component-range checks,
not clipping. Checked builds use `-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`. Successful calls following failures are
explicitly tested in both standalone and production paths.

## Susceptibility stress and remaining conditioning

The original primary alpha-particle state is
`tau_d=2, mu=.25, alpha=1.3, Q=.5, D=.05, u=.4`, with
`w=.73+.13i, kpar_d=-.31, kperp_d=3, vA/c=.02`.
The production reference state is warm protons/electrons/alphas from
`stress_systems()[1]` (`beta_ref=.4`, `alpha_ref=1.7`).

| Evaluation at the primary state | Scaled alpha susceptibility discrepancy |
| --- | ---: |
| Historical Z, high fixed count | approximately 7.624e-4 |
| Faddeeva, fixed 1024, SciPy Bessel/Stix reference | 6.113e-14 |
| Faddeeva, adaptive default, SciPy Bessel/Stix reference | 1.924e-12 |
| Faddeeva, fixed 1024, identical Bessel coefficients in Stix | 3.195e-14 |
| Faddeeva, adaptive, identical Bessel coefficients in Stix | 3.158e-14 |

Thus replacing Z resolves the audited asymptotic-transition defect to the
accuracy established by J1.B.1. Tests require <2e-13 for the fixed-1024 alpha
response and for the matched-coefficient alpha comparison in either mode.
The native standalone fixed-coefficient intervention also gives 6.113e-14;
it uses compiled Fortran/C without a Python Z callback.

Every complex tensor element is compared across both signs of kpar and both
frequency half-planes, kinetic-only, warm-reference mixed, cold-reference mixed,
both normalizations and independent Z/harmonic choices. **The new Z evaluator
does not guarantee machine-precision susceptibility**: the whole three-species
stress tensor has maximum error 1.527e-10 (electron response), although the
original problematic alpha tensor is accurate to 6.113e-14.

The broader new production grid has maximum raw SciPy-Bessel/Stix discrepancy
**3.619e-6**, in fixed128 cold-reference warm-electron/alpha calculations at
`kpar_d=-.31, kperp_d=3, w=.73-.08i, vA/c=.02`. This is not concealed by
requiring a loose reference tolerance: raw discrepancies are measured and
recorded separately. Replacing the independent reference's Gamma coefficients
with the **actual unchanged production coefficients**, while retaining an
independent Stix formulation and SciPy Z, reduces the maximum across that grid
to **1.663e-10**. The matched-coefficient assertion uses `rtol=2e-10,
atol=2e-9`; these are conditioning-aware tensor tests, distinct from the
1e-13 direct-Z criterion. Fixed128 retains the historical Bessel approximations;
the vector/adaptive Bessel implementation and cancellation are likewise
unchanged. Improving those approximations or summations is outside this task.
No existing accuracy tolerances were weakened to accommodate this observation.

Weakly collisional tests use the already-corrected zz response and independent
Gaussian moments with SciPy Z, holding production Gamma coefficients identical.
The rate is PLUME's existing `nu_s=sqrt(mu_s*tau_d)/(sqrt(2)*Kn)`; its definition
is unchanged. Fixed128 and adaptive checked/production builds pass the moment
comparison (`rtol=1e-10, atol=1e-9`). Historical correction tests also remain
intact. No collisional derivation or tensor equation was revisited.

## Independent roots and dispersion residuals

SciPy solves the independent finite-mass Stix/cold wave matrix at identical
physical inputs (`vA/c=.01`, normally `kpar_d=.06, kperp_d=.05`; beam uses
`.3,.1`, parallel uses `kperp_d=0`). PLUME uses its unchanged secant update and
convergence tolerances. Both production and checked builds find the same
branches. Differences below are complex-frequency absolute differences; full
real/imaginary components, legacy roots, parameters and normalized matrix
residuals are in `FADDEEVA_BACKEND_RESULTS.json` and the test case definitions.

| Case | Accurate PLUME root | Independent root | Absolute difference | `abs(det(Lambda/max(abs(Lambda))))` |
| --- | --- | --- | ---: | ---: |
| kinetic | 0.060987152361867-1.2709968e-05i | 0.060987152362065-1.2709968e-05i | 1.975e-13 | 6.582e-20 |
| warm_reference_mixed | 0.061538419583637-5.5643255e-06i | 0.061538419583664-5.5643255e-06i | 2.692e-14 | 4.733e-24 |
| cold_reference_mixed | 0.059167349121515-1.5887687e-06i | 0.059167349122274-1.5887688e-06i | 7.650e-13 | 1.821e-19 |
| parallel_cold_reference | 0.057642296896783+1.7002145e-13i | 0.057642296898146+6.0312797e-13i | 1.430e-12 | 1.981e-19 |
| growing_beam | 0.30574584309027+0.079649334i | 0.3057458430903+0.079649334i | 3.855e-14 | 1.643e-22 |

Maximum checked/production root difference from the independent solution is
1.431e-12. The growing alpha-beam root has imaginary part +.0796493343;
Landau-damped oblique roots have negative imaginary parts. A small determinant
residual alone is not asserted to establish root accuracy: branch agreement and
separate real/imaginary frequencies are also checked. Accurate and historical
roots need not match exactly; both are recorded, with exact compatibility
verified separately in the pinned legacy regression.

## Failure handling

| Situation | Standalone status | Accurate production behavior |
| --- | ---: | --- |
| Finite representable value | 0 | Ordinary susceptibility |
| NaN/Inf argument or k | 1 | Assigned NaN, invalid response |
| `kpar=0` Z call | 2 | Assigned NaN, invalid response |
| Output outside binary64 component range | 3 | Assigned NaN, numerical failure |
| Backend/fenv/nonfinite numerical failure | 4 | Assigned NaN, numerical failure |

The status-returning API is `calc_z_faddeeva[_signed]`. The unchanged value-only
`zet_in` interface encodes failure as assigned quiet NaNs without mutating the
physical state. Evaluation guards restore caller traps/flags and check returned
determinants and arithmetic exceptions. Existing dispersion status meanings are
retained: invalid frequency 1, zero frequency 2, cold poles 2–4, nonfinite
numerical response 5. Range and backend failure both map to 5 in production;
they remain separately inspectable with the standalone API. No failed Z value
is replaced by zero or a legacy approximation.

All-kinetic, warm-reference and cold-reference failed trials followed by valid
trials recover cleanly, including checked builds. Secant initial/trial nonfinite
values and iteration exhaustion return negative failure flags for accurate
runs. Root refinement refuses those results. Map searches mask isolated invalid
points, exclude failed refinements and continue with valid points; an all-warm
map containing omega=0 still recovers the independently checked root. Existing
mixed warning text and legacy failure behavior are preserved. Invalid direct
single-frequency runs skip supplemental eigenfunction/heating diagnostics.
Fully cold tests instrument the C value entry and verify **zero Z calls**.

## Native performance

All timings below are CPU seconds, median of three samples on the toolchain
recorded in the JSON artifacts. They are representative measurements, not
portable performance guarantees. No Python callback occurs in a timed Z loop.
The bridge's fenv save/default/restore cost is intentional: it handles traps
and the fast-math FTZ/DAZ behavior safely.

### Standalone checkpoint, before changing the default

600,000 scalar calls per distribution (native loops over small, radius-3/
continued-fraction and large arguments):

| Distribution | Legacy seconds | Faddeeva seconds | Ratio |
| --- | ---: | ---: | ---: |
| small | 0.040132 | 0.196784 | 4.903 |
| transition | 0.039712 | 0.202319 | 5.095 |
| large | 0.024901 | 0.204638 | 8.218 |

Native single-species susceptibility loops, 200 calls with identical Gamma and N:

| Case | N | Legacy seconds | Faddeeva seconds | Ratio |
| --- | ---: | ---: | ---: | ---: |
| small | 8 | 0.001612 | 0.002478 | 1.537 |
| stress | 128 | 0.020048 | 0.032641 | 1.628 |
| large | 512 | 0.076535 | 0.126370 | 1.651 |

### Final production endpoint

250 cached full-determinant evaluations, fixed counts held equal between
backends, no call-count instrumentation, no Python callbacks, and the production
fast-math link flags (strict C/wrapper objects). This includes guard overhead.

| Case | N | Legacy seconds | Faddeeva seconds | Ratio |
| --- | ---: | ---: | ---: | ---: |
| small | 128 | 0.033513 | 0.053204 | 1.588 |
| moderate | 128 | 0.030655 | 0.059006 | 1.925 |
| large | 512 | 0.168221 | 0.288027 | 1.712 |
| cold_reference | 128 | 0.219668 | 0.247955 | 1.129 |

Performance depends on harmonic count, argument distribution and time spent in
Bessel/susceptibility arithmetic. These fixed-count timings isolate the backend
choice; changing adaptive cutoffs can independently change total runtime.
The older audit's Python-callback timings are not used for this comparison.

## Reproduction and final regression

Toolchain: GNU Fortran and C 11.4.0, Linux x86_64; detailed platform, build flags
and tool versions are stored in the result artifacts. SciPy/NumPy are used only
for testing. The frozen mpmath references avoid requiring mpmath at runtime.

```bash
make clean
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v

# Reproduce final numerical tables and uninstrumented native CPU timings:
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/measure_faddeeva_backend.py \
  --output /tmp/FADDEEVA_BACKEND_RESULTS.json

# Reproduce the standalone checkpoint benchmark (keep historical artifact intact):
python3 tests/benchmark_faddeeva.py --output /tmp/FADDEEVA_STANDALONE_BENCHMARK.json
```

The clean build succeeds. **125 tests pass**, including all **110 tests present
at task start**, 5 new standalone tests and 10 new production tests. The original
request's earlier 92-test count predates the additional audit/zz-correction
coverage present in this working tree. No expected numerical data was replaced.
The pinned exact legacy comparison passes, as do checked bounds/FPE builds,
independent collisional tests, both examples and Python backend round-tripping.
`git diff --check` passes. All work remains local and uncommitted.

## Remaining limits

- Binary64 overflow of the analytic continuation is a real limitation, not a
  removable algorithm failure. Returned failure values are explicit.
- Accurate Z does not fix Bessel approximations, unresolved harmonic tails,
  kinetic cancellation, determinant conditioning or root branch selection.
- Exact historical reproduction is established for the pinned fixtures and GNU
  toolchain, not a cross-compiler bitwise guarantee. Other C/Fortran compilers
  were not available for execution tests; standard C99/ISO_C_BINDING is used.
- Explicit legacy mode intentionally retains uninitialized-result and accuracy
  defects. Use it for reproduction, with the documented safeguards applying to
  the accurate backend.
- Further Bessel/summation conditioning work should begin with the recorded
  matched-coefficient error budget, without changing this stage's physics or
  making a Z tolerance into a susceptibility-error guarantee.

## Exact file inventory for this stage

This inventory is relative to the working tree captured at task start, rather
than `git diff` against HEAD (which also contains earlier uncommitted stages).
Pre-existing unrelated modifications, the collisional include/correction report,
Z audit report/JSON, cold kernel, kinetic includes, harmonic selector and Bessel
cache are preserved. Their pre-task SHA256 values were checked.

Modified files:

- `.github/workflows/doc.yml`
- `.gitignore`
- `Makefile`
- `README.md`
- `ford_project.md`
- `input.md`
- `linfpclib/linfpc.py`
- `output.md`
- `src/adaptive_kinetic.f90`
- `src/disprels.f90`
- `src/functions.f90`
- `src/vars.f90`
- `tests/adaptive_core_probe.f90`
- `tests/auxiliary_kernel.py`
- `tests/auxiliary_probe.f90`
- `tests/benchmark_bessel.py`
- `tests/kinetic_core_probe.f90`
- `tests/plasma_dispersion_audit.py`
- `tests/plasma_dispersion_probe.f90`
- `tests/plasma_z_direct_probe.f90`
- `tests/test_adaptive_dispersion.py`
- `tests/test_auxiliary_normalization.py`
- `tests/test_cold_reference_mixed.py`
- `tests/test_collisional_zz.py`
- `tests/test_fully_cold_dispersion.py`
- `tests/test_kinetic_core.py`
- `tests/test_mixed_dispersion.py`
- `tests/test_normalization.py`
- `tests/test_plasma_dispersion_audit.py`

Added files:

- `src/plasma_z_backend.f90`
- `src/plasma_z_faddeeva.f90`
- `src/faddeeva_bridge.c`
- `src/plasma_z_legacy.f90`
- `tests/faddeeva_build.py`
- `tests/FADDEEVA_STANDALONE_BENCHMARK.json`
- `tests/faddeeva_count_probe.c`
- `tests/test_faddeeva_standalone.py`
- `tests/faddeeva_direct_probe.f90`
- `tests/FADDEEVA_BACKEND_REPORT.md`
- `tests/faddeeva_selected_probe.f90`
- `tests/faddeeva_harness.py`
- `tests/faddeeva_production_probe.f90`
- `tests/measure_faddeeva_backend.py`
- `tests/test_faddeeva_backend.py`
- `tests/benchmark_faddeeva.py`
- `tests/FADDEEVA_BACKEND_RESULTS.json`
- `tests/faddeeva_probe.f90`
- `third_party/faddeeva/LICENSE`
- `third_party/faddeeva/Faddeeva.h`
- `third_party/faddeeva/Faddeeva.cc`
- `third_party/faddeeva/README.md`
- `inputs/example/faddeeva_accuracy.in`
- `inputs/example/legacy_numerics.in`
