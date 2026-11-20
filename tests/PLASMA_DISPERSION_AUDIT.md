# Stage J1.B.1: independent plasma dispersion function and conditioning audit

Historical audit snapshot: its measured JSON and source fingerprints are retained.
The collisional zz and signed-harmonic inconsistencies identified below have
subsequently been corrected; see [the correction report](COLLISIONAL_ZZ_CORRECTION.md).
The Z-evaluation findings remain applicable and the Z backend is unchanged.

## Findings and audit boundary

**The existing Z evaluator is the dominant source of the largest collisionless
J.1 stress discrepancy. Its asymptotic branch, rather than the scalar stopping
tolerance alone, is responsible.** Holding all Bessel coefficients, harmonic
counts, input parameters and susceptibility expressions fixed, substituting
SciPy Z values reduces the reproduced discrepancy from **7.62447135e-4** to
**5.19e-14**. Tightening the existing internal tolerance to 1e-15 leaves the
large discrepancy essentially unchanged.

Additional findings:

- Maximum valid-domain scale-aware Z error on this audit grid: **2.268e-4**,
  immediately outside `|zeta|=3`. A near-real-axis continuity defect accompanies
  the asymptotic approximation there.
- Tightening the tolerance improves the small-argument series. The continued
  fraction encounters its unchanged 30-iteration cap and can remain inaccurate
  at **2.61e-9** despite much tighter requested tolerances.
- Some lower-half-plane inputs return an **unassigned function result** rather
  than an explicit failure. Several of their true values are representable.
- The independent NHDS/Stix reference is confirmed by arbitrary precision at the
  largest stress point. A smaller residual, about 1.9e-10 at another point, does
  involve binary64 conditioning in both implementations.
- A separate anisotropic **zero-collision-rate inconsistency** exists between
  the collisionless and collisional `zz` formulas. It is traced below to different
  parentheses and is not a Z error or the cause of the collisionless J.1 result.

This stage adds only audit files under `tests/`. No production source, equations,
normalizations, cache semantics, defaults, root method, or pre-existing files
were edited. No commit or push was made. Stage J1.B.2 is not implemented.

The machine-readable record is
[PLASMA_DISPERSION_AUDIT_RESULTS.json](PLASMA_DISPERSION_AUDIT_RESULTS.json).
Its source hashes identify the actual working-tree implementation, including
pre-existing uncommitted J.1 code. Timing samples are measurements, not CI limits.

## 1. Definition, conventions and source call paths

For Im(zeta)>0 the mathematical definition is

\[
 Z(\zeta)=\frac1{\sqrt\pi}\int_{-\infty}^{\infty}
           \frac{e^{-t^2}}{t-\zeta}\,dt
         =i\sqrt\pi\,w(\zeta),\qquad
 w(z)=e^{-z^2}\operatorname{erfc}(-iz).
\]

The last expression gives the entire Landau analytic continuation into the lower
half-plane. There is no real-axis singularity of the continued Z function.
PLUME uses `exp(i k.r - i omega t)`, with positive Im(omega) representing growth.

`disprels::zet_in(kpar,zeta)` implements

\[
 Z_k(\zeta)=\begin{cases}Z(\zeta),&k_\parallel>0,\\
                         -Z(-\zeta),&k_\parallel<0.\end{cases}
\]

The independent oracle must therefore use
`sign(kpar)*i*sqrt(pi)*wofz(sign(kpar)*zeta)`; comparing negative-k inputs directly
with `i*sqrt(pi)*wofz(zeta)` would be incorrect. Positive thermal normalization
factors preserve the sign of kpar. The kinetic argument is physically
`(omega-kpar_d*u_s-n*mu_s/Q_s)/(kpar_d*sqrt(mu_s*tau_s_d))`.
The collisional expression adds `+i*nu_ns` to its numerator.

Useful identities for the underlying entire function are

\[
 Z(-z)+Z(z)=2i\sqrt\pi e^{-z^2},\quad
 Z(-z^*)=-Z(z)^*,\quad
 Z(z^*)=Z(z)^*+2i\sqrt\pi e^{-(z^*)^2}.
\]

Ordinary conjugation alone is not the Landau continuation. The audit tests both
signs of kpar and the above identities. Passing a symmetry test does not establish
absolute accuracy: the inaccurate asymptotic values can still satisfy reflection.

| Source/path | Z implementation used |
| --- | --- |
| `src/disprels.f90`, `zet_in` / private `zetout` (around lines 4724–4844) | Only active production Z implementation. Wrapper chooses the signed continuation. |
| `disp0`, `kinetic_species.inc` | Every collisionless harmonic calls `zet_in(kpar,tsi(n))`. |
| `disp1`, `kinetic_collisional_species.inc` | Same Z implementation, evaluated at collision-shifted arguments. |
| `kinetic_core`, `adaptive_kinetic` | Production passes the same `disprels::zet_in` callback. Includes share the existing species equations. |
| Warm-reference `rho_ref` / `d_ref` calculations | Adapter changes thermal coordinate state, not the Z backend. |
| Cold-reference mixtures, `cold_reference::calc_warm_auxiliary` | Auxiliary normalization changes coordinates; the same callback is supplied for warm components. |
| `disp_mixed_parallel` | Same `zet_in` for the analytical parallel warm response. Fully cold components use no Z. |
| Root refinement, maps, eigenfunction/heating derivatives | Evaluate `disp`, and hence the same backend; there is no separate root-specific Z algorithm. |
| `fpc::wparth_from_ratio` | Imports `zet_in` but does not call it. Other FPC velocity-space expressions are not additional Z implementations. |
| `src/gk_disp_gf.f90` | Contains calls to external `gk_funcs::zetvec`, but is excluded from the default Makefile. Neither `gk_funcs` nor its Z implementation is supplied in this tree. It is not an active PLUME executable path and is not characterized here. |

The mathematical reference definitions are documented by
[SciPy `wofz`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wofz.html).
The independent arbitrary-precision check evaluates the complex complementary
error function using [mpmath](https://mpmath.org/doc/current/functions/expintegrals.html#erfc).
No external implementation code was copied into PLUME.

## 2. Existing algorithm and numerical contracts

The exact source is compiled unchanged in a temporary build. An additional,
test-only extracted copy permits changing `error` and recording branches/stops.

Before choosing a numerical branch, lower-half-plane and real-axis arguments
set `hold=i*sqrt(pi)*exp(-zeta**2)` when
`abs(Re(zeta**2))<174` and `abs(Im(zeta**2))<50000`.
Outside this region, positive Re(zeta**2) sets `hold=0`; otherwise the function
sets a local unused `fail` variable and returns **without assigning its result**.
`zet_in` likewise has no assigned result for kpar=0. Warm exact-perpendicular
propagation remains unsupported; fully cold calculations bypass Z altogether.

The numerical branches are:

1. **Series:** `|zeta|<=3` and `abs(Im(zeta))<1`. Iteration stops when changes in
   both components of an auxiliary running sum fall below `error`. This is an
   absolute change criterion, not a guaranteed final relative Z error.
2. **Asymptotic derivative expansion:** `|zeta|>3`. The routine sums an auxiliary
   third-derivative series and reconstructs lower derivatives and Z algebraically.
   It stops either on an absolute change in that auxiliary sum or when the
   one-norm `abs(Re(term/previous))+abs(Im(term/previous))` exceeds one. The latter
   stop is independent of the requested tolerance and can occur before terms
   increase in complex modulus. Exponential continuation terms use different
   factors for Im(zeta)>0, =0 and <0.
3. **Continued fraction:** remaining points inside radius three, with
   `abs(Im(zeta))>=1`. A componentwise relative iterate-change test is used,
   but iteration stops after 30 steps even if that test is unsatisfied. The
   ordinary API returns no convergence status.

The historical `error=1.e-7` is therefore neither a uniform absolute nor a
uniform relative accuracy guarantee for Z.

### Continuity and failure examples

At the nearly real point x=3.01:

| zeta | Historical Im Z | Reference Im Z |
| --- | ---: | ---: |
| 3.01 - 1e-12 i | 4.11958526774e-4 | 2.05979263318e-4 |
| 3.01 | 2.05979263456e-4 | 2.05979263456e-4 |
| 3.01 + 1e-12 i | 1.38710e-13 | 2.05979263594e-4 |

The real component also has about 5.286e-5 error. The continued function is smooth;
the asymptotic prescription creates these jumps. Tightening `error` does not
supply the missing uniformly accurate behavior near this transition.

At zeta=-14i, the true value is approximately `4.691626e85 i`, representable in
binary64, yet the source takes the unassigned-return path. The raw production
return can be finite register-dependent garbage; a finite-value check alone is
insufficient. The instrumented copy marks this path with a NaN and an audit status,
without repairing production. Other tested examples include -26i and 200-200i.
At -27i and -30i the true ordinary Z value exceeds binary64 range; those are genuine
representation failures and are not counted as finite-reference accuracy errors.

## 3. Harness, references and test domain

`plasma_dispersion_audit.py` builds temporary shared libraries with the unchanged
source modules. The Fortran C-ABI probe has three diagnostic choices:

- **0:** actual `disprels::zet_in` from the unchanged module;
- **1:** instrumented extracted algorithm, with explicitly supplied tolerance;
- **2:** reference callback evaluating SciPy at the **exact complex arguments
  computed by the Fortran susceptibility expressions**.

For susceptibility comparisons, the probe includes `kinetic_species.inc` and
`kinetic_collisional_species.inc` unchanged. Gamma_n arrays are supplied explicitly
and identically to each Z choice. Local arrays avoid the production Bessel cache.
This intervention changes only the Z values, not their arguments, Bessel values,
harmonic cutoff, physical parameters or algebra. The independent Stix expression
is implemented separately and its affine-in-Z harmonic terms support an error
budget. The existing reference in `mixed_reference.py` is not edited.

A direct executable additionally links the already-built production objects,
independent of the root finder. Its results agree with the shared-library production
entry within `2e-14` absolute/relative on the comparison set. Checked builds use
bounds checks and floating-point traps; audit wrappers save/restore exception
flags and halting modes so failures can be recorded rather than aborting the process.

The recorded grid contains **32,783 Z inputs**, including **32,573 unique signed
argument records recovered from the J.1 stress grid**. Recovery covers all
`n=-1024..1024` calls, rather than just low-order harmonics. Additional points cover:

- Positive/negative real arguments; imaginary parts above, on and below the axis.
- Nearly real arguments, zeta=0, and small resonance arguments.
- Radii immediately inside/outside three and imaginary parts around +/-1.
- Large magnitudes through 1e6, lower-half-plane range limits and large phases.
- Points close to a numerically located Z zero near
  `1.9914668428338795 - 1.3548101281120064 i`.

The stress reconstruction covers five J.1 plasma layouts, both parallel signs,
`kperp_d=.07,3,30`, and frequencies `.73+.13i` and `.73-.08i`: 60 plasma
configurations and 108 warm-species evaluations. This includes anisotropy,
drifts, protons/electrons/alphas, warm and cold references, and both normalizations'
thermal coordinate choices. Auxiliary scales `.01,.1,1,10` are checked separately.

References use SciPy `wofz` / `ive`. Selected Z values are also checked with
**90 decimal digits** using mpmath; full susceptibility references use **80 digits**.
A 128-harmonic high-precision sum resolves the primary stress point, while 16
harmonics suffice for the small-lambda electron precision control. Optional
arguments allow promoting the very same binary64 Gamma and Z values to high
precision, isolating arithmetic from special-function and argument-rounding errors.

Metrics are `abs(actual-reference)`, ordinary relative error where the reference
is useful, and the scale-aware value
`abs(actual-reference)/max(1,abs(reference))`. Tensor maxima are componentwise.
Relative error is intrinsically unhelpful very close to a Z zero. Nonfinite or
unassigned-domain returns are counted separately, not included in finite error
maxima. Large absolute errors on exponentially huge Z values must likewise be
interpreted with their relative/scaled errors.

## 4. Z accuracy by regime and internal tolerance

Entries below are maximum **scale-aware** errors; the first column calls actual
production. Other columns change only the instrumented algorithm's tolerance.

| Argument regime | Historical 1e-7 | Copy 1e-10 | Copy 1e-13 | Copy 1e-15 |
| --- | ---: | ---: | ---: | ---: |
| Small series interior | 1.051e-8 | 7.136e-12 | 1.824e-15 | 1.398e-15 |
| Continued fraction | 9.119e-8 | 2.607e-9 | 2.607e-9 | 2.607e-9 |
| General asymptotic sample | 1.855e-4 | 1.855e-4 | 1.855e-4 | 1.855e-4 |
| Radius-three transition | 2.268e-4 | 2.268e-4 | 2.268e-4 | 2.268e-4 |
| Nearly real, outside radius three | 2.268e-4 | 2.268e-4 | 2.268e-4 | 2.268e-4 |
| Large arguments away from transition | 1.735e-18 | <=2.2e-19 | <=2.2e-19 | <=2.2e-19 |
| Near a Z zero | 1.057e-8 | 3.334e-11 | 7.261e-14 | 7.261e-14 |
| Recovered J.1 stress arguments | 1.040e-4 | 1.040e-4 | 1.040e-4 | 1.040e-4 |

The continued fraction hits its cap at ten grid points in its dedicated regime
for each tighter tolerance; further cap events occur in transition/near-zero and
sampled-argument groups. The exact counts, absolute/relative errors, exception flags
and worst arguments are in the JSON record. The largest ordinary relative error
on the radius-three transition is about **6.361e-4**. Near a true zero, a relative
error of order 1e7 can accompany an absolute error of only 1e-8.

No invalid/divide-by-zero/overflow flags occurred in the 108 valid stress tensor
evaluations. Underflow was recorded in 96, associated with negligible high-order
terms. The separate domain/range group records five undefined or nonfinite cases.
Some finite exponentially large values are very accurate: at -13i the absolute
error is about 2.51e58 but the relative error is only 2.85e-16.

SciPy versus the 90-digit Z reference has maximum scale-aware discrepancy
**4.45e-16** on the selected independent checks, including the primary problematic
arguments and the representable -14i case.

## 5. Reproducing and isolating the J.1 susceptibility discrepancy

The largest case is the drifting warm alpha component of the following system,
with tuples `(tau_d,mu,alpha,Q,D,u)`:

```
protons:   (.4, 1,    1.7,  1,  1,    0)
electrons: (.1, 1836,  .8, -1,  1.1, -.2)
alpha:     (2,  .25, 1.3,  .5,  .05,  .4)
omega=.73+.13i, kpar_d=-.31, kperp_d=3, v_A/c=.02
```

The alpha argument is lambda=11.7. The exact positive weighted tail beyond N=128
is approximately **5.94e-121**, so 128 already resolves the Bessel truncation here.
Frozen Gamma and N=1024 give the same physical conclusion.

The worst element is chi_zz. Its high-precision value is approximately
`17.68204462178117 + 148.41347438925067 i`.
The historical evaluator gives
`17.73727297 + 148.51315481 i`, an absolute discrepancy of **0.113958**.

The dominant harmonic is n=3:

```
Fortran zeta =  2.947038584816192 - 0.5930573003500077 i
Z effective input = -2.947038584816192 + 0.5930573003500077 i
```

It lies just outside radius three. The asymptotic expansion stops after five
terms because of its ratio criterion. The returned signed Z is
`-.34284053922727314 - .07985883510876199 i`, versus the reference
`-.34283323133103083 - .07975504344539389 i`.
Its absolute error is **1.04048617e-4** and relative error **2.95603e-4**.
The associated chi_zz error contribution is approximately
`0.05508547257 + 0.09925429944 i`, accounting for almost all of the total error.

| Candidate source | Controlled measurement | Conclusion |
| --- | --- | --- |
| Z evaluation | Frozen Gamma/N, historical vs SciPy substitution: scaled chi error 7.62447e-4 -> 5.19e-14 | Dominant source established by direct intervention. |
| Internal tolerance alone | 1e-7 -> 1e-15 gives 7.62447e-4 -> 7.62451e-4 | Insufficient; asymptotic ratio termination remains active. |
| Adaptive harmonic truncation | Actual adaptive vs frozen N=1024: 1.92e-12 maximum scaled component difference | Not responsible. |
| Higher fixed-harmonic evaluation | Actual fixed1024 vs frozen Gamma evaluator: 3.42e-15 | Same discrepancy survives adequate fixed summation. |
| Historical Bessel coefficients | Fixed128 vs accurate frozen Gamma/default Z: 1.14e-6 maximum scaled component difference | Secondary coefficient error; it does not explain the 7.62e-4 error. |
| Auxiliary normalization | b*=.01,.1,1,10 retains 7.62447e-4, varying by <4e-15 in the reported scaled maximum; substituted results <=5.2e-14 | No physical normalization error found. |
| Harmonic/algebraic conditioning | Harmonic condition 9.38; raw constant-plus-Z condition 194; relative Z condition 94.4 | Amplifies input Z error, but does not generate the large error with accurate Z. |
| Binary64 assembly at fixed Gamma/Z | Fortran vs promoted 80-digit arithmetic: 1.81e-14 scaled at the primary point | Roundoff is negligible compared with the original discrepancy. |
| Independent Stix formulation | Double-precision Stix vs 80-digit sum: 5.11e-14 scaled | Reference validated at the largest stress point. |
| Wave-matrix conditioning | Stress matrix cond2 about 528; chi discrepancy exists before matrix assembly | Cannot be the origin of the species susceptibility error. |

The susceptibility is affine in the individual Z_n values. The audit constructs
independent coefficients `chi=sum(C_n+B_n*Z_n)+drift`, then predicts
`delta_chi=sum(B_n*delta_Z_n)`. This closes the observed full tensor difference to
**9.34e-13 absolute**, providing an elementwise error budget rather than an
attribution based only on determinants or on the presence of a loose tolerance.

### Smaller precision floor and the reference calculation

Across all 108 species evaluations, the largest discrepancy remaining after Z
substitution is **1.914e-10**, in a small-lambda drifting-electron case:
`tau=.1,mu=1836,alpha=.8,Q=-1,D=1.1,u=-.2`, omega=.73-.08i,
kpar=-.31, kperp=.07, a=.02. An independent 80-digit calculation finds:

- Double Stix reference error: **1.324e-10** scaled.
- Fortran with SciPy Z error: **5.900e-11** scaled.
- Fortran versus promoted binary64 Gamma/Z arithmetic: **7.627e-11** scaled.
- Promoted Gamma/Z versus exact special-function calculation: **1.727e-11** scaled.

Thus the residual includes ordinary finite-precision cancellation in both
implementations and amplification of binary64 special-function/argument errors.
It should not be mislabeled as a residual large Z algorithm failure or evidence
that the independent reference is exact to machine precision at every component.
No susceptibility reformulation is justified by the primary 7.62e-4 discrepancy.

## 6. Collision controls and a distinct zero-rate formula issue

Using identical collisional equations and Gamma arrays while replacing only Z
still changes susceptibility at nonzero nu: maximum scaled historical-versus-
reference-Z differences are approximately **3.19e-4**, **3.49e-4**, and **1.36e-4**
for nu=1e-4, .01 and .1 at the alpha stress state. This isolates Z sensitivity
within the collisional implementation; it does not independently validate its
collision operator or full physical tensor.

A stronger control forces the **collisional** expression at exactly nu=0.
With accurate Z, its anisotropic chi_zz differs from the collisionless expression by

```
20926.47156817147 - 6396.507320744184 i
```

The scaled discrepancy is **146.406**, while the isotropic alpha=1 control agrees
within **5.11e-15** scaled. The original production collisional evaluator with
`collision_type=1,Kn=1e30` confirms the test probe within **3.45e-15** scaled.

The exact source difference is in the nonzero-harmonic zz polynomial. The
collisionless expression places `n^2*g^2*(1-alpha)` **inside** multiplication by
`n*g*(Z_n-Z_-n)`. The collisional expression places it **outside** that product.
At zero rate the resulting difference is

\[
 \Delta\chi_{zz}=N_{zz}C\sum_{n=1}^N\Gamma_n n^2g^2(1-\alpha)
                  [1-ng(Z_n-Z_{-n})],
\]

with `C=1/(kpar_d*sqrt(mu*tau_d))` and
`N_zz=2D/(a^2*w^2*Q^2*kpar_d^2*alpha*tau_d)`.
This predicts the measured complex difference within **4.55e-12 absolute**.
It vanishes for alpha=1. The same mismatch persists as nu tends to zero.

This is a separately identified formula-review issue, independent of the Z
backend and of the original collisionless J.1 stress discrepancy. It requires a
review of the intended weak-collision derivation before any correction. This
audit deliberately makes no algebraic fix.

## 7. Limited dispersion-matrix and root sensitivity

The ordinary dielectric and electromagnetic wave matrices use identical cold
responses and physical parameters; only warm-species Z values change. Roots in
the first three rows are obtained using the existing production secant method
and supplied guesses, with fixed1024 harmonics. Matrices are compared at the
**same resulting complex frequency**. Define
`eta=abs(det(Lambda))/max(abs(Lambda))**3`, using the same reference scale for both
values in each comparison.

| Case / historical frequency | Max scaled element difference | eta historical | eta reference Z | cond2 reference matrix |
| --- | ---: | ---: | ---: | ---: |
| Stable cold control: .059725410002175 + 9.13e-16 i | 0 | 2.61e-24 | 2.61e-24 | 5.28e17 |
| Damped cold-p/warm-e: .059170218376222 - 1.5890036172e-6 i | 1.21e-8 | 1.23e-18 | 6.84e-18 | 4.67e13 |
| Growing warm-alpha beam: .305745843090275 + .079649334334368 i | 5.92e-15 | 1.20e-22 | 1.20e-22 | 1.07e17 |
| Cancellation stress at .73+.13i (not a root) | 6.18e-6 | 1.66038894e-4 | 1.66038231e-4 | 5.28e2 |

Large condition numbers at a root are expected because the wave matrix is
singular there. They do not independently establish a large root error. Likewise,
a tiny normalized determinant alone does not establish the accuracy of omega.

As a local diagnostic, central differences of the reference-Z determinant at
steps 5e-7, 1e-6 and 2e-6 times max(1,abs(omega)) give a stable linear correction
`(det_historical-det_reference)/(d det_reference/d omega)`. The damped root's
estimate is approximately **1.267e-12 - 9.15e-17 i**. The growing root's estimate
is about **9.44e-17 - 5.98e-16 i**. These are local estimates, not independently
refined accurate roots or certified error bounds. The stress point is not a root;
its recorded frequency-equivalent determinant perturbation is explicitly labeled
as such and is not presented as a root shift. No root algorithm is changed.

## 8. Runtime of tighter tolerances

Measurements use medians of three CPU-time samples. Each Z sample performs
600,000 native Fortran evaluations, excluding Python/C call overhead. The six
arguments cover series, continued-fraction, asymptotic and problematic stress
regimes. Each tensor sample evaluates one species 100 times with fixed N=1024
and identical Gamma arrays. The test-only instrumented variant has measurable
instrumentation overhead, so tolerance cost should primarily be compared against
its own 1e-7 baseline.

| Evaluator / tolerance | Z CPU seconds | Z ratio to copy 1e-7 | Tensor CPU seconds | Tensor ratio to copy 1e-7 |
| --- | ---: | ---: | ---: | ---: |
| Production historical 1e-7 | 0.05132 | 0.899 | 0.07577 | 0.985 |
| Instrumented 1e-7 | 0.05708 | 1.000 | 0.07693 | 1.000 |
| Instrumented 1e-10 | 0.07163 | 1.255 | 0.07804 | 1.014 |
| Instrumented 1e-13 | 0.07643 | 1.339 | 0.07965 | 1.035 |
| Instrumented 1e-15 | 0.08195 | 1.436 | 0.08036 | 1.044 |

A tolerance-only tightening costs about **44%** in this Z microbenchmark and
about **4.4%** in this fixed1024 tensor workload, while failing to fix the dominant
asymptotic error. These timings cannot predict the cost of a future robust
Fortran Faddeeva backend. The Python reference callback is only a diagnostic
intervention and is not a proposed production implementation or fair production
performance comparator. Results depend on the argument distribution and hardware.

## 9. Reproducibility, tests and preserved behavior

Recorded toolchain: GNU Fortran 11.4.0, Python 3.10.12, NumPy 2.2.6, SciPy 1.15.3,
mpmath 1.4.1, x86-64 Linux. Probe libraries compile in temporary directories with
`-fPIC -fdefault-real-8 -O3 -funroll-loops -ffast-math`; modules requiring IEEE
semantics and the probe wrapper use trailing `-fno-fast-math`, matching the normal
module policy. Checked builds use `-O0 -g -fcheck=all` and
`-ffpe-trap=invalid,zero,overflow`. No production objects are replaced by audit
objects. The actual production-object comparison executable is linked separately.

Commands:

```sh
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_plasma_dispersion_audit.py' -v
```

The baseline **92 tests pass**, including the exact pinned pre-normalization
fixed128 regression. Ten additional characterization tests cover the audit
interventions and discovered defects. The final complete run passes **102 tests**.
The tests document legacy defects rather than demand reference accuracy from the
current legacy algorithm. Future accurate-backend tests must use their own quality
criteria while these checks remain associated with compatibility mode.

For the full measured audit, install mpmath in an importable environment and run:

```sh
MPLCONFIGDIR=/tmp/plume-mpl python3 tests/run_plasma_dispersion_audit.py \
  --high-precision --output /tmp/plume-z-results.json \
  --arguments-csv /tmp/plume-z-arguments.csv
```

For this run mpmath was installed only into `/tmp/plume-z-python`, and the command
used `PYTHONPATH=/tmp/plume-z-python`. The ordinary unittest suite does not require
mpmath. Omitting `--high-precision` omits only the arbitrary-precision measurements.
The optional CSV exports every recovered signed complex argument at binary64
round-trip precision; the JSON retains the dominant arguments, parameter states,
source hashes, summaries and budgets. An audit run does not rewrite the report
or production files.

Added files:

- `tests/plasma_dispersion_audit.py`: temporary builds, instrumented source extraction,
  exact-argument reference callback, frozen-input susceptibility probe API, independent
  affine terms and arbitrary-precision references.
- `tests/plasma_dispersion_probe.f90`: C-ABI diagnostic entry points, IEEE recording,
  and direct inclusion of the unchanged species equations.
- `tests/plasma_z_direct_probe.f90`: direct calls linked against original objects.
- `tests/run_plasma_dispersion_audit.py`: deterministic grid/stress reconstruction,
  runtime measurements, production checks, root sensitivity and result export.
- `tests/test_plasma_dispersion_audit.py`: ten standalone characterization tests.
- `tests/PLASMA_DISPERSION_AUDIT_RESULTS.json`: measured data and source fingerprints.
- `tests/PLASMA_DISPERSION_AUDIT.md`: this report.

A before/after checksum comparison verifies all pre-existing tracked files and
all source files unchanged. Pre-existing uncommitted and unrelated modifications
are retained. No expected regression values or old tests are edited.

## 10. Recommendation for Stage J1.B.2

Proceed with a **targeted Z-backend accuracy and failure-contract correction**:

1. Preserve a selectable, bitwise legacy pathway. A separately compiled accurate
   backend should be tested through the existing callback interface before any
   default change. An independent Z accuracy control is preferable to coupling
   Z accuracy implicitly to the Bessel truncation mode. Retain the exact fixed128
   baseline as a compatibility requirement.
2. Use a robust Faddeeva/continued-Z algorithm across the entire required domain,
   with particular validation near radius three, the real axis, the continued-
   fraction cap, zeros and the lower-half-plane limits. A tolerance-only edit is
   ruled out by the measurements. If extending the current algorithms is chosen,
   both their regime boundaries and convergence/range logic require correction;
   changing the scalar tolerance alone cannot meet an accuracy contract.
3. Define a scale-aware target (initially about 1e-13 on representable ordinary
   Z values) and return explicit status for invalid k, nonconvergence and true
   range failure. Distinguish mathematical overflow from unsupported-domain exits.
   Do not add damping, floors or offsets. A finite unassigned return must never
   pass as a valid response.
4. Repeat the fixed-Gamma intervention and complete production regression before
   integration. Benchmark a real compiled backend: the instrumented-copy and
   Python callback timings are not replacements for that measurement. Review
   licensing before adopting an external implementation.
5. Retain the residual electron conditioning control. Its 1e-10-scale component
   sensitivity warrants a later targeted cancellation study if required by
   scientific accuracy goals, but does not justify rewriting the susceptibility
   equations to repair the primary 7.62e-4 discrepancy.

The collisional zero-rate parenthesis inconsistency needs a **separate derivation
and scope approval** before a formula correction. It should have an explicit
anisotropic nu->0 acceptance test; the present audit only demonstrates and
isolates it. No such correction, numerical replacement, or Stage J1.B.2 work
has been performed here.
