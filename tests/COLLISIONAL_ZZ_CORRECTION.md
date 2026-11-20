# Collisional zz susceptibility and harmonic correction

## Result and scope

The longitudinal response of the existing non-conservative Krook model
`C[f1]=-nu*f1` is corrected in `src/kinetic_collisional_species.inc`.
Both the historical fixed128 and dynamic/adaptive evaluators include this source.
The correction covers the total chi_zz and its signed n=0,+1,-1 diagnostics.

The original audit identified misplaced parentheses in the anisotropic term.
Planning then tested the proposed grouping correction against a direct kinetic
integral: grouping alone recovered nu=0 but left absolute errors of about 9.9
and 93.7 at nu=.01 and .1 in the planning example. Additional finite-rate
coefficients and the zero harmonic therefore required correction.

Only collisional zz production expressions are changed. The other tensor
components, collisionless expressions, Z evaluator, normalization, Bessel
selection/cache, collision-rate conversion, root method, and public input/API
remain unchanged. Collisional numerical results intentionally change; retaining
the incorrect results is not a compatibility mode. The exact pinned collisionless
regression remains required and passes. All changes remain local and uncommitted.

The [original Z audit](PLASMA_DISPERSION_AUDIT.md) and its
[measurement JSON](PLASMA_DISPERSION_AUDIT_RESULTS.json) are historical evidence;
the JSON has not been regenerated or edited.

## Independent longitudinal derivation

Let all velocities be in reference Alfvén units and define

\[
 v_t^2=\mu_s\tau_s^{(d)},\quad g=\mu_s/Q_s,\quad
 A_s=D_s\mu_s/(Q_s^2a^2),\quad W=w+i\nu,\quad k=k_\parallel d_R.
\]

The equilibrium parallel Gaussian has mean u and thermal width v_t. For an
applied E_z, Faraday's law uses the physical frequency w, whereas the kinetic
resolvent uses W. Keeping both the electric force and the perturbed magnetic
force gives the signed-harmonic response

\[
 \chi_{zz,n}=-\frac{2A_s\Gamma_n}{wv_t^2}
 \left\langle
 \frac{v\left[(v-u)+\frac{ng}{w}\{(\alpha^{-1}-1)v+u\}\right]}
 {W-kv-ng}
 \right\rangle,
\]

where angle brackets denote the normalized Gaussian average. This retains the
full equilibrium-frame current velocity v and the physical drift. No frequency
shift is inserted into Faraday's law or the current prefactor.

For r=W-ng and zeta=(r-ku)/(kv_t), the independent moment oracle uses

\[
 M_0=-Z_k(\zeta)/(kv_t),\qquad
 M_1=(rM_0-1)/k,\qquad M_2=(rM_1-u)/k,
\]
\[
 \chi_{zz,n}=-\frac{2A_s\Gamma_n}{wv_t^2}
 \left[\left(1+\frac{ng}{w}(\alpha^{-1}-1)\right)M_2
       +\left(\frac{ng}{w}-1\right)uM_1\right].
\]

These recurrences follow directly by dividing velocity polynomials by the
resolvent. They provide an independently structured reference rather than a
second implementation of the expanded production polynomial. Tests also
integrate the original velocity-space expression with 128- and 256-point
Gauss-Hermite quadrature where Im(W)>0. For the damped continuation, the moment
oracle uses the same signed Landau convention as PLUME. No real-axis quadrature
is used as a lower-half-plane continuation.

## Corrected paired and zero-harmonic expressions

Define V=ku, m=ng and C=1/(kv_t), and retain the existing normalization

\[
 N_{zz}=\frac{2D_s}{a^2w^2Q_s^2k^2\alpha\tau_s^{(d)}}.
\]

The positive/negative pair is

\[
 N_{zz}\Gamma_n[P+C\{(Z_n+Z_{-n})A+m(Z_n-Z_{-n})B\}],
\]
\[
 P=2w(\alpha W-V)-2(1-\alpha)m^2,
\]
\[
 A=\alpha wW(W-V)+m^2[\alpha w-2(1-\alpha)W-\alpha V],
\]
\[
 B=(1-\alpha)W^2-2\alpha wW+\alpha V(W+w)+(1-\alpha)m^2.
\]

The zero-harmonic bulk term is

\[
 N_{zz}\Gamma_0[w(\alpha W-V)+CZ_0\alpha wW(W-V)],
\]

followed by the unchanged separate drift term N_zz*w*V.

The source keeps its existing expanded arithmetic structure. The changes are:

- Move the anisotropic m^2 term inside the odd-Z bracket.
- Use i*nu*alpha, rather than i*nu*(alpha-1), in the leading zz constant and
  leading even-Z frequency factor, including their n=0 counterparts.
- Change the odd-Z coefficient from -i*nu*w*alpha to -2*i*nu*w*alpha.

At nu=0 the bulk tensor recovers the collisionless expression. No correction is
applied to the analogous factors in other components; those formulas are
outside the authorized scope of this change.

## Signed harmonic diagnostics and drift allocation

For signed m=ng, the n=0,+1,-1 diagnostic now evaluates

\[
 N_{zz}\Gamma_{|n|}\left[
 (\alpha w+(1-\alpha)m)(W-m)+mV
 +CZ_n\{(\alpha w+(1-\alpha)m)(W-m)^2
       +\alpha(m-w)V(W-m)\}\right].
\]

This is the complete individual longitudinal moment response. It corrects the
duplicated +1 assignment, wrong negative-harmonic Z argument, and nonresonant
splitting. All three zz diagnostic slots are explicitly assigned at every call.
The signed gyrofrequency, drift and Z argument remain signed throughout.

These diagnostics contain each harmonic's full nonresonant terms. In contrast,
the bulk implementation analytically regroups a complete drift contribution
outside its truncated harmonic sum. With

\[
 S_N=\Gamma_0+2\sum_{n=1}^N\Gamma_n,
\]

the exact reconstruction for the supplied finite coefficient array is

\[
 \chi_{zz}^{\rm bulk}=
 \sum_{n=-N}^N\chi_{zz,n}+N_{zz}wV(1-S_N).
\]

Thus the diagnostic sum must not receive a second full drift correction.
A deliberately truncated N=1 test keeps the remainder large, making accidental
omission or double counting detectable. This identity also accounts for the
small normalization error of historical fixed128 Bessel coefficients.

Existing collisionless diagnostic formulas are not modified. Their nonresonant
allocation need not equal these corrected raw collisional entries at complex
frequency. Their real-frequency zz dissipative parts have the same zero-rate
limit: differences in the nonresonant real part do not contribute there.

## Quantitative checks

For species `(tau_d,mu,alpha,Q,D,u)=(2,.25,1.3,.5,.05,.4)`, use
w=.73+.5i, kpar=-.31, kperp=3, a=.02 and N=128. Both old and corrected probes
use identical SciPy Gamma coefficients. The last column substitutes SciPy Z
only to isolate the corrected algebra from the unchanged Z approximation.

| nu | Historical chi_zz with production Z | Corrected chi_zz with production Z | Corrected chi_zz with SciPy Z |
| ---: | --- | --- | --- |
| 0 | 6288.34507699 - 13230.1530178 i | 40.83817923 + 65.13660389 i | 40.8377509160477 + 65.1365605870392 i |
| .01 | 6281.57541720 - 13204.7039446 i | 40.66472153 + 64.84596392 i | 40.6647603155415 + 64.8456538388319 i |
| .1 | 6218.78173726 - 12985.0911603 i | 39.17959004 + 62.38458678 i | 39.1896297507456 + 62.3759399710877 i |

The independent quadrature values are respectively
`40.8377509160470+65.1365605870382i`,
`40.6647603155409+64.8456538388307i`, and
`39.1896297507452+62.3759399710870i`.
Absolute differences after Z substitution are at most **1.42e-12** in this table.
The remaining differences with production Z are inherited Z-evaluation error;
this change does not implement Stage J1.B.2.

Across the focused automated cases, componentwise scaled errors use
`abs(actual-reference)/max(1,abs(reference))`:

| Check | Maximum observed discrepancy |
| --- | ---: |
| Total zz versus independent velocity moments | 1.825e-12 scaled |
| Total zz / moment versus direct quadrature checks | 1.633e-12 scaled |
| Signed harmonics versus moments/quadrature | 1.036e-10 scaled |
| Full tensor zero-rate limit | 1.204e-14 scaled |
| Production dispatch, coefficient arrays and normalization comparisons | 3.154e-14 scaled |
| Eight other elements versus frozen original audit values | 3.442e-15 scaled |

The larger signed-harmonic discrepancy occurs in cancellation-sensitive
binary64 comparisons, including electron parameters. The generic scalar and
harmonic checks use rtol=2e-10 and atol=2e-11; the isolated audit's zero-rate
check retains rtol=2e-13 and atol=2e-11. The independent quadrature is checked
for convergence before it is used as a reference.

An additional before/after capture of nine tensors (three species, three rates)
finds **identical binary64 values** for all eight other elements and their
harmonic diagnostics. This is separate from the frozen-audit comparison above,
whose different normalization/compiler contexts can differ by last bits.

### Root and heating precision

Zero-rate production tests compare converged frequencies, fields, density,
velocity and polarization with the collisionless calculation. Both fixed and
adaptive modes and checked builds are exercised. The secant method is unchanged;
its final iteration count can differ for equivalent floating-point expressions.

The existing wave-energy normalization uses a frequency finite difference of
1e-6. Equivalent zz polynomials differing at roundoff can therefore produce a
small difference in normalized heating. The measured maximum in the new
zero-rate comparison is **2.970e-6 relative**. Its heating-only acceptance is
rtol=1e-5, atol=1e-10; all other records use the unchanged normalization-test
comparison. No original regression tolerance is relaxed. The exact pinned
collisionless legacy comparison remains bitwise exact.

## Tests, integration and files

Eight focused tests in `test_collisional_zz.py` cover:

- Independent moments and converged velocity quadrature, including n=0,+1,-1.
- Isotropic/anisotropic distributions, proton/electron/alpha parameters,
  both charge signs, signed drifts and both parallel-wavevector signs.
- Complex frequencies in both half-planes; zero, decreasing and finite rates.
- Identical frozen Bessel/Z inputs and a deliberate N=1 drift reconstruction.
- Poisoned diagnostic arrays, finite slots, sign distinction and repeat calls.
- Actual fixed128 and adaptive dispatch in rho_ref and d_ref, using the
  coefficient arrays from the production evaluators.
- Zero-rate full tensor, root and supplementary-diagnostic comparisons.
- Other components against frozen pre-correction values.

Production and bounds-checked/FPE-trapping builds are tested. The diagnostic
Fortran probe's new low-output array is test-only; it does not change a PLUME
production API or input format. The historical audit's intentional zz defect
assertion is replaced with a tight zero-rate equality assertion. All Z defect
and accuracy checks remain in place.

Commands:

```sh
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_collisional_zz.py' -v
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
```

Focused result: **8 tests pass**. Final full result: **110 tests pass** in
51.1 seconds, including the exact pinned collisionless legacy regression.
GNU Fortran 11.4.0 and the existing production flags are used; checked probes
use `-O0 -g -fdefault-real-8 -fcheck=all -ffpe-trap=invalid,zero,overflow`.

Files changed or added in this correction:

- `src/kinetic_collisional_species.inc`: total and low-order collisional zz fixes.
- `tests/collisional_zz_reference.py`: independent moment and velocity-integral oracles.
- `tests/collisional_zz_probe.f90`: production-path tensors, diagnostics and coefficients.
- `tests/test_collisional_zz.py`: focused mathematical and production tests.
- `tests/plasma_dispersion_probe.f90`, `tests/plasma_dispersion_audit.py`: expose
  and poison test-only diagnostic arrays, retaining the existing audit interfaces
  at Python call sites.
- `tests/test_plasma_dispersion_audit.py`: require the corrected zero-rate limit.
- `tests/run_plasma_dispersion_audit.py`: label the old parenthesis prediction
  as a historical counterfactual when run against corrected source.
- `tests/PLASMA_DISPERSION_AUDIT.md`, `output.md`: link this correction and
  explain the signed diagnostic allocation.
- `tests/COLLISIONAL_ZZ_CORRECTION.md`: this derivation and implementation report.

Remaining limits are unchanged: the legacy Z approximation can dominate physical
accuracy; other collisional tensor components are not independently rederived
here; cold-containing collisional calculations remain unsupported. No new
collision model, root algorithm, normalization, or scan capability is introduced.
