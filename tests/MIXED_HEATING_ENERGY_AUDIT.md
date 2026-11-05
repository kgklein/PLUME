# Stage 4B.3A: heating and wave-energy audit

## Executive summary

**Audit only. No production source, heating formula, input default, output writer,
Python loader, or earlier regression expectation is changed.** The numerical
probes call the existing compiled routines or evaluate independently structured
references. The recommendations below are for a subsequent implementation task.

The historical species diagnostic is a **dimensionless absorption/emission
quadratic form divided by a wave-energy-like normalizer**, not dimensional
power or an energy fraction per wave period. Its seven-channel implementation
is recoverable without changing the susceptibility physics. At a common real
frequency and with initialized harmonic storage, the four n=0 energy channels
sum to the n=0 total. The n=+1 and -1 outputs use only transverse electric
fields; they and n=0 do **not** exhaust the species total. Summing all seven
columns also double-counts n=0.

Independent checks confirm these implementation issues:

1. The historical energy normalizer lacks the vacuum electric contribution
   |E|^2 because it differentiates sum(chi), not I+sum(chi).
2. Its one-sided relative-1e-6 derivative is subject to truncation and
   cancellation; a smaller step does not reliably improve the energy contraction.
3. The species total uses chi at omega_r, while the seven channels use harmonic
   arrays overwritten at 1.000001*omega_r during the derivative calculation.
4. Fixed128 collisionless n=0 storage has five unassigned entries per species;
   its full n=0 quadratic form can read uninitialized/stale values. Adaptive
   evaluation clears these entries. This is separate from Bessel truncation.
5. `Ew` in some scan/radial output lists is not the local `ewave` and is not
   assigned on successful evaluations. Heating-only/low_n=false also has a
   format-reversion problem that writes the extra `Ew` on a separate record.
6. The current Python loader reads seven channels with a stride of six, and its
   heating-only total-power starting offset is wrong. Some written documentation
   has the same six-versus-seven inconsistency; the Fortran low_n=true files
   themselves retain seven channels in the correct species order.
7. A separate finite-rate collisional check finds discrepancies in oblique xz/yz
   response for anisotropic/drifting populations against direct Krook Vlasov
   quadrature. The corrected longitudinal row agrees, as do stationary isotropic
   controls. This does not invalidate the earlier zz correction, and nothing in
   that correction is revisited or modified. The remaining off-diagonal issue
   needs a separately scoped operator/tensor audit before claiming full
   collisional heating validation.

Warm-reference mixtures already retain warm harmonic arrays away from parallel
propagation, but cold-reference and exact-parallel mixed dispatch deliberately
leave their low-harmonic arrays unavailable. Input, scan, and diagnostic guards
also reject mixed heating. Restoring it therefore requires species-resolved
availability, full-plasma fields and energy, and delivery of **all existing
channels**, not merely enabling the old `calc_eigen` guard.

## 1. Existing implementation and dependency map

Locations below are routine/include names, avoiding fragile revision line ranges.

```text
functions::read_in_params -> spec_read -> normalization adapter / cold validation
    low_n -> allocate susc_low(species,3,3,-1:1)
    scan_read -> heating/eigen flags
                 |
disprels::disp -> disp0 / disp1 / disp_mixed
    kinetic_species.inc / kinetic_collisional_species.inc
    adaptive_kinetic::evaluate (same includes, per-species cache)
    cold_reference::calc_warm_auxiliary / disp_mixed_parallel / cold kernel
                 |
    susc(species,row,column): ordinary chi_s, no vacuum identity
    susc_low(species,row,column,-1:1): selected harmonic contributions
    lam = I + sum(chi_s) + n n^T - n^2 I
                 |
    calc_eigen(omega,...,eigen_L,heat_L) [all-kinetic only]
       fields at complex omega from current lam
       heat_L -> disp(omega_r) -> species chi_A -> raw species Q
       heat_L -> disp(1.000001*omega_r) -> finite difference chi_H
       local ewave -> Ps and seven Ps_split values
                 |
    test_disp: prints Ps and only first four Ps_split channels
    om_scan / om_double_scan: write Ps and all seven channels if low_n=true
    om_radial: separate older total-heating layout (excluded from restoration)
    fpc callers: consume calc_eigen internally (excluded from this audit's changes)
```

`calc_eigen` aliases the thermal solver state. Its field/moment portions and
real-frequency derivative calls share global `susc`, `susc_low` and `lam`.
Heating leaves those global tensors at the shifted real frequency, not at the
input complex eigenfrequency. Scan tensor output explicitly reevaluates `disp`
at the mode when requested, but other consumers must not assume an unchanged
snapshot after heating. The new Stage 4B.1/4B.2 field/moment APIs are separate:
`calc_mixed_eigen` currently rejects heating and assigns its heating outputs NaN.
Cold-reference states cannot safely enter the old thermal-reference code.

| Diagnostic | Definition / physical interpretation | Dependencies and output | Current status / mixed obstruction |
| --- | --- | --- | --- |
| Species total | Q_s/H_old, Q_s=E† chi_A,s(omega_r) E; normalized lab-frame absorption/emission | Species chi; full-plasma E/B; full-plasma denominator. `Ps(s)` is real dimensionless | Calculated for all-kinetic inputs; mixed guards and old thermal aliases block use |
| Total over populations | sum_s Ps(s); not a separate independent calculator | Common denominator and all species Ps | Available by summing legacy species columns; summing must not ignore missing species silently |
| Four n=0 channels | yy, yz, zy, zz terms described below | Only that species' harmonic tensor and full-plasma E, normalized by common H | Existing all-kinetic path, subject to shifted-frequency/initialization caveats |
| n=0 total | E† chi_A,s,0 E/H_old | `Ps_split(5,s)` | Full quadratic form reads all nine harmonic entries |
| Signed cyclotron outputs | E_perp† chi_A,s,+1/-1 E_perp/H_old | `Ps_split(6:7,s)`; signed harmonic ordering fixed | Transverse projection only; cold-reference/parallel harmonic arrays not exposed |
| Energy normalizer | |B_complex|^2 + E† D_forward E, D_forward approximates d[w sum(chi_H)]/dw | Local `ewave`, used to divide Ps/Ps_split | Not a complete dielectric wave energy; not returned to caller |
| Written `Ew` | Caller scalar intended to represent wave energy | Some low_n=false and radial writes | Unassigned on successful paths; not an exported `ewave` |
| FPC | Uses calc_eigen and has additional output/thresholding code | `src/fpc.f90` | Only call/availability mapping audited; no FPC implementation or validation added |

There is no entropy-production calculation, thermalization operator, or explicit
background bulk/internal-energy split in these species `Ps`. In drifting or
anisotropic populations negative absorption represents available free energy;
calling every signed value irreversible thermodynamic heating would overstate
what is computed. The existing names/column positions must nevertheless be
preserved and their definitions clarified.

## 2. Mathematical definition and normalization

Coordinates are B0 along z and k=(kperp,0,kpar). The Fourier convention is
exp(i*k.r-i*omega*t), with positive Im(omega) growth. Write

    w = omega/Omega_R, kappa=k*d_R, a=vA,R/c,
    E = deltaE/deltaEx, B = deltaB/deltaEx, Ex=1.

PLUME stores ordinary dimensionless susceptibility, with

    epsilon = I + sum_s chi_s
    chi_H = (chi + chi†)/2
    chi_A = (chi - chi†)/(2i).

At real frequency wr=omega_r/Omega_R, the time-averaged lab-frame species work
is

    P_s = Omega_R*wr*|deltaEx|^2/(8*pi) * Q_s,
    Q_s = E† chi_A,s(wr,k) E.

Here `Q_s` in the audit equations is a scalar absorption form, **not** the input
charge ratio named `Qs`. Anti-Hermitian refers to this standard Hermitian
absorptive matrix, not to entrywise imaginary parts of every tensor element.
For arbitrary complex E, Q_s=Im(E† chi_s E) is real. Individual cross channels
can be signed and can cancel.

For a weakly damped/growing isolated mode, the conventional leading-order
real-frequency wave energy is

    W = |deltaEx|^2/(16*pi) * H,
    H = |B_r|^2 + |E|^2
        + sum_s E† d[wr*chi_H,s(wr,k)]/dwr E,
    B_r = (kappa cross E)/(a*wr).

The derivative holds physical k, densities, temperatures, anisotropies, drifts,
magnetic field and reference definitions fixed. Differentiate laboratory
`w*chi`, not Doppler-shifted `delta*chi`: d(delta)/dw=1 at fixed drift, and the
drift couplings are already part of chi. It does not differentiate the
mode eigenvector or B with respect to frequency. Since omega=Omega_R*w,
`d(omega*epsilon)/domega=d(w*epsilon)/dw`; no extra Omega_R, beta or thermal
factor belongs in H. Contributions from cold species must be included in the
sum even when their kinetic heating columns are unavailable.

Historical rho_ref uses K=k*rho_R with the **perpendicular** gyroradius:

    K=kappa*sqrt(beta_R*alpha_R), vtp_solver=a*sqrt(beta_R),
    K/(vtp_solver*sqrt(alpha_R))=kappa/a.

Thus the same physical field and derivative normalization applies to rho_ref
and d_ref. Cold-reference auxiliary B=b_star and V=a*sqrt(b_star) likewise
cancel in the physical chi, kappa/a, and frequency derivative. The auxiliary
scale is not an energy-bearing species and contributes no extra term to W.

For nondrifting cold populations, an independent test verifies the positive
field-plus-fluid kinetic-energy identity

    H = |E|^2 + |B_r|^2 + sum_s D_s/(mu_s*a^2) * |v_s|^2,
    v_s = deltaU_s/(c*deltaEx/B0).

This follows from the time-averaged fluid kinetic energy n_s*m_s*|deltaU_s|^2/4.
It independently checks vacuum, magnetic and cold susceptibility derivative
coefficients. It is **not** imposed as a general positive-energy decomposition
for streaming populations; their canonical wave energy can be negative.

### What the historical normalizer means

The executable expression is

    H_old = |B_complex|^2 + E† D_old E,
    D_old = [1.000001*sum(chi_H(1.000001*wr))
                          -sum(chi_H(wr))]/1e-6.

The missing vacuum |E|^2 cannot be recovered by interpreting `susc` as epsilon:
the source and tensor reconstruction tests establish that it contains chi only.
This omission is often small in the low-frequency/low-vA/c regime and could
serve as an approximation there, but it is not identically zero. No explicit
full-dielectric justification for dropping it was found in the implementation. The
isolated fully cold moderate-a test has a vacuum contribution about half of H.
No fully cold wave energy is currently advertised by production; that test
shows why simply reusing H_old for it would be wrong.

The existing B is obtained at the full complex mode w, while the derivative
uses wr. For fixed E,
`|B_complex|^2/|B_r|^2=wr^2/(wr^2+gamma_hat^2)`.
This is second order in gamma/wr for weak modes but substantial for the reactive
benchmarks. The complex eigenvector can be used as a leading-order phasor in a
weak-mode estimate; this does not make the result an exact conserved energy of
an arbitrary complex-frequency mode.

Historical output is `Ps(s)=Q_s/H_old`. With the complete H, one would have

    gamma_hat approximately -wr * sum_s(Q_s/H),
    gamma_hat/wr approximately -sum_s(Q_s/H).

This explains the existing plot labels gamma_s/|omega_r|. It is neither
normalized dimensional power nor an energy loss fraction per period. For a
positive period T=2*pi/|omega_r|, the latter is approximately
`4*pi*sign(wr)*Q_s/H`. Comments claiming power "in one wave period" need this
clarification; introducing a 2*pi or omega factor into the existing Ps columns
would change their meaning and is not recommended.

## 3. Harmonic and energy-channel decomposition

Let C_s,n=chi_s,n(wr), and E_perp=(Ex,Ey,0). The seven outputs, in their actual
Fortran order, are:

| Slot | Source interpretation | Raw numerator before the common denominator |
| --- | --- | --- |
| 1 | TTD yy | abs(Ey)^2 Im(C_0,yy) |
| 2 | TTD yz | Im(Ey* Ez C_0,yz) |
| 3 | LD zy | Im(Ez* Ey C_0,zy) |
| 4 | LD zz | abs(Ez)^2 Im(C_0,zz) |
| 5 | n=0 total | E† (C_0-C_0†)/(2i) E |
| 6 | n=+1 transverse cyclotron | E_perp† (C_+1-C_+1†)/(2i) E_perp |
| 7 | n=-1 transverse cyclotron | E_perp† (C_-1-C_-1†)/(2i) E_perp |

Stars in the cross terms denote complex conjugation. The yy/yz grouping gives
transverse-electric n=0 work conventionally called TTD; zy/zz gives
parallel-electric n=0 work conventionally called LD. Interference cross terms
are essential; these are not two independently positive entropy-production
rates. At zero-initialized n=0 transverse-x entries, slots 1+2+3+4 equal slot 5.
These operations require only each kinetic species' harmonic tensor plus the
full-plasma eigenvector. Other species enter the field and common denominator,
not the species' raw harmonic decomposition.

Slots 5+6+7 omit higher |n| and omit the longitudinal/cross-field part of n=+/-1.
They therefore need not sum to the total. A measured off-root oblique proton
case has 28.4487% of total work outside these three displayed totals: 2.76547%
from higher harmonics and 25.6832% from excluded +/-1 longitudinal/cross terms.
Summing slots 1 through 7 is worse: it counts n=0 twice. No restoration should
replace this documented projected decomposition by a different full-harmonic
closure just to force equality.

The collisionless signed zz diagnostic allocates nonresonant real terms
differently from an independent Stix harmonic grouping. Such terms do not
contribute to the absorptive quadratic form at **real** frequency. The tests
compare the complete absorptive matrix per signed harmonic, not an unsupported
claim that every complex harmonic tensor has a unique nonresonant allocation.
The existing collisional signed zz allocation has its own prior derivation;
that correction remains unchanged. No harmonic heating interpretation is
extended to arbitrary complex-frequency harmonic tensors.

### Frequency and storage consistency

`calc_eigen` first forms species totals at wr, then calls disp at
1.000001*wr for D_old, then reads `susc_low` for every channel. The observed
seven-channel relative bias from this shift is 3.53592e-6 in the chosen
three-warm-species case. Preserving the mathematical decomposition requires
snapshots at one explicitly selected real frequency before derivative sampling.
This is a correction distinct from mixed-species availability.

In the fixed128 collisionless include, n=0 yy, yz, zy and zz are assigned, while
xx, xy, xz, yx and zx are left untouched. The allocation in `functions` supplies
no initializer. A controlled probe fills scratch with 2i and demonstrates those
five entries survive; n=0 channel closure is wrong by 4.2 in the test quadratic
form. Adaptive evaluation initializes its local low array to zero and passes.
The collisional include explicitly clears the unused n=0 entries. Neither fresh
OS memory being zero nor a benign prior call is a sufficient numerical contract.

## 4. Mixed-species semantics and ingredient availability

| Composition/path | Species total chi | Warm selected-harmonic arrays now | Needed restoration |
| --- | --- | --- | --- |
| All-kinetic rho_ref/d_ref, oblique | Available | Available, with fixed scratch caveat | Retain compatibility; separately approve corrections |
| Warm-reference mixed, oblique | Available, including cold contributions | Warm arrays available; cold arrays NaN | Per-species heating availability, shared full-plasma fields/energy |
| Cold-reference mixed, oblique | Available | Production auxiliary path disables/invalidates low arrays | Expose compute_low in the shared auxiliary kinetic interface |
| Any mixed exact-parallel path | Available | Entire low array deliberately NaN | Populate existing analytical n=0 and +/-1 limits, without a kperp floor |
| Fully cold | Available | Kinetic harmonics not applicable | Kinetic heating NaN; energy-only diagnostic separately qualified |

An isolated call to the existing auxiliary/adaptive API with `compute_low=true`
recovers the warm arrays, including interleaved species and signed harmonics,
and agrees with the independent reference. No production flag or dispatch was
changed to obtain this evidence. Fixed128's auxiliary core currently hardcodes
`low_n=false`, so its stable future interface must expose optional low outputs
without duplicating equations or perturbing the compiler-sensitive legacy unit.

Future output must use the full plasma E/B for every warm species. Cold species'
**kinetic** total, harmonic and channel values must be NaN/unavailable in their
existing slots, not zero. Cold susceptibility derivatives remain part of W.
Sum raw real-frequency work over all species for an energy-consistency test;
collisionless cold chi is Hermitian away from poles, so its absorptive work is
then computed zero. This mathematical fact does not turn unavailable cold
kinetic heating columns into measured zeros. At complex frequency, cold work
includes reversible storage/growth terms and cannot be relabeled kinetic heat.

For example, direct complex-frequency work contains
`wr*E†chi_A(w)E + gamma_hat*E†chi_H(w)E`, not simply wr times the absorptive
quadratic form. The audit uses the real-frequency expression for weak-mode
energy consistency, rather than confusing those quantities.

## 5. Independent energy-consistency checks

At a real, weakly non-Hermitian isolated mode,

    E† epsilon_H E approximately |B_r|^2,
    H approximately wr * E† (d Lambda_H/dwr) E.

The electromagnetic n n^T-n^2 I term is real at wr. Expanding the complex
wave equation to first order in gamma_hat therefore gives

    0 approximately E† sum(chi_A,s) E
                    + gamma_hat * E† (d Lambda_H/dwr) E,
    gamma_hat approximately -wr * sum(Q_s)/H = -P_total/(2*W*Omega_R).

The last equality uses mutually consistent **dimensional** P and W defined
above. The normalized equation is gamma_hat=-wr*sum(Q_s)/H. No measured gamma
is used to set H or Q: roots come from the existing solver, independent SciPy
root refinement checks the same physical branch, Q comes from susceptibility,
and H comes from independent derivatives.

The real-frequency relation includes all relevant species' field-particle work.
Collisionless cold species contribute reversible inertia/streaming terms to H,
and no absorptive work away from real poles. An explicit kinetic mask is needed
for kinetic-heating output; using an indiscriminate `nansum` is not sufficient.
It would turn a fully cold plasma, or a failed warm diagnostic, into a misleading
physical zero. A missing warm contribution must make a claimed total unavailable.

The benchmark roots agree with the independent determinant roots within
1.14e-12 in absolute complex frequency. Their maximum independent scaled field
null residual is 6.77e-13. These are supporting checks, not a claim that a small
determinant residual by itself establishes correct heating or root accuracy.

The following table uses high-precision derivatives and the real-frequency
magnetic field. Values are dimensionless omega/Omega_R and the energy H.

| Case | Complex root | H | Relative weak-mode gamma error |
| --- | --- | ---: | ---: |
| kinetic | 0.0609871523619 -1.270997e-05 i | 20342.6801 | 2.192e-07 |
| warm_reference_mixed | 0.0615384195836 -5.564326e-06 i | 20384.8001 | 8.716e-08 |
| cold_reference_mixed | 0.0591673491216 -1.588769e-06 i | 20516.8792 | 1.657e-09 |
| fully_cold | 0.059722520022 -9.097903e-15 i | 20607.6296 | nondissipative: gamma is roundoff |
| interleaved | 0.0348316977935 -2.655359e-05 i | 115849.972 | 4.314e-06 |
| weak_growth | 0.322862839067 +0.001674958 i | 59263.7354 | 6.965e-05 |
| reactive_mixed | 0.30574584309 +0.07964933 i | -9993.5887 | 1.000e+00 |
| reactive_cold | 0.304781926478 +0.07916839 i | -9651.504 | 1.000e+00 |
| high_a_cold | 7.46110067295 -2.125421e-15 i | 53675.3078 | nondissipative: gamma is roundoff |

For the four weakly damped kinetic/mixed cases, the largest relative discrepancy
is 4.32e-6. The weakly growing anisotropic proton mode has gamma/wr=0.00519 and
relative error 6.97e-5. These are approximate energy balances, not exact imposed
identities. The fully cold stable root has gamma about -9.1e-15, consistent
with a real-frequency nondissipative mode; a relative error against that noise
is not a meaningful accuracy metric.

The two reactive beam cases have gamma/wr about 0.26. Their real-frequency
absorptive work is zero or exponentially tiny, while their growth is finite.
The formal H is negative. This does not establish a conserved negative energy
for the strongly unstable mode: the weak-mode expansion itself is inapplicable.
The gamma relation deliberately fails by order unity, as it should for these
benchmarks; no damping floor or closure adjustment is applied.

Raw species Q agrees with the arbitrary-precision independent susceptibility
reference to maximum `abs(diff)/max(1,abs(reference))=7.27e-16` over these cases.
Thus the principal ordinary heating-normalization issues here are not errors
in the absorptive quadratic form. The moderate-a fully cold benchmark is an
additional energy-normalization stress case, not a claim that production already
exports a cold wave-energy diagnostic. Its missing-vacuum fraction is 0.5000003.

The energy pieces below use the original species order, with labels in the
machine-readable case data. They are contributions to the real-frequency
quadratic expression, not independently positive thermodynamic energies for
general warm/streaming populations.

| Case | Vacuum | Magnetic | Species contributions in order |
| --- | ---: | ---: | --- |
| kinetic | 1.0012563 | 9708.5585 | 10812.267, -179.14646 |
| warm_reference_mixed | 1.0020029 | 9538.4899 | 10839.827, 5.4814841 |
| cold_reference_mixed | 1.0049271 | 10380.475 | 10323.932, -188.53348 |
| fully_cold | 1.0072072 | 10216.371 | 10384.752, 5.4995463 |
| interleaved | 1.4279675 | 51714.835 | 16187.324, 2600.7555, 22934.407, 22411.223 |

For example, cold reference protons supply about 10324 of H=20517 in the
cold-reference mixture. Omitting their derivative merely because their kinetic
heating is unavailable would lose about half the energy normalizer.

## 6. Derivative convergence and conditioning

The independent warm derivative differentiates the NHDS/Stix harmonic tensor
analytically at 65 decimal digits, using Z'(z)=-2[1+z Z(z)] before rounding to
binary64. Cold derivatives use high-precision differentiation of the rational
single-species tensor. A 60-versus-90-digit and 24-versus-32-harmonic check
has zero discrepancy after binary64 rounding in its representative cases.
All energy benchmarks have small Bessel arguments; the selected arbitrary-
precision harmonic range is not presented as sufficient for arbitrary lambda.

Centered differences at relative steps 1e-3, 1e-4, 1e-5, 1e-6 and 1e-7 use
**unchanged production chi** at prescribed real frequencies. Comparison is made
both before and after contraction with the mode E. This distinction matters:
a global matrix norm can be dominated by a large longitudinal entry whose
contribution to energy is suppressed by tiny Ez. Matrix convergence does not
by itself establish convergence of E† D E.

The historical comparison reproduces its actual frequency product
`wr*1.000001`, Hermitian sum-before-difference ordering and common denominator.
The raw historical denominator for all-kinetic cases is also inferred from
Q_s/Ps(s), so no newly rearranged formula is mistaken for exact legacy arithmetic.
For cold-containing cases the same expression is a **hypothetical reuse** of
H_old, not an existing enabled heating output.

| Case | Historical total derivative matrix error | Historical particle-energy contraction error | Best centered species-matrix error |
| --- | ---: | ---: | ---: |
| kinetic | 6.213e-07 | 6.308e-06 | 4.924e-10 |
| warm_reference_mixed | 9.971e-07 | 4.725e-08 | 3.400e-11 |
| cold_reference_mixed | 1.115e-07 | 2.382e-05 | 1.041e-09 |
| fully_cold | 1.000e-06 | 2.002e-08 | 2.301e-11 |
| interleaved | 1.212e-06 | 1.762e-05 | 2.423e-09 |

A representative cold-reference mixture illustrates the energy contraction's
step sensitivity (relative to its independent particle contribution):

| Relative step | Centered error | One-sided error |
| --- | ---: | ---: |
| 0.001 | 1.772e-08 | 3.904e-05 |
| 0.0001 | 8.827e-08 | 4.111e-06 |
| 1e-05 | 1.865e-06 | 3.963e-07 |
| 1e-06 | 4.021e-06 | 2.382e-05 |
| 1e-07 | 6.847e-05 | 1.879e-04 |

The one-sided 1e-6 energy error is about 2.38e-5 here; simply reducing the step
to 1e-7 increases it to about 1.88e-4. Centered steps near 1e-3--1e-4 are much
better for this projected quantity, even when a smaller step wins the matrix
norm comparison. Cancellation in the kinetic expressions and subtraction of
nearby large susceptibilities explain why a universal small step is unsafe.
No production derivative or special-function method is changed by this audit.

## 7. Additional collisional finding: isolate the operator before correcting it

The prior zz correction used the literal nonconservative Krook resolvent
`omega+i*nu-k.v` with electric/magnetic forcing at the **laboratory** omega.
This audit independently solves the full linearized Vlasov equation in
velocity/gyroangle space under the same convention. B is obtained from
Faraday at the laboratory frequency, not by replacing it with omega+i*nu.
The current includes the drift-density contribution. Thus the comparison
checks more than changing a Z argument or reusing the old polynomial.

A full three-electric-basis comparison at w=0.73, kpar=0.1, kperp=0.3,
a=0.02 and finite collision rates finds:

| Population | nu/Omega_R | Largest scaled xz/yz difference | Remaining seven entries | Reference grid difference |
| --- | ---: | ---: | ---: | ---: |
| Anisotropic stationary proton, tau=0.4, alpha=1.3 | 0.1 | 4.692e-2 | 9.72e-13 | 7.37e-14 |
| Anisotropic drifting electron, tau=1e-4, alpha=0.8, u=-0.2 | 0.06775 | 9.383e-2 | 2.49e-10 | 3.23e-9 |
| Isotropic stationary proton control | 0.1 | 8.64e-14 | 2.82e-13 | 7.91e-14 |
| Isotropic stationary electron control | 0.06775 | 2.87e-10 | 2.38e-10 | 2.10e-10 |

Scaled element error is abs(diff)/max(1,abs(production entry)). The comparison
uses 48x48/40-angle and 64x64/48-angle quadratures. Fixed128 and adaptive paths
show the same discrepancies. Representative complex matrices, not just these
maxima, are retained in the measurements JSON. A separate 32-versus-48-node
work calculation at two rates demonstrates that the effect survives quadrature
refinement. The observed scalar absorption discrepancies for its chosen E are
1.57e-4--2.31e-4 relative, even though component errors are larger.

The source always imposes xz=zx and yz=-zy. The independent finite-rate
anisotropic/drifting response does not obey those equalities; the calculated
longitudinal row and transverse 2x2 block agree while the copied oblique column
does not. This is compelling evidence for a residual **finite-rate
symmetry/operator inconsistency**, not a regression in the corrected zz, a
Bessel cutoff problem, or an inaccurate Z backend. It disappears in the
stationary isotropic controls. The finite-rate test initially failed a full
agreement assertion; subsequent basis/grid/control isolation established this
finding, and the final audit tests explicitly characterize it rather than
silently treating the mismatch as acceptable accuracy.

No new off-diagonal formulas are proposed as a drop-in patch. A literal
nonconservative Krook operator also raises a charge-continuity/reservoir closure
question. Any intended additional model terms must be stated before interpreting
it as a different collision model. Furthermore, `disp1`'s six-component tensor
assembly and compressed determinant assume the same symmetries. Correcting
only a heating contraction or one mirrored entry would leave inconsistent
root and energy calculations. A dedicated full-operator/tensor/assembly audit
is needed before a correction, with legacy compatibility explicitly retained.
That is broader than collisionless mixed-heating availability and remains
outside this stage. Cold-containing collisional calculations remain unsupported.

## 8. Output, compatibility and availability audit

### The actual low_n=true file schema

Let S=nspec, j=1..S, r=1..7, and

    A = 18+8*S when eigen=true; A=6 otherwise.
    species total P_j:                column A+j
    signed/channel slot r, species j: column A+S+7*(j-1)+r
    parameter slot t=1..6, species j: column A+8*S+6*(j-1)+t.

These are one-based indices. Existing eigen+heating rows have 18+22*S values;
heating-only rows have 6+14*S. `Ps_split(1:7,1:S)` is written in Fortran array
order, all seven slots per species, preserving +1 then -1. Negative charge
does not swap the signed-harmonic column labels. Correcting an interpretation
must not reorder those columns or remove the cross terms.

The current reader's seven values advance with `6*(j-1)` instead of seven.
For heating-only files its first P starts 2*S too late. Tests use actual serial
and 2D file layouts and then unique markers to establish the offsets
unambiguously, without relying on physically tiny/zero channels:

| Layout, S=3 | True P1 / read P1 | True species-2 yy / read value |
| --- | --- | --- |
| Eigen+heating, 84 columns | 1000 / 1000 | 1010 / 1009 |
| Heating only, 48 columns | 1000 / 1006 | 1010 / 1015 |

The same error appears in both scan types. The species parameter starts in the
current low_n=true reader are separate and must not be shifted along with the
channel fix. Automatic flag/species-count inference is inherently ambiguous;
explicit nspec/eigen/heating arguments remain advisable. Existing gnuplot
examples using hardcoded two-species Fortran columns have the correct seven
slots and label them as partial damping/growth contributions. They need NaN
masking/availability awareness for mixed output, not a column reorder.

### Single-point, energy-only and map behavior

- Single-point all-kinetic `test_disp` calls `calc_eigen` with both eigen/heat
  flags true. It prints total Ps and **only four** Ps_split entries per species,
  although seven exist internally. With low_n=false those channel outputs are
  not assigned by the calculator, yet the print block still executes; they
  must not be interpreted as calculated values. Cold-containing single-point calls currently
  print no heating. Preserve existing stdout labels and four-entry block if
  additional signed-harmonic output is later added separately.
- Serial and 2D scans call heating after each accepted output root, using that
  mode's complex E/B and current real-frequency work. Their near-zero-real-
  frequency overrides are inconsistent: serial uses +/-0.1 at a 1e-15 threshold;
  2D uses -0.1/+1 at 1e-7. These are historical placeholders, not physical work.
  A new mixed diagnostic must return unavailable instead of copying them.
- The local `ewave` is not returned by `calc_eigen`. `Ew` declared in scan and
  radial callers is unassigned on successful calculations. An audit-only build
  with `-finit-real=nan` demonstrates NaN in the intended legacy energy slot.
  It changes no repository source. The ordinary uninitialized value is not
  predictable and must not be interpreted as zero or a calculated wave energy.
- For eigen+heating/low_n=false, S=2 gives 49-value rows with the unassigned
  energy last. For heating-only/low_n=false the format accommodates 20 values
  but the write list supplies 21, giving alternating 20-value rows and one-value
  unassigned-energy records through format reversion. A loader must not treat
  those one-value records as independent roots.
- Radial output has the same caller-energy issue and a different older total-
  heating layout; radial/parameter restoration is explicitly deferred.
- Repeated maps (`option=4`) call `map_scan -> map_search`, not `calc_eigen`.
  All-kinetic heating flags are accepted but produce no heating product. Tests
  confirm unchanged five-column maps and 6+6*S root files (18 for S=2), without
  `.mode` heating rows. Mixed heating flags currently fail preflight. A future
  implementation must explicitly postprocess accepted roots or state that map
  heating remains unavailable; it must not silently append fields to map files.
- FPC consumers call the old diagnostics and apply further output thresholds.
  They must not be enabled just because mixed scalar heating becomes available.
  No FPC routines, arrays, thresholds or output are changed or newly validated.

Future cold-species heating blocks should be assigned NaN in **all eight**
existing per-species heating values (one total plus seven channels). Warm
species retain their calculated values in original order. Unavailable and
numerically failed warm species must be distinguishable in metadata/statuses;
neither should overwrite valid E/B, moments or other species. Total kinetic
heating is unavailable for an all-cold input. Formal wave energy is a separate
API/status and may be meaningful there, notably stationary weakly damped modes.

## 9. Exact follow-up implementation sequence

The next implementation should be divided into these independently reviewable
steps. Nothing in this list is activated by this audit.

1. **Raw calculator and state isolation (support only).** Add a small explicit-
   state heating postprocessor receiving a successful full-plasma E/B, native
   kappa/a, frequency, model/status masks, and snapshots of each chi and selected
   harmonics at one real frequency. Compute every existing raw total and all
   seven raw channels before normalization. Use the existing formulas and
   ordering. Never enter `solver_*` thermal variables for a cold reference.
   Protect/restore mode-frequency `susc`, `susc_low`, `lam`, and dispersion status
   around derivative sampling; do not leave root continuation with shifted data.
2. **Expose missing warm harmonic ingredients (support only).** Extend the
   separate auxiliary kinetic interface to return optional low tensors using
   the shared includes, and populate the already existing exact-parallel limits
   n=0 and +/-1. Keep cold low tensors unavailable. Fully initialize scratch in
   the new interface and preserve species ordering. Do not create a second
   kinetic physics implementation, alter harmonic/cache policy, or introduce
   a kperp/kpar floor. The compiler-sensitive historical fixed128 unit stays
   on its existing compatibility path.
3. **Qualified full wave energy (mathematical correction).** Include vacuum
   |E|^2, B at the chosen real frequency, and derivatives of **every** cold and
   kinetic susceptibility at fixed physical parameters. Start with centered
   differences and an explicit step-convergence check on both matrices and
   E†D E; use analytical cold derivatives where practical. An analytic kinetic
   derivative may be a separate numerical task if no finite-difference plateau
   is resolved. Do not divide by zero/unresolved H, wr=0, a real cold pole, or
   nonfinite results. Return unavailable with a reason, without invalidating a
   previously accepted mode. Negative H is not itself a numerical error.
   Strongly complex/degenerate modes need applicability metadata; small gamma/wr
   alone is not proof of the weak-mode approximation near a branch merger.
4. **Make compatibility explicit before changing old results.** Preserve the
   original all-kinetic `calc_eigen` execution for exact historical replay.
   Correcting its missing vacuum term or snapshot frequency changes Ps even
   with identical roots. It is impossible to silently make both corrections
   and promise byte-for-byte old heating. Recommended architecture: an explicit
   energy/diagnostic convention at the new API, with historical and complete-
   dielectric choices; leave existing all-kinetic default replay unchanged
   pending a deliberate input/default migration decision. Mixed restored
   physical diagnostics should use the qualified complete denominator and
   state that convention in metadata. An eventual public setting's name/default
   needs approval in the implementation specification; no new setting exists
   here. Z and Bessel choices must remain independent of that decision.
5. **Enable collisionless mixed diagnostics (support only).** Remove mixed
   heating guards only after the new backend is validated. For warm species
   calculate existing totals/channels with the full-plasma field; for cold
   species write NaN/status, not zero. Keep fully cold kinetic heating
   unavailable even when its separately requested wave energy is valid. Serial
   and 2D scans reuse the existing continuation and output rows. Preserve
   valid fields/moments on a species heating failure. No temperature-transition,
   radial, parameter, perpendicular-scan or cold-collision work is included.
6. **Output and reader corrections (separate compatibility review).** Correct
   the documented/reader channel stride and heating-only start index without
   changing Fortran column order or precision. Return/assign energy explicitly
   if writing an energy field; distinguish it from the old unassigned Ew.
   Repair low_n=false record reversion only as an acknowledged formatting
   correction, retaining an explicit legacy-replay route if needed. Do not
   insert a wave-energy column into current low_n=true files. An energy sidecar
   or explicit API is preferable. Repeated-map root heating, if enabled, should
   use an optional companion file in the existing diagnostic row format; keep
   `.map` and `.roots` schemas unchanged. Otherwise reject that request clearly.
7. **Existing kinetic corrections requiring separate gates.** Fix shifted
   harmonic snapshots and uninitialized fixed scratch in the corrected path,
   with targeted regression evidence; do not silently force those changes into
   historical replay. The finite-rate off-diagonal/operator finding needs its
   own derivation and tensor/assembly validation, not an isolated heating patch.
   Collisionless mixed restoration can proceed independently; collisional
   heating must not be newly advertised as fully validated until that is done.
8. **FPC remains deferred.** Stage 4B.3B must separately audit FPC normalization,
   operator/current moments and integrated energy transfer. It must not inherit
   unavailable cold values as zeros or a new normalization without validation.

These steps distinguish support plumbing from mathematical changes and parsing
corrections. The operator/closure ambiguity and the public compatibility-policy
decision are recorded rather than resolved by inventing physics or defaults.

## 10. Numerical acceptance criteria for a follow-up implementation

The following are measured-domain gates, not universal guarantees from a
Bessel-tail tolerance or a small determinant residual:

- Preserve exact historical point and full-file regressions under the explicit
  historical diagnostic path, legacy Z and fixed128 settings. Do not replace
  expected values with corrected outputs.
- Check every species raw total and all seven channel numerators against an
  independent implementation, with representative componentwise rtol 3e-8 /
  atol 1e-8 or tighter where conditioning permits. The observed signed-harmonic
  absorptive-matrix scaled discrepancy is 8.27e-15. Raw modal Q agrees to
  7.27e-16 scaled in the high-precision cases.
- At one selected real frequency, require first-four-channel closure to n=0
  within rounding tolerance, but test/record the higher-harmonic and excluded
  longitudinal remainder rather than asserting false total closure. The observed
  normalized closure error is 6.94e-18. Keep low-array initialization tests.
- Validate the complete stationary-cold field/fluid energy identity to about
  5e-12 relative in well-conditioned points. Include the moderate-a vacuum-term
  test, where ignoring |E|^2 causes roughly 50% error.
- Resolve both the derivative matrix and its field contraction. A candidate
  corrected energy should agree with high precision to about 1e-6 relative in
  the ordinary benchmark set and show a stable step window. If cancellation
  prevents this, report energy/normalized heating unavailable; do not silently
  select an arbitrarily tiny step or set H to a finite floor.
- Check weak-mode growth/damping independently, with a benchmark acceptance
  around 2e-4 relative in gamma here; retain the observed 6.97e-5 weak-growth
  and <=4.32e-6 weak-damping errors as reference measurements. Use absolute
  errors when damping is at roundoff. Do not enforce this relation for the
  reactive cases or near an unresolved/degenerate energy denominator.
- Preserve physical rho_ref/d_ref equivalence and auxiliary invariance using
  physical parameters and fixed reference definitions. The direct normalization
  comparison is bit-identical for Ps/Ps_split in this run across both Z backends
  and fixed/adaptive modes. That does not validate the historical Z evaluator's
  absolute physical accuracy near its known defects.
- Require original species and signed-harmonic order, NaNs/status for every
  unavailable cold heating slot, explicit failure of incomplete warm totals,
  and successful subsequent evaluations after diagnostic failures. Preserve
  valid fields, roots and moments on a heating-only failure.
- Test reader columns with unique markers as well as actual files. Exact marker
  equality must hold at the schema above. Do not accept small numerical channels
  masking a stride error. Keep single-point, serial, 2D and map behaviors explicit.
- Keep collisionless mixed availability separate from the additional finite-rate
  operator issue; new collisional accuracy claims require all nine response
  entries, work, the appropriate continuity/reservoir model and wave assembly
  to agree, not just zz or a stationary isotropic test.

## 11. Reproducibility, tests, limitations and explicit failures

The pre-audit standard build and full **169-test** suite passed (108.249 s).
With the optional arbitrary-precision reference importable, full discovery
passed **183 tests** (123.637 s), including the exact historical compatibility
regressions and all 14 new audit tests. Standard-environment discovery also
passed with one optional high-precision test skipped: **182 passed, 1 skipped**
(183 discovered, 122.226 s). The remaining new tests use independent SciPy
fourth-order real differences and analytically differentiated cold linear
solves when mpmath is unavailable. No production package/runtime dependency
is added.

The measurement command runs the 14 audit tests and writes only observations,
not golden data. Production source hashes are recorded before the probes and
verified after them. Test probes use the normal production objects and a
separate complete checked build with O0, `-fcheck=all`, and
`-ffpe-trap=invalid,zero,overflow`. The Ew experiment recompiles unchanged
`disprels.f90` in a temporary directory with `-finit-real=nan`; it does not edit
or replace the ordinary executable. No checked-build numerical failure was
concealed with a floor, damping, or offset.

Toolchain: GNU Fortran 11.4.0 on Linux x86_64; NumPy 2.2.6, SciPy 1.15.3,
and optional mpmath 1.4.1. Production remains in its existing compiler-promoted
binary64 real/complex kinds. Other compilers were not executed. Exact build
flags remain in the unchanged Makefile and the isolated harness. The standard
Python loader emits existing non-fatal unclosed-file ResourceWarnings during
its characterization; no file-lifecycle refactor is included here.

```bash
make -j2
MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v
# Optional arbitrary precision in any importable test environment:
PYTHONPATH=/tmp/plume-z-python MPLCONFIGDIR=/tmp/plume-mpl OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s tests -p 'test_*.py' -v
PYTHONPATH=/tmp/plume-z-python MPLCONFIGDIR=/tmp/plume-mpl OPENBLAS_NUM_THREADS=1 \
  python3 tests/run_mixed_heating_energy_audit.py --require-high-precision \
  --output /tmp/MIXED_HEATING_ENERGY_AUDIT_RESULTS.json
```

`/tmp/plume-z-python` is the already available optional package location used
for this run, not a hardcoded production or test import. If it does not exist,
make mpmath importable through a test environment instead. Without it, omit
`--require-high-precision`; the artifact explicitly labels the reference method.

Additional failures/uncertainties are not hidden:

- The initial full-agreement assertion for finite-rate anisotropic/drifting
  work failed; the component/grid/control investigation in section 7 establishes
  the mismatch and its limits. The final tests are characterization tests, as
  appropriate for an audit, not a claim that the existing collisional tensor
  passes that physical accuracy requirement.
- An exploratory high-frequency **all-warm** trial with tau=0.01 for each species,
  masses mu=1 and 4, charges +/-1, a=0.3, kpar=0.2, kperp=0.1 did not refine
  from the supplied 8+0.01i guess in production. The independent calculation
  located a nearby real frequency about 7.461253, with sampled normalized determinant
  residuals of order 1e-11 in binary64 refinements (not a proven error floor). This trial is recorded in JSON
  and is not counted as a validated weak-mode benchmark. Its cause was not fully
  isolated; no root or susceptibility change is inferred from it. The
  independently validated fully cold moderate-a benchmark is distinct.
- A strict-zero assertion in the optional binary64 cold reference was replaced
  by a machine-roundoff bound: an independent linear solve can leave a tiny
  anti-Hermitian residue. The positive cold energy identity still agrees to
  rounding accuracy. No production or earlier regression tolerance changed.
- One interleaved benchmark has nonzero equilibrium charge/current. Existing
  warning-only behavior is preserved; physical interpretation assumes any
  required external background/current support. Its algebraic response/root
  comparison does not establish the energy balance of an unspecified external
  reservoir.
- Real-velocity quadrature is used only where Im(omega+i*nu)>0; it is not used
  as an incorrect real-contour oracle for damped collisionless Landau response.
  The latter uses independent analytically continued susceptibility. The
  arbitrary-precision energy references sample real frequency.
- A real projection at an exact cold pole, wr=0, a vanishing/resolution-limited
  energy denominator, strong complex mode or branch degeneracy needs an explicit
  applicability/failure policy. No universal exact conserved energy or entropy
  heating is inferred for those cases.

No existing regression failed. The audit ends at diagnosis and an implementation
plan. In particular, it does not enable mixed heating or implement any FPC.

## 12. Complete added/modified-file inventory

Only these **seven new files** are added by this audit:

| File | Role |
| --- | --- |
| `tests/heating_energy_probe.f90` | Direct production snapshots, historical heating calls, controlled scratch poisoning and isolated recovery of existing auxiliary low harmonics. |
| `tests/heating_energy_reference.py` | Independent quadratic forms/channels, signed Stix harmonics, high-precision analytic derivatives and optional independent fallback. |
| `tests/heating_energy_audit.py` | Read-only normal/checked harness, temporary Ew canary build, derivative sweeps and independent root/energy comparisons. |
| `tests/test_heating_energy_audit.py` | Fourteen deterministic numerical/structural characterization tests. |
| `tests/run_mixed_heating_energy_audit.py` | Reproducible measurement runner; no expected-value updates. |
| `tests/MIXED_HEATING_ENERGY_AUDIT_RESULTS.json` | Machine-readable toolchain, physical cases, complex roots/fields, energies, derivative sweeps, channel/reader/initialization/collision findings, exploratory failure and source hashes. |
| `tests/MIXED_HEATING_ENERGY_AUDIT.md` | This audit, equations, evidence, limitations and follow-up plan. |

No pre-existing source, test, documentation or tracked input file is modified
by this stage; normal build products are regenerated by the test commands. All task-start production source,
existing tests, reports, measurements, documentation, build configuration,
plotting/wrapper code, and user changes are preserved. In particular the staged
Stage 4B.2 changes are left as they were. No commits or pushes are made.

Prerequisite conventions and the unchanged moment/field implementations are
recorded in [the Stage 4B.2 report](COLD_FLUID_EIGENFUNCTION_REPORT.md),
[the field report](MIXED_EIGENFUNCTION_REPORT.md), and
[the susceptibility design](COLD_SUSCEPTIBILITY_DESIGN.md). This report supersedes
none of their historical numerical evidence; it provides the audit boundary for
the next heating implementation specification.
