# Cold susceptibility: stages 0 and 1

Historical Stage 0/1 audit and validation. Production integration was added in
[Stage 2](MIXED_DISPERSION_REPORT.md); the cold kernel equations and numerical
contract below remain unchanged.

## Audit completed before implementation

The audit uses the normalization report and input documentation, and the
following production routines (line numbers refer to revision `58f3aad`):

| Location | Convention established by the expressions |
| --- | --- |
| `src/vars.f90`, `specie`, lines 63–103 | `mu_s=m_R/m_s`, signed `Q_s=q_R/q_s`, `D_s=n_s/n_R`, `vv_s=U_s/v_A,R`. |
| `src/vars.f90`, `sync_solver_state`, lines 390–438; `src/functions.f90`, `read_in_params` | Input validation and conversion require positive kinetic temperatures. The new routine must not call them. |
| `src/disprels.f90`, `disp`, lines 3676–3700 | Synchronizes thermal solver state, then selects collisionless or weakly collisional kinetics. No cold dispatch is present or added. |
| `disp0`, lines 4560–4645 and 4858–4937 | Resonance arguments contain `om-kpar*vv_s/sqrt(betap*alpha_R)-n*mu_s/Q_s`. The susceptibility prefactor includes `betap*mu_s*D_s/(om**2*vtp**2*Q_s**2)`. |
| `disp1`, lines 3723–4360 | Uses the same tensor normalization, with collision terms in the kinetic response. |
| `disp0`, lines 4884–4929; `disp1`, lines 4306–4352 | `susc(species,row,column)` stores ordinary susceptibility, not dielectric identity or a frequency-weighted response. Each species is stored separately. |
| `calc_eigen`, lines 2422–2433 | Faraday's law is `B1=c*k cross E1/omega`, consistent with the Fourier convention below. No change to this routine is needed in this stage. |

Cartesian coordinates have `B0` along positive z and `k=(kperp,0,kpar)`.
The susceptibility and signed resonance denominators are consistent with
`exp(i*k.r-i*omega*t)`; positive imaginary frequency means growth. In particular,
the stationary cold limit of the kinetic xy expression is
`-i*A*g/[w*(w**2-g**2)]`, fixing the Hall sign independently of array symmetry.
Species gyrofrequency is signed: `g=mu_s/Q_s` (negative for electrons).
Frequency is `w=omega/Omega_R`, with the usual positive reference gyrofrequency
scale and signed species charges.

The nine stored elements are rows `(xx,xy,xz)`, `(yx,yy,yz)`, `(zx,zy,zz)`.
The existing relationships are `yx=-xy`, `zx=xz`, `zy=-yz`, with no complex
conjugation. The wave matrix is

\[
\Lambda=I+\sum_s\chi_s+\mathbf n\mathbf n^T-n^2 I,
\qquad \mathbf n=c\mathbf k/\omega.
\]

Thus its diagonal electromagnetic terms are `-nz**2`, `-nx**2-nz**2`,
`-nx**2`, and its xz/zx terms include `+nx*nz`. Both routines return
`det(Lambda)`. The historical comment claiming multiplication by `w**2` is
inconsistent with the current executable expressions. Some diagnostics apply
additional output scaling; that does not change `susc` or `lam`.

PLUME's historical gyroradius is perpendicular:
`rho_R=v_tperp,R/abs(Omega_R)`. The adapter uses
`k_rho=k_d*sqrt(beta_R*alpha_R)`, `vtp_solver=a*sqrt(beta_R)` and
`tau_solver=beta_R/tau_d`, preserving drift `u_s=vv_s`. Consequently
`betap/vtp_solver**2=1/a**2`, and the kinetic Doppler shift becomes
`kpar_d*u_s`. Omitting `alpha_R` would change the physical wavevector.

## Independent cold-fluid equations

Use [Verscharen & Chandran (2013), Section 2.1, equations (3)–(10)](
https://arxiv.org/pdf/1212.5192), retaining finite electron mass. Define

\[
a=v_{A,R}/c,\quad w=\omega/\Omega_R,\quad
\delta=w-\kappa_\parallel u_s,\quad b=\kappa_\perp u_s,\quad
g=\mu_s/Q_s,\quad A=D_s\mu_s/(Q_s^2a^2).
\]

The single-species terms (without any vacuum identity) are

\[
R_s=-\frac{A\delta}{w^2(\delta+g)},\quad
L_s=-\frac{A\delta}{w^2(\delta-g)},\quad
J_s=-\frac{Ab\delta}{w^2(\delta^2-g^2)},\quad
M_s=\frac{iAbg}{w^2(\delta^2-g^2)},
\]
\[
P_s=-\frac{A}{\delta^2}-\frac{Ab^2}{w^2(\delta^2-g^2)},\qquad
\chi_s=\begin{pmatrix} S_s&-iH_s&J_s\\iH_s&S_s&M_s\\J_s&-M_s&P_s\end{pmatrix},
\quad S_s=(R_s+L_s)/2,\ H_s=(R_s-L_s)/2.
\]

`H_s` denotes the Hall coefficient, distinct from the input density `D_s`.
There is no transpose, conjugation, or sign conversion to PLUME. For numerical
evaluation, factor `F=(delta-g)*(delta+g)` and use
`S_s=-A*delta**2/(w**2*F)`, `H_s=A*delta*g/(w**2*F)` instead of subtracting
nearly equal circular responses. Split `P_s` as above to cancel its redundant
`w**2` factors analytically. Keep all frequency arithmetic complex.

[NHDS susceptibility](https://github.com/danielver02/NHDS/blob/master/src/suscep.f90)
and [wave matrix](https://github.com/danielver02/NHDS/blob/master/src/disptensor.f90)
return/assemble the scaled species response `a**2*w**2*chi_s`. ALPS's
`ALPS_NHDS::calc_chi_cold` does likewise on its `kperp_norm` path;
`ALPS_fns` divides by `om**2*vA**2` to form `chi0`. These scaled quantities
must not be substituted directly into PLUME. The
[NHDS license](https://github.com/danielver02/NHDS/blob/master/LICENSE) and
[ALPS license](https://github.com/danielver02/ALPS/blob/master/LICENSE) were
reviewed: BSD-style redistribution requires retaining notices. This kernel is
an independent implementation of the analytical equations; no implementation
code is copied from either project.

## Interface and numerical contract

```fortran
subroutine calc_chi_cold(species, omega, kpar_d, kperp_d, va_ref_c, chi, ierr)
```

The new `cold_susceptibility` module accepts `type(specie)` and explicit
inertial-normalized inputs, and returns a full complex `chi(3,3)`. Only `mu_s`,
`Q_s`, `D_s`, and `vv_s` are read. Temperature, anisotropy, collision rate,
global beta, global `vtp`, species arrays and normalization mode are unused.
The calculation works with zero reference beta and no allocated global species
state. It neither enforces equilibrium neutrality/current balance for a single
species nor chooses a species model. It is collisionless and nonrelativistic.

Use `nrtype::SP/SPC`, which match PLUME's default real/complex kinds (64-bit
real under `-fdefault-real-8`). `nrtype::DP=kind(1.d0)` becomes 128-bit under
that GNU flag and is not the production arithmetic kind. No precision changes
are made to existing modules.

All used inputs must be finite, with `mu_s>0`, `a>0`, `Q_s!=0`, `D_s>=0`.
After validating these parameters, zero density returns an exactly zero tensor,
including at frequencies that would be poles for a populated species.
The required status is:

| Code | Public constant | Meaning |
| --- | --- | --- |
| 0 | `cold_ok` | Valid response. |
| 1 | `cold_invalid_input` | Nonfinite input or invalid physical parameter. |
| 2 | `cold_zero_frequency` | Populated species at `w=0`. |
| 3 | `cold_zero_doppler` | Populated species at `delta=0`. |
| 4 | `cold_cyclotron_pole` | Populated species at `delta=+g` or `-g`. |

Checks occur in this order; all failure outputs have quiet NaN real and imaginary
parts. Callers must inspect `ierr`. No arbitrary neighborhood of a pole is
rejected and no floor, damping or collision frequency is added. Finite inputs
must also be in the representable arithmetic range; this is not arbitrary
precision evaluation.

At `w=0` the ordinary susceptibility generally has poles. For a stationary
species the transverse diagonal has a removable cancellation, but the Hall and
longitudinal components still diverge. At nonzero `w` and `delta=0`, transverse
components and drift couplings can remain finite or vanish, while `P_s` has a
second-order pole. At `delta=+/-g`, cyclotron poles remain in the full tensor,
although `b=0` removes oblique couplings. A full tensor at these exact poles is
therefore not returned as a finite response. Nearby complex frequencies are
evaluated directly. Hermitian symmetry holds at nonsingular real frequencies,
not generally at complex frequencies; transpose identities above always hold.

GNU's production `-ffast-math` can optimize away IEEE finite tests. Compile only
the new kernel with a trailing `-fno-fast-math` (Intel: `-fp-model precise`) to
preserve its validation contract. Existing kinetic flags and code stay intact.

## Independent validation matrix

A direct-call Fortran probe links only the new module, `vars`, and `nrtype`.
Python unittest builds both the production configuration and a checked GNU
configuration with bounds checks and floating-point traps. All nine complex
elements are compared with two independent oracles:

1. Literal Python circular-response equations from the paper, forming `R,L,J,M,P`.
2. A NumPy linear solve of the cold momentum and continuity equations. For each
   electric-field basis vector, obtain `B1=kappa cross E/(a*w)` and
   `F=E+a*u*z_hat cross B1`. Solve
   `[-i*delta,-g,0; g,-i*delta,0; 0,0,-i*delta]*h=F` and form
   `chi*E=i*A/w*(h+u*z_hat*(kappa.h)/delta)`. The second term is the drift
   density current; dropping it or the magnetic force loses oblique couplings.

Cases cover stationary species, both drift directions, parallel/perpendicular
and oblique propagation, signed wavevectors, growing/damped and real
frequencies, proton/electron/alpha parameters, multiple mass/density ratios and
Alfven speeds, zero beta/thermal parameters, thermal-field invariance, charge
reversal, density scaling, exact vanishing entries, exact poles and nearby
complex points. Both oracles compare elements, not determinants. Well-conditioned
points use componentwise `atol=5e-13, rtol=5e-13`; near-pole tests use separately
documented tolerances to allow conditioning of the independent linear solve:
imaginary offsets `+/-1e-3`, `+/-1e-4`, `+/-1e-5` from each pole use
`atol=5e-13, rtol=1e-9`. These offsets are test inputs, not kernel regularization.

Run `make -j2` and
`MPLCONFIGDIR=/tmp/plume-mpl python3 -m unittest discover -s tests -p 'test_*.py' -v`.
The preimplementation baseline passes all 15 existing normalization tests.
Existing expected numerical results are not changed.

## Implementation and validation results

With GNU Fortran 11.4.0, `make -j2` succeeds and full unittest discovery passes
**27 tests**: 12 new cold-susceptibility tests and all 15 existing normalization
and kinetic regression tests. Each new test exercises both the production
build (with the kernel's IEEE-preserving override) and the checked build.
Both probes confirm 64-bit real storage and 15 decimal digits of precision.
The Intel build configuration was not executed because `ifort` is unavailable.

Maximum measured `abs(actual-reference)/max(1,abs(reference))`, over every
complex tensor element and both builds:

| Comparison | Maximum scaled error |
| --- | --- |
| Paper equations, well-conditioned cases | `2.062e-13` |
| Independent fluid solve, well-conditioned cases | `8.129e-16` |
| Both oracles, near-pole cases | `3.374e-12` |

All comparisons pass their componentwise tolerances. The paper oracle uses
the literal circular-response subtraction and is less well-conditioned than
the combined expressions for some mass/frequency ratios. The unchanged
normalization suite reports maximum scaled difference `4.668e-13` and maximum
frequency-component difference `6.939e-18`; its pinned legacy probe comparison
remains exact. `git diff --check` passes.

Files added or changed by this task:

| File | Change |
| --- | --- |
| `src/cold_susceptibility.f90` | Independent single-species tensor and explicit status contract. |
| `Makefile` | Build/link the unused module; add module dependencies and kernel-only precision flags. |
| `tests/cold_susceptibility_probe.f90` | Direct-call probe with zero global beta/thermal speed and no allocated species arrays. |
| `tests/test_cold_susceptibility.py` | Two independent oracles and production/checked-build validation. |
| `tests/COLD_SUSCEPTIBILITY_DESIGN.md` | Preimplementation audit, equations, numerical contract and results. |

No existing production Fortran source or existing test/expected value was
modified. Source-call inspection finds the kernel invoked only by the new test
probe. Pre-existing workspace changes, including `plotter/example/cycle.sh`,
are outside this task and were preserved.

## Subsequent integration (excluded here)

Future work must dispatch species before the thermal adapter, handle a cold
reference without fictitious thermal scales, and specify cold density/velocity
eigenfunctions, heating, collision and FPC behavior. Production validation,
root finding, dispatch, kinetic susceptibility, scans and normalization behavior
remain unchanged in these stages. Linking an unused standalone module does not
activate cold species in normal runs.
