<a id="numerical-methods-and-scientific-checks"></a>

# Numerical methods and scientific checks

[Manual home](index.md) · [Normalization](normalization.md) · [Limitations](limitations.md)

This page describes the numerical algorithms. Source constants such as the
secant tolerance, derivative ladder and FPC harmonic range are not public
namelist fields. Record the source revision and all selectable settings with
each result.

<a id="tensors-and-the-determinant"></a>

## Tensors and the determinant

For Fourier convention exp(i k.r-i omega t), B0 along z and k=(kperp,0,kpar),
PLUME stores the ordinary dimensionless per-population tensor
`susc(species,row,column)=chi_s`. With w=omega/Omega_R,
kappa=k d_R and a=v_A,R/c, it forms

\[
 \epsilon=I+\sum_s\chi_s,\qquad
 \mathbf n=\frac{\boldsymbol\kappa}{a w},\qquad
 \Lambda=\epsilon+\mathbf n\mathbf n^T-n^2I,\qquad
 \det\Lambda=0.
\]

The tensor n n^T has ordinary transpose, not conjugation; complex n is allowed.
Vacuum identity and electromagnetic terms appear exactly once, outside species
kernels. Tensor component order is xx,xy,xz; yx,yy,yz; zx,zy,zz. The gyrotropic
response has chi_yx=-chi_xy, chi_zx=chi_xz, chi_zy=-chi_yz, **without conjugation**
at complex frequencies. Hermitian symmetry at a real nondissipative point is
a separate physical statement.

The API `disp` and map files use the raw determinant.
Direct stdout historically multiplies it by solver_vtp^6 for a warm reference,
whereas cold-reference direct evaluation reports it raw. This display scaling
does not redefine chi, epsilon or Lambda. Historical comments about multiplying
expressions by w^2 do not describe the implemented matrix.

A small determinant can result from scaling or cancellation. Check the
normalized determinant, matrix/nullvector residual, finiteness and root stability
under controlled resolution changes. A small raw determinant alone does not
establish frequency accuracy.

<a id="kinetic-response-and-normalization-adapter"></a>

## Kinetic response and normalization adapter

The thermal susceptibility follows Stix, *Waves in Plasmas*, Chapter 10,
equations 66–73, for the [bi-Maxwellian model](models.md#kinetic-populations).
Warm-reference d_ref inputs are converted to the historical thermal coordinates
with the reference anisotropy factor. Cold-reference mixtures use the independent
auxiliary state, with physical species anisotropies untouched. Physical
prefactors, Doppler shifts and thermal/Bessel arguments are invariant under that
coordinate choice. The adapter reuses the shared kinetic expressions.
Fully cold response bypasses thermal state.

<a id="adaptive-bessel-harmonics-and-caching"></a>

## Adaptive Bessel harmonics and caching

For a warm species define the scaled modified Bessel coefficient
\(\Gamma_n(\lambda)=e^{-\lambda}I_{|n|}(\lambda)\), where

\[
 \lambda_s=\frac{\kappa_\perp^2 Q_s^2\alpha_s\theta_s}{2\mu_s}
 =\frac{K_\perp^2 Q_s^2\alpha_s}
        {2\mu_s\tau_s^{(\rho)}\alpha_R}.
\]

theta_s is d_ref tauS; the second form applies to a thermal reference. This is
one physical argument, independent of auxiliary scale. Cold species have no
Bessel arrays. Harmonic count N means **maximum absolute signed order**: -N..N
contains 2N+1 terms, with additional adjacent orders stored for combinations.

`bessel_mode='adaptive'` is the default. The selector bounds both the omitted
ordinary weight and the normalized second moment:

\[
 E_0=2\sum_{n>N}\Gamma_n,\qquad
 E_2=\frac{2}{\lambda}\sum_{n>N}n^2\Gamma_n\quad(\lambda>0).
\]

The selector uses conservative logarithmic bounds from the generating function
\(M(t)=\exp[\lambda(\cosh t-1)]\). For m=N+1 and
\(t=\operatorname{asinh}(m/\lambda)>0\),

\[
 E_0\le B_0=2e^{\lambda(\cosh t-1)-mt},\qquad
 E_2\le B_0(\cosh t+\lambda\sinh^2t).
\]

It chooses a species-specific N satisfying both bounds at bessel_tol (default
1e-10), subject to bessel_nmax (default 512). Bounds are evaluated stably in
logarithmic form, avoiding subtraction of nearly converged sums from one and
premature termination for broad distributions.
For lambda=0 the selector handles the weights analytically; supported exact
parallel response uses its analytic harmonic limit. This does **not** enable the
currently unsupported all-kinetic kperp=0 dispatcher.

The tail tolerance does not guarantee total susceptibility or root accuracy. Resonance
amplification, cancellation, Bessel coefficient accuracy and Z accuracy also matter.
Compare all relevant tensor elements and mode frequency while tightening
bessel_tol and increasing the cap, especially for broad lambda distributions.

If the cap is reached before the bounds meet tolerance, PLUME warns once for the
unchanged request and continues with the maximum allowed N. The .bessel metadata
records the convergence as **unresolved**.
A nonfinite Bessel calculation is numerical failure, distinct from a finite but
under-resolved sum. There is no silent physical substitution or damping.

The per-species cache tracks physical lambda, request, active N,
allocation capacity, coefficients, tail bounds and convergence. Frequency,
parallel wavevector, drift and density alone do not require new physical Bessel
coefficients. Perpendicular wavevector, temperature, anisotropy, mass or charge
changes that alter lambda do. Exact computed-argument matching avoids unbounded
approximate reuse; roundoff in forming an equivalent argument can cause a miss,
not reuse of incorrect data. Allocation grows as needed without shrinking at
every small step. Species-count/order changes cannot reuse the wrong argument.

`bessel_mode='fixed'` uses bessel_nfixed; `'fixed128'` requires nfixed=128.
The historical fixed-128 compiled arithmetic and its old cache path are retained
for exact compatibility. The historical path does not report per-species cache
statistics. Changing the auxiliary scale may require recomputation. Independent
choices of Z/heating conventions are required for
[legacy reproduction](reproducibility.md).

FPC distributions use a separate **internal fixed nbesmax=15**, not this adaptive
susceptibility policy. Neither bessel_tol nor an adaptive root calculation is
evidence that a velocity-space FPC reconstruction is harmonically converged.

<a id="plasma-dispersion-function"></a>

## Plasma dispersion function

For Im zeta>0 the kinetic Z function is

\[
 Z(\zeta)=\frac{1}{\sqrt\pi}\int_{-\infty}^{\infty}
                 \frac{e^{-t^2}}{t-\zeta}\,dt.
\]

The Landau analytic continuation defines it elsewhere. With Faddeeva function
\(\mathcal w(\zeta)=e^{-\zeta^2}\operatorname{erfc}(-i\zeta)\),
\(Z(\zeta)=i\sqrt\pi\,\mathcal w(\zeta)\). This mathcal w is not the normalized
wave frequency w. The signed parallel convention used by the production
`zet_in(kpar,zeta)` interface is Z(zeta) for kpar>0 and -Z(-zeta) for kpar<0.
It is not naive lower-half-plane conjugation.

`z_backend='faddeeva'` (default) uses bundled Steven G. Johnson code with a C99
bridge and Fortran ISO_C_BINDING. The bridge uses upstream evaluation at machine
precision; there is no user Z-tolerance setting, and no Python runtime is required.
Its continued
fraction/Algorithm-916/Taylor branches are numerical techniques, not changes to
the Landau prescription. The historical radius-three approximation error motivated
replacement; the audit isolated its susceptibility amplification independently.

The Faddeeva backend assigns a result on every call. Invalid/nonfinite arguments,
unsupported kpar=0 calls, arithmetic failures and nonrepresentable outputs
propagate invalid status/NaNs through guarded dispersion paths. It does not clamp
Z, add imaginary frequency or fall back silently to the old algorithm. Some
lower-half-plane Z values exceed binary64 range; this is distinct from an
implementation failure. Subsequent regular evaluations must still succeed.
Fully cold dispersion does not call Z.

`z_backend='legacy'` preserves historical zetout arithmetic/stopping behavior,
including known accuracy and unassigned-result limitations; it is not an accurate
fallback. The Faddeeva bridge uses c_double/c_double_complex matching the promoted
Fortran real/complex storage. Under these flags, do not assume `nrtype::DP` is
interchangeable with these C kinds.
[Provenance and license](https://github.com/kgklein/PLUME/blob/main/third_party/faddeeva/README.md)
records the upstream mirror commit, checksums, MIT attribution and unchanged
numerical sources. No external shared third-party library is required.

<a id="drifting-cold-susceptibility"></a>

## Drifting cold susceptibility

For one species let delta=w-kappa_parallel u, b=kappa_perp u,
g=mu/Q signed and A=D mu/(Q^2 a^2), using the physical
[normalization definitions](normalization.md). The species circular and coupling
terms from [Verscharen & Chandran, Section 2.1](https://arxiv.org/pdf/1212.5192) are

\[
 R_s=-\frac{A\delta}{w^2(\delta+g)},\quad
 L_s=-\frac{A\delta}{w^2(\delta-g)},\quad
 J_s=-\frac{Ab\delta}{w^2(\delta^2-g^2)},\quad
 M_s=\frac{iAbg}{w^2(\delta^2-g^2)},
\]
\[
 P_s=-\frac{A}{\delta^2}-\frac{Ab^2}{w^2(\delta^2-g^2)},\quad
 S_s=(R_s+L_s)/2,\quad D_{H,s}=(R_s-L_s)/2,
\]
\[
 \chi_s=\begin{pmatrix}
 S_s&-iD_{H,s}&J_s\\ iD_{H,s}&S_s&M_s\\ J_s&-M_s&P_s
 \end{pmatrix}.
\]

D_H is the Hall coefficient, not the input density ratio D. These are ordinary
species chi terms with **no vacuum identity**. PLUME uses factored denominators
(delta-g)(delta+g) and combined transverse formulas to reduce
subtractive cancellation. Finite electron mass is retained. At drift zero the
cross couplings vanish; parallel stationary transverse roots recover n^2=R or L
for R=1+sum R_s, L=1+sum L_s. Perpendicular stationary O/X relations are n^2=P
and n^2=RL/S with the total dielectric coefficients, away from poles.

Exact w=0, delta=0 and delta=+/-g are explicitly classified for nonzero density;
zero density gives zero susceptibility. A full tensor may be invalid at a pole
even if some individual component has a removable cancellation. No pressure,
thermal correction, collision floor or artificial damping is added.

<a id="root-solving-derivatives-and-numerical-checks"></a>

## Root solving, derivatives and numerical checks

The production root method is a complex secant iteration, with internal
tolerance 1e-13, nearby initial fractional offsets 1e-7 and a maximum of
75 iterations. These constants are not user inputs. See [Modes](modes.md) for
map minima, deduplication/filter order and continuation. Guarded scan points
also check status, finiteness and a determinant divided by the cube of the
maximum matrix magnitude against 1e-8. This is a scaling/acceptance safeguard,
not a universal frequency error estimate or a new branch tracker.

Corrected energy derivatives use centered finite differences with Richardson
extrapolation and convergence checks on matrices **and their field contractions**,
together with analytic cold response derivatives. Their physical meaning and
resolution limits are on [Heating](heating.md#derivatives-convergence-and-isolation). Failed supplementary
energy must not change accepted-mode matrix state or root continuation.

For publishable calculations, vary map resolution or supplied guesses, scan
step size, Bessel tolerance/cap and backend. Reconstruct Lambda
from species tensors and check field/current identities where relevant. Compare
forward/reverse traversal with attention to physical mode crossings. Retain
unresolved cap warnings, invalid status and unavailable diagnostic distinctions
when analyzing or plotting results.
