<a id="kinetic-absorptionemission-and-electromagnetic-wave-energy"></a>

# Kinetic absorption/emission and electromagnetic wave energy

[Manual home](./) · [Eigenfunctions](eigenfunctions.md) · [Energy file schema](output.md#energy-companion)

<a id="what-is-calculated"></a>

## What is calculated

`heating_convention='corrected'` is the default, independently of Z backend,
Bessel mode and input normalization. It supports **collisionless** kinetic,
cold and mixed plasmas, including cold references, at single-wavevector diagnostic
calls and serial/2D wavevector scans. Each warm population retains its full species
absorption/emission diagnostic and the original seven selected channels. Every
physical species, including cold ones, contributes to formal wave energy.

These are weak-mode linear response diagnostics, not dimensional thermodynamic
heating power, entropy production or an exact conserved energy for arbitrary
unstable waves. Cold-fluid reversible energy exchange is not kinetic heating.
Cold species' Ps and seven channels are NaN; a fully cold run has no available
kinetic-heating total. This differs physically from a calculated zero for a
warm population at a nonabsorbing point.

<a id="the-accepted-mode-and-real-frequency-snapshots"></a>

## The accepted mode and real-frequency snapshots

Use the normalized mode \(w=w_r+i\gamma/\Omega_R\),
\(\boldsymbol\kappa=\mathbf k d_R\), \(a=v_{A,R}/c\), and
\(\mathbf E=\delta\mathbf E/\delta E_x\), Ex=1, from the full-plasma accepted
**complex-frequency** eigenfunction. In this section w_r means omega_r/Omega_R;
gamma in formulas involving w is likewise normalized. The susceptibility
snapshots below are evaluated at **real laboratory w_r**, with the same fixed
physical wavevector and plasma. We do not replace that laboratory frequency by
a species Doppler frequency for the energy derivative.

Define
\(\chi_{H,s}=(\chi_s+\chi_s^\dagger)/2\) and
\(\chi_{A,s}=(\chi_s-\chi_s^\dagger)/(2i)\). The dagger is complex conjugate
transpose; this is distinct from the ordinary tensor transpose symmetries.
The unnormalized absorption/emission quadratic form is

\[
 Q_{\mathrm{abs},s}=\mathbf E^\dagger\chi_{A,s}(w_r,\mathbf k)\mathbf E
 =\operatorname{Im}[\mathbf E^\dagger\chi_s(w_r,\mathbf k)\mathbf E].
\]

The companion column Q_sN denotes this numerator for species N, **not** the
signed input charge ratio Qs. It is evaluated locally without overwriting the
accepted-mode global susceptibility or wave matrix.

<a id="complete-formal-wave-energy"></a>

## Complete formal wave energy

For the conventional real-frequency homogeneous-plasma weak-mode expression,

\[
 W=\frac{|\delta E_x|^2}{16\pi}\,H,\qquad
 H=H_E+H_B+\sum_s H_s,
\]
\[
 H_E=\mathbf E^\dagger\mathbf E,\quad
 \mathbf B_r=\frac{\boldsymbol\kappa\times\mathbf E}{a w_r},\quad
 H_B=\mathbf B_r^\dagger\mathbf B_r,\quad
 H_s=\mathbf E^\dagger
       \frac{\partial[w_r\chi_{H,s}(w_r,\mathbf k)]}{\partial w_r}\mathbf E.
\]

W is in Gaussian physical field-energy units for the chosen physical deltaEx.
H and its terms are dimensionless. H_E is the **vacuum electric-field** term
from differentiating the identity in epsilon. H_B uses real w_r even when the
mode has small nonzero gamma; it is not obtained by blindly reusing the magnetic
field evaluated with complex w. H_s includes the full species dielectric
response, drifts and anisotropy. Every cold/kinetic species contributes once.
No auxiliary normalization scale is a population or energy term.

The derivative fixes k, species densities/temperatures/anisotropies/drifts,
background B0 and reference scales. Dimensionless differentiation introduces no
extra Omega_R factor because d(omega chi)/domega=d(w chi)/dw. The expression
is conventional for weakly damped/growing isolated modes; its continuation as a
formal quadratic form is not an exact strong-mode conservation law.

<a id="signed-fractions-and-the-electric-magnitude-measure"></a>

For stationary collisionless cold species at real frequency away from poles,
the same expression has an independent bulk-energy identity:

\[
 H_s=\frac{D_s}{\mu_s a^2}\,\mathbf v_s^\dagger\mathbf v_s,
 \qquad
 W_{\mathrm{bulk},s}=\frac{n_{0s}m_s}{4}|\delta\mathbf U_s|^2.
\]

Here v_s is the normalized fluid velocity defined on [Eigenfunctions](eigenfunctions.md).
The second expression is time-averaged oscillatory bulk kinetic energy for
complex-amplitude perturbations; multiplying H_s by abs(deltaEx)^2/(16 pi)
gives it. This is reversible cold bulk energy, not kinetic heating. The simple
positive identity assumes zero equilibrium drift; drifting formal energy need
not be positive and must use the full dielectric derivative.

### Signed fractions and the electric magnitude measure

The companion exposes H_total, H_E, H_B, H_species_total and each H_s, together
with f_E=H_E/H, f_B=H_B/H, f_species=sum(H_s)/H and f_s=H_s/H. These are **signed
fractions of the formal quadratic form**. Negative-energy/drifting modes may
have negative terms or fractions outside [0,1]. They are not clamped. When H
is resolved, f_E+f_B+f_species=1 within numerical precision.

A separate noncancelling contribution measure is

\[
 f_{E,\mathrm{abs}}=
 \frac{H_E}{|H_E|+|H_B|+\sum_s|H_s|}.
\]

It measures electric contribution relative to individual term magnitudes; it
is **not** a physical energy fraction. It is useful when signed terms nearly
cancel. A zero/nonfinite denominator is unavailable. Negative H is not itself
failure; unresolved H or near-zero H relative to its numerical uncertainty
makes signed fractions and normalized Ps unavailable, without a small floor.

<a id="species-normalization-and-energy-consistency"></a>

## Species normalization and energy consistency

The corrected dimensionless species diagnostic is

\[
 P_s=\frac{Q_{\mathrm{abs},s}}{H}.
\]

There is no extra omega, wave period, 2pi or factor of two. For a weak isolated
mode, mutually consistent field/particle energy transfer gives approximately

\[
 \frac{\gamma}{\omega_r}\simeq-\sum_{s\ \mathrm{kinetic}}P_s.
\]

Cold collisionless real-frequency response is Hermitian away from poles; cold
species supply reversible energy in H, not a kinetic absorption channel in that
sum. This approximate relation is independently **checked**, not imposed to
create closure. Reactive/degenerate instabilities, strong damping, branch
interference or unresolved responses can invalidate the interpretation.
A good energy-balance residual alone does not prove applicability.

<a id="the-seven-selected-channels"></a>

## The seven selected channels

Each numerator uses the **same** w_r snapshots and the accepted mode's E, then
is divided by the **same H** as the full species Ps. Here chi_ij,n denotes the
single signed harmonic susceptibility at real w_r, E_perp=(Ex,Ey,0), and
conjugate field factors are explicit.

| Channel | Unnormalized numerator | Historical interpretation |
|---|---|---|
| 1 | abs(Ey)^2 Im chi_yy,0 | n=0 TTD |
| 2 | Im[conj(Ey) Ez chi_yz,0] | n=0 TTD cross |
| 3 | Im[conj(Ez) Ey chi_zy,0] | n=0 LD cross |
| 4 | abs(Ez)^2 Im chi_zz,0 | n=0 LD |
| 5 | Im[E-dagger chi_0 E] | Full n=0 |
| 6 | Im[E_perp-dagger chi_+1 E_perp] | Transverse n=+1 cyclotron |
| 7 | Im[E_perp-dagger chi_-1 E_perp] | Transverse n=-1 cyclotron |

TTD denotes transit-time damping and LD Landau damping in this diagnostic
naming; the signed quantities can also describe emission. At the corresponding
n=0 tensor structure, channels 1–4 sum to channel 5. Channels 5–7 do **not**
generally sum to species-total Ps: higher harmonics and +/-1 longitudinal/cross
components are omitted intentionally. Summing all seven double-counts n=0.
The full species total remains separate. These selected diagnostics are not
redefined to force total closure.

Corrected snapshots initialize every harmonic tensor entry, including the five
unused n=0 entries, and retain exact-parallel analytic n=0,+/-1 limits where the
mixed/cold dispersion supports them. There are no kperp/kpar floors. Derivative
calls cannot overwrite the central total/harmonic snapshots. `low_n` controls
whether seven channels are written in the old mode block; the corrected companion
still contains all seven. FPC harmonic evaluation is a different calculation.

<a id="derivatives-convergence-and-isolation"></a>

## Derivatives, convergence and isolation

Kinetic energy derivatives use centered differences of w_r chi_H and Richardson
extrapolation. For a species let F(w)=w chi_H(w), fixed k/plasma, and h>0:

\[
 C(h)=\frac{F(w_r+h)-F(w_r-h)}{2h},\qquad
 R(h)=\frac{4C(h)-C(2h)}{3}.
\]

The implementation tries 13 halving relative steps from .064 to .000015625,
requiring successive estimates to resolve both matrix Frobenius norm and the
E-contracted scalar at tolerance 1e-6, with comparison scales max(1,norm/value).
It selects the best stable window, not necessarily the smallest step. Exported
dmatrix_s, dproject_s and dstep_s describe that window; uncertainty is an
empirical numerical resolution estimate, not a rigorous error bound. Poles,
Bessel/Z errors and cancellation can defeat convergence. No fallback silently
substitutes the historical normalizer.

Cold derivatives use the analytic differentiation of the normalized fluid
momentum/continuity response. Exact real-frequency poles remain failures.
The finite-difference sampling of kinetic species does not unnecessarily evaluate
cold responses at off-center frequencies. Corrected diagnostic work is local:
accepted-mode susc, susc_low, lam and dispersion status remain at the converged
mode. Failed energy/heating does not invalidate root, field or species moments.

<a id="availability-applicability-and-statuses"></a>

## Availability, applicability and statuses

Corrected diagnostics require collision_type=0. They reject finite-rate physics
rather than making new collisional claims. Legacy all-kinetic collisional
behavior can be requested explicitly, with [known limitations](limitations.md).

Failures include w_r=0, real cold poles, nonfinite response, unresolved derivative,
unresolved H and unavailable Ex=1 field. A valid complex dispersion root may
therefore have unavailable real-frequency energy/heating. A cold singularity at
w_r is not concealed by the mode's finite imaginary part.

Applicability bit 1 is set for abs(gamma/w_r)>.05. Bit 2 flags weak-balance or
near-reactive inconsistency. Formal energy can remain resolved (energy_status=0)
while an applicability bit makes normalized warm Ps/channels NaN. Raw Q and H
terms remain informative but must not be presented as validated thermodynamic
heating. A reactive mode is not repaired by dividing by an artificial H floor.
Numerical energy status and physical applicability are separate fields; the
[output status tables](output.md#energy-companion) list both.

<a id="running-and-reading-existing-examples"></a>

## Running and reading existing examples

```bash
mkdir -p data/example
./plume.e inputs/example/heating_corrected_kinetic.in
./plume.e inputs/example/heating_corrected_warm_reference.in
./plume.e inputs/example/heating_corrected_cold_reference.in
./plume.e inputs/example/heating_corrected_cold.in
```

| Existing input | First printed mode frequency | Formal H near first point |
|---|---|---:|
| heating_corrected_kinetic | .06098715 - 1.270997e-5 i | 20342.68 |
| heating_corrected_warm_reference | .06153842 - 5.564326e-6 i | 20384.80 |
| heating_corrected_cold_reference | .05916735 - 1.588769e-6 i | 20516.88 |
| heating_corrected_cold | .05972252 plus roundoff imaginary part | 20607.63 |

These bounded inputs and numbers are verified against current output and the
independent heating tests; they are not universal branch frequencies. Inspect
statuses as well as values. For a separate moderate-a test with artificial mass
ratio 4, the independent suite finds H_E/H about .500000281. That is a test
benchmark, **not** the physical mass ratio in the supplied cold example; the
curated public moderate-a input is still outstanding.

```python
from linfpclib.linfpc import load_plume_energy
r = load_plume_energy('data/example/heating_corrected_cold_reference.energy')
print(r['H_total'], r['H_E'], r['f_E'], r['f_E_abs'])
print(r['energy_status'], r['applicability'])
print(r['H_s1'], r['Ps_s1'])  # cold proton energy exists; kinetic Ps is NaN
```

New energy columns are in the versioned `.energy` companion, not inserted in
the old low_n=true heating blocks or in the old unassigned `Ew` slot. Output
association and all status fields are specified in [Output](output.md).

<a id="legacy-convention"></a>

## Legacy convention

`heating_convention='legacy'` preserves historical all-kinetic calc_eigen
arithmetic, old normalizer (without vacuum electric energy), one-sided derivative,
shifted harmonic snapshots and legacy scratch/output behavior. It does not
silently inherit corrected mathematics, and cold/mixed legacy heating is not
available. It is for reproducibility, not a second accurate complete energy
convention. See [Historical reproduction](reproducibility.md) for all three
independent switches and the separate role of a source revision.

Scientific context: [PLUME (Klein, Howes & Brown, 2025)](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
and [Huang et al. (2024)](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/velocityspace-signature-of-transittime-damping/60F34E118AFFC3701071477C52A81174).
Independent mathematical checks and historical implementation choices are linked
from the [developer/history guide](development.md), not substituted for this
current user definition.
