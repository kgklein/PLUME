<a id="normalization-and-physical-units"></a>

# Normalization and physical units

[Manual home](./) · [Input reference](input.md) · [Species models](models.md)

This page defines the physical symbols used throughout the manual. Temperatures
are energies (Boltzmann's constant is absorbed into T). Electromagnetic formulas
use Gaussian units. The equilibrium is homogeneous, with
\(\mathbf B_0=B_0\hat{\mathbf z}\) and
\(\mathbf k=(k_\perp,0,k_\parallel)\). Time dependence is
\(\exp(-i\omega t)\); \(\omega=\omega_r+i\gamma\) grows when \(\gamma>0\).

<a id="reference-scales"></a>

## Reference scales

PLUME uses `species_1` as reference R; there is no reference-index input.
Use the usual positive reference gyrofrequency scale (a positive reference
charge and B0), retaining signed charge ratios for other populations. Define

\[
 \Omega_R=\frac{q_R B_0}{m_R c},\qquad
 \omega_{pR}^2=\frac{4\pi n_R q_R^2}{m_R},\qquad
 v_{A,R}=\frac{B_0}{\sqrt{4\pi n_R m_R}},\qquad
 d_R=\frac{c}{\omega_{pR}}=\frac{v_{A,R}}{\Omega_R}.
\]

Here c is light speed, m_R, q_R and n_R set the reference mass, charge and
number-density scales. For the thermal reference define

\[
 v_{t\parallel,R}=\sqrt{\frac{2T_{\parallel,R}}{m_R}},\quad
 \alpha_R=\frac{T_{\perp,R}}{T_{\parallel,R}},\quad
 \rho_R=\frac{v_{t\perp,R}}{|\Omega_R|},\quad
 \beta_R=\frac{8\pi n_R T_{\parallel,R}}{B_0^2}.
\]

The historical gyroradius uses the **perpendicular** thermal speed. Thus
\(v_{t\parallel,R}/v_{A,R}=\sqrt{\beta_R}\) and
\(\rho_R/d_R=\sqrt{\beta_R\alpha_R}\) for a finite-temperature reference.
PLUME's Alfvén speed uses reference density/mass, not the mixture's total mass
density. Do not substitute a total-mixture Alfvén speed without converting
all normalized inputs.

<a id="inputs-in-the-two-conventions"></a>

## Inputs in the two conventions

Use `normalization='rho_ref'` (default) or `'d_ref'` in `&params`. Let
\(w=\omega/\Omega_R\), \(a=v_{A,R}/c\),
\(K_j=k_j\rho_R\), and \(\kappa_j=k_jd_R\), with j parallel or perpendicular.
For a physical population s, use \(\mu_s=m_R/m_s\),
\(Q_s=q_R/q_s\) (signed), \(D_s=n_s/n_R\),
\(u_s=U_s/v_{A,R}\), and \(\alpha_s=T_{\perp,s}/T_{\parallel,s}\).

| Input | `rho_ref` | `d_ref` |
|---|---|---|
| `betap` | Reference parallel beta, beta_R | Same; zero for a cold reference |
| `kperp`, `kpar` | K_perp, K_parallel | kappa_perp, kappa_parallel |
| `vtp` | v_tparallel,R/c | a=v_A,R/c |
| `tauS` | tau_rho,s=T_parallel,R/T_parallel,s | theta_s=2T_parallel,s/(m_R v_A,R^2) |
| `muS` | m_R/m_s | Same |
| `Qs` | q_R/q_s, signed | Same |
| `Ds` | n_s/n_R | Same |
| `alphS` | T_perp,s/T_parallel,s | Same for warm species; ignored for cold |
| `vvS` | u_s=U_s/v_A,R | Same; never thermal-speed-normalized |
| `g_om`, `g_gam` | omega_r/Omega_R, gamma/Omega_R | Same |
| `Kn` | Neutral mean free path/rho_R | Neutral mean free path/d_R |

Use theta_s as a symbol for d_ref `tauS` in this manual to avoid confusing a
physical temperature with the inverse ratio in rho_ref. Both are the same
Fortran input name. Ordinary self ratios for a physical reference population
are mu_1=Q_1=D_1=1. The solver does not comprehensively enforce self-ratio or
reference-density consistency; unusual component decompositions require a
clearly specified reference density. They do not cause automatic renormalization.

For d_ref the physical thermal and species frequency relations are

\[
 \frac{v_{t\parallel,s}^2}{v_{A,R}^2}=\mu_s\theta_s,\quad
 \beta_{\parallel,s}=D_s\theta_s,\quad
 \frac{\Omega_s}{\Omega_R}=\frac{\mu_s}{Q_s},\quad
 \frac{\omega_{ps}^2}{\Omega_R^2}=\frac{D_s\mu_s}{Q_s^2a^2}.
\]

The Doppler-shifted frequency is \(\delta_s=w-\kappa_\parallel u_s\).
The gyrofrequency is signed; the plasma-frequency coefficient uses Q_s squared.
These definitions apply to kinetic and cold species. theta_s=0 means cold fluid,
not a small positive kinetic temperature.

<a id="converting-equivalent-physical-states"></a>

## Converting equivalent physical states

When beta_R>0 and alpha_R>0,

\[
 K_j=\kappa_j\sqrt{\beta_R\alpha_R},\qquad
 \mathrm{vtp}_{\rho}=a\sqrt{\beta_R},\qquad
 \tau_s^{(\rho)}=\frac{\beta_R}{\theta_s}.
\]

Keep mu_s, Q_s, D_s, alpha_s, u_s and w unchanged. In particular, converting
only the wavevector by sqrt(beta_R) is wrong unless alpha_R=1. Matching collision
frequencies additionally requires
\(\mathrm{Kn}_d=\mathrm{Kn}_\rho\sqrt{\beta_R\alpha_R}\).
The existing derived rates are

\[
 \frac{\nu_s}{\Omega_R}=\frac{1}{\sqrt{2}\,\mathrm{Kn}_\rho}
       \sqrt{\frac{\mu_s}{\tau_s^{(\rho)}\alpha_R}}
   =\frac{\sqrt{\mu_s\theta_s}}{\sqrt{2}\,\mathrm{Kn}_d}.
\]

These rates belong only to the restricted all-kinetic Krook implementation;
none of these conversions enables cold collisions.

<a id="worked-equivalent-inputs"></a>

### Worked equivalent inputs

The unchanged
[normalization_rho_ref.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/normalization_rho_ref.in)
and [normalization_d_ref.in](https://raw.githubusercontent.com/kgklein/PLUME/main/inputs/example/normalization_d_ref.in)
represent a proton/electron system with beta_R=.4, alpha_R=1.5,
a=1e-4, alpha_e=1, u_e=.02, and theta_e=.2.

| Quantity | rho_ref input | d_ref input |
|---|---:|---:|
| Initial kperp | .03872983346207417 | .05 |
| Initial kpar | .046475800154489 | .06 |
| vtp | .00006324555320336759 | .0001 |
| Proton tauS | 1 | .4 |
| Electron tauS | 2 | .2 |
| Final scanned kpar | .06196773353931867 | .08 |

The length ratio is sqrt(.4*1.5)=sqrt(.6), not sqrt(.4).
Run both from the repository root:

```bash
mkdir -p data/example
./plume.e inputs/example/normalization_rho_ref.in
./plume.e inputs/example/normalization_d_ref.in
```

Both first mode rows give, at printed precision,
\(w=0.06024424-0.0004186511i\). All three written roots match between the two
representations. Their native wavevector and vtp columns differ by the factors
above. Convert those columns before overlaying curves. These example files
omit the numerical selectors and therefore use current defaults, not a
historical-release numerical convention.

<a id="scan-paths"></a>

### Scan paths

A beta-only d_ref scan fixes native kappa, a and other theta_s; changing beta_R
also updates theta_1. A reference theta_1 scan updates beta_R in return. Varying
both independently in a 2D scan is rejected. A beta-only rho_ref scan instead
holds legacy K, vtp and inverse temperature ratios fixed. Thus a beta-only scan
in the two conventions traces different physical systems. Convert every point,
including endpoint and interpolation rules, before making a physical comparison.

<a id="cold-references-and-the-auxiliary-thermal-scale"></a>

## Cold references and the auxiliary thermal scale

Cold-containing calculations require d_ref. For a cold reference use
`betap=0` and species_1 `tauS=0`. All species may be cold, or some may be kinetic.
No physical reference thermal speed or gyroradius exists in that limit; d_R,
v_A,R and Omega_R remain finite.

For warm species with a cold reference, the existing kinetic algebra uses an
internal positive scale b_star selected from the **first warm species'** theta.
It constructs auxiliary beta B=b_star, thermal velocity over c, V=a sqrt(b_star),
wavevector K_star=kappa sqrt(b_star), inverse temperature t_s=b_star/theta_s,
and normalization anisotropy alpha_star=1. The physical alpha_s is unchanged.
These factors cancel from physical resonance arguments, Bessel arguments,
frequency prefactors and output amplitudes. This is a coordinate transformation,
not an added species, fictitious physical temperature, new pressure or new
reference frame. Fully cold calculations bypass thermal conversion entirely.

There is no public b_star/auxiliary-temperature namelist control. Reordering warm
non-reference populations can change the internally chosen scale; results should
remain physically equivalent within the documented numerical conditioning.

<a id="diagnostic-units"></a>

## Diagnostic units

Let E=deltaE/deltaEx, B=deltaB/deltaEx and e0=deltaEx/B0. The Ex=1 gauge fixes
the shared phase/amplitude convention. Then

\[
 w\mathbf B=\frac{\boldsymbol\kappa\times\mathbf E}{a},\qquad
 \mathbf v_s=\frac{\delta\mathbf U_s}{c e_0},\qquad
 \eta_s=\frac{\delta n_s/n_{0s}}{e_0}.
\]

E and B are Gaussian field ratios; species velocities are **not** v_A-normalized
or thermal-speed-normalized. Density uses each species' own equilibrium n0s.
The physical moments follow by multiplying v_s by c e0 and eta_s by n0s e0.
All amplitudes can be complex. For drifting populations the current is
q_s(n0s deltaU_s+delta n_s U_s zhat), not just q_s n0s deltaU_s.

Corrected heating uses dimensionless Ps_s=Im(E-dagger chi_s E)/H. Formal energy
is W=abs(deltaEx)^2 H/(16 pi); H has separate vacuum electric, magnetic and
species terms. Their full definitions and applicability live on
[Heating and wave energy](heating.md), rather than being repeated in every
input/output table. FPC grids/header variables retain historical thermal units
and are a separate contract; see the [FPC appendix](README-JETPLUME.md).

The earlier normalization implementation audit is preserved as
[historical evidence](https://github.com/kgklein/PLUME/blob/main/tests/NORMALIZATION_IMPLEMENTATION_HISTORY.md).
