<a id="normalization-and-physical-units"></a>

# Normalization and physical units

[Manual home](index.md) · [Input reference](input.md) · [Species models](models.md)

Here, we define the temporal and spatial scales used in PLUME.
For brevity, temperatures are in energy units (Boltzmann's constant is absorbed into T).
Electromagnetic formulas are in Gaussian units. 

The equilibrium is homogeneous, with
mean magnetic field
\(\mathbf B_0=B_0\hat{\mathbf z}\) and
wavevector
\(\mathbf k=(k_\perp,0,k_\parallel)\). 
Time dependence is given by \(\exp(-i\omega t)\).
The solution \(\omega=\omega_r+i\gamma\) grows when \(\gamma>0\).

<a id="reference-scales"></a>

## Reference scales

PLUME always treats `species_1` as reference R.
This species has a positive reference gyrofrequency scale (a positive reference charge and B0), 
with the other speices retaining signed charge ratios for other populations. 

For completeness, we define

\[
 \Omega_R=\frac{q_R B_0}{m_R c},\qquad
 \omega_{pR}^2=\frac{4\pi n_R q_R^2}{m_R},\qquad
 v_{A,R}=\frac{B_0}{\sqrt{4\pi n_R m_R}},\qquad
 d_R=\frac{c}{\omega_{pR}}=\frac{v_{A,R}}{\Omega_R}.
\]

Here c is the speed of light, and
\(m_R\), \(q_R\), and \(n_R\) set 
the reference mass, charge, and number density. 
For the thermal reference define

\[
 v_{t\parallel,R}=\sqrt{\frac{2T_{\parallel,R}}{m_R}},\quad
 \alpha_R=\frac{T_{\perp,R}}{T_{\parallel,R}},\quad
 \rho_R=\frac{v_{t\perp,R}}{|\Omega_R|},\quad
 \beta_R=\frac{8\pi n_R T_{\parallel,R}}{B_0^2}.
\]

The gyroradius uses the **perpendicular** thermal speed, while plasma beta uses the **parallel** reference temperature.
This yields \(v_{t\parallel,R}/v_{A,R}=\sqrt{\beta_R}\) and
\(\rho_R/d_R=\sqrt{\beta_R\alpha_R}\).
PLUME's Alfvén speed uses reference density/mass, **not** the mixture's total mass density.
*Do not substitute a total-mixture Alfvén speed without converting all normalized inputs.*

<a id="inputs-in-the-two-conventions"></a>

## Inputs in the two conventions

The normalizations are selected using `normalization='rho_ref'` (default) or `'d_ref'` in `&params`.
Let
\(w=\omega/\Omega_R\), \(a=v_{A,R}/c\),
\(K_{\rho,j}=k_j\rho_R\), and \(K_{d,j}=k_jd_R\), with j parallel or perpendicular.
For a physical population s, use \(\mu_s=m_R/m_s\),
\(Q_s=q_R/q_s\) (signed), \(D_s=n_s/n_R\),
\(u_s=U_s/v_{A,R}\), and \(\alpha_s=T_{\perp,s}/T_{\parallel,s}\).

| Input | `rho_ref` | `d_ref` |
|---|---|---|
| `betap` | Reference parallel beta, beta_R | Same; zero for a cold reference |
| `kperp`, `kpar` | \(k_\perp \rho_R\), \(k_\parallel \rho_R\) | \(k_\perp d_R\), \(k_\parallel d_R\) |
| `vtp` | \(v_{t,\parallel,R}/c\) | \(a=v_{A,R}/c\) |
| `tauS` | \(\tau_{\rho,s}=T_{\parallel,R}/T_{\parallel,s}\) | \(\theta_s=2T_{\parallel,s}/(m_R v_{A,R}^2)\) |
| `muS` | \(m_R/m_s\) | Same |
| `Qs` | \(q_R/q_s\), signed | Same |
| `Ds` | \(n_s/n_R\) | Same |
| `alphS` | \(T_{\perp,s}/T_{\parallel,s}\) | Same for warm species; ignored for cold |
| `vvS` | \(u_s=U_s/v_{A,R}\) | Same; never thermal-speed-normalized |
| `g_om`, `g_gam` | \(\omega_r/\Omega_R\), \(\gamma/\Omega_R\) | Same |
| `Kn` | Neutral mean free path/\(\rho_R\) | Neutral mean free path/\(d_R\) |

We use \(\theta_s\) as a symbol for the temperature definition `d_ref` `tauS` in this manual
to avoid confusing a physical temperature with the inverse ratio in `rho_ref`.
*Both use the same Fortran input name.*
Self ratios for the reference population are `mu_1`=`Q_1`=`D_1`=1.
The solver does not enforce self-ratio or reference-density consistency;
unusual component decompositions require a clearly specified reference density.
**There is no automatic renormalization.**

For `d_ref` the physical thermal and species frequency relations are

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

When \(beta_R>0\) and \(alpha_R>0\),

\[
 k_{\rho,j}=k_{d,j}\sqrt{\beta_R\alpha_R},\qquad
 v_{t,\parallel,R}/c=a\sqrt{\beta_R},\qquad
 \tau_{\rho,s}=\frac{\beta_R}{\theta_s}.
\]

\(\mu_s\), \(Q_s\), \(D_s\), \(\alpha_s\), \(u_s\) and \(w\) are unchanged.
Matching collision frequencies additionally requires
\(\mathrm{Kn}_d=\mathrm{Kn}_\rho\sqrt{\beta_R\alpha_R}\).
The existing derived rates are

\[
 \frac{\nu_s}{\Omega_R}=\frac{1}{\sqrt{2}\,\mathrm{Kn}_\rho}
       \sqrt{\frac{\mu_s}{\tau_s^{(\rho)}\alpha_R}}
   =\frac{\sqrt{\mu_s\theta_s}}{\sqrt{2}\,\mathrm{Kn}_d}.
\]

These rates belong only to the restricted all-kinetic Krook implementation,
**which is a work in progress.**

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

The length ratio is sqrt(.4*1.5)=sqrt(.6).
Run both from the repository root:

```bash
mkdir -p data/example
./plume.e inputs/example/normalization_rho_ref.in
./plume.e inputs/example/normalization_d_ref.in
```

Both first mode rows give, at printed precision,
\(w=0.06024424-0.0004186511i\).

**Show that these two solutions match.**
Their native wavevector and vtp columns differ by the factors above.
Convert those columns and then plot overlaying curves. 

<a id="scan-paths"></a>

### Scan paths


*Caveat:* One of the useful aspects of PLUME is the ability to serially or in parallel
scan dispersion relation solutions over wavevector and plasma parameters.
Generally, the user can select any sets of parameters for these scans.
The exception is for a beta_R and theta_s `d_ref` scans.
By definition for the `d_ref` normalization, changing beta_R also updates theta_1.
A reference theta_1 scan updates beta_R in return.
Varying both independently in a 2D scan is rejected.
A beta-only `rho_ref` scan traces different a different physical system than that in `d_ref`.
Convert every point, including endpoint and interpolation rules, before making a physical comparison.

<a id="cold-references-and-the-auxiliary-thermal-scale"></a>

## Cold references and the auxiliary thermal scale

Cold-containing calculations require selecting the `d_ref` normalization.
For a cold reference use `betap=0` and species_1 `tauS=0`.
All species may be cold, or some may be kinetic.
No physical reference thermal speed or gyroradius exists in that limit,
but \(d_R\), \(v_{A,R}\) and \(\Omega_R\) remain finite.

For warm species with a cold reference,
the existing kinetic algebra uses an internal positive scale \(b_{\rm star}\) selected from the **first warm species'** theta.
It constructs auxiliary beta \(B=b_{\rm star}\),
thermal velocity over c, \(V=a \sqrt(b_{\rm star})\),
wavevector \(k_{\rm star}=K_{d,R} \sqrt(b_{\rm star}\)),
inverse temperature \(t_s=b_{\rm star}/\theta_s\),
and normalization anisotropy \(\alpha_{\rm star}=1\).
The physical \(\alpha_s\) is unchanged.
This is simply a coordinate transformation.
Fully cold calculations bypass this thermal conversion entirely.

There is no public \(b_{\rm star}\)/auxiliary-temperature namelist control.
Reordering warm non-reference populations can change the internally chosen scale,
but results should remain physically equivalent within the documented numerical conditioning.

<a id="diagnostic-units"></a>

## Diagnostic units

Let \( E=\delta E/\deltaE_x\), \(B=\delta B/\delta E_x\) and \(e_0=\delta E_x/B_0\).
We choose the \(E_x=1\) gauge to fix the shared phase and amplitude eigenvalue convention.
This lets us write:

\[
 w\mathbf B=\frac{\boldsymbolK_{d}\times\mathbf E}{a},\qquad
 \delta \mathbf v_s=\frac{\delta\mathbf U_s}{c e_0},\qquad
 \eta_s=\frac{\delta n_s/n_{0s}}{e_0}.
\]

**Double check velocity perturbation normalization in `d_ref` units.**

E and B are Gaussian field ratios.
Species velocity fluctuations are **not** v_A-normalized
or thermal-speed-normalized.
Density uses each species' own equilibrium \(n_0s\).
The physical moments follow by multiplying \(v_s\) by \(c e_0\) and \(\eta_s\) by \(n_0 e_0\).
All amplitudes are generally complex.
For drifting populations the current is
\(q_s(n_{0s} \delta U_s+\delta n_s U_s \zhat)\),
not just \(q_s n_{0s} \delta U_s\).

Heating uses dimensionless

\[
Ps_s=\Im\left(E^{\dagger} \chi_s E\right)/H.
\]

Formal energy is

\[
W=\abs(\delta E_x)^2 H/(16 \pi)
\]

where H has separate vacuum electric, magnetic and species terms.
Their full definitions and applicability live in [Heating and wave energy](heating.md).
FPC grids and header variables retain historical thermal units;
see the [FPC appendix](README-JETPLUME.md).
