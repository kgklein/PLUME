<a id="plasma-species-models-and-capabilities"></a>

# Plasma species, models and capabilities

[Manual home](./) · [Normalization](normalization.md) · [Modes and scans](modes.md)

<a id="kinetic-populations"></a>

## Kinetic populations

A warm population s is a drifting gyrotropic bi-Maxwellian. With number density
n0s, parallel drift U_s and thermal speeds
\(v_{t\parallel,s}=\sqrt{2T_{\parallel,s}/m_s}\),
\(v_{t\perp,s}=\sqrt{2T_{\perp,s}/m_s}\), its equilibrium distribution is

\[
 f_{0s}(v_\parallel,v_\perp)=
 \frac{n_{0s}}{\pi^{3/2}v_{t\perp,s}^{2}v_{t\parallel,s}}
 \exp\!\left[-\frac{v_\perp^2}{v_{t\perp,s}^2}
              -\frac{(v_\parallel-U_s)^2}{v_{t\parallel,s}^2}\right].
\]

The magnetic field and drift are along z. Kinetic anisotropy alpha_s is
T_perp,s/T_parallel,s>0. In d_ref theta_s=`tauS`>0, however small, selects this
model. In rho_ref `tauS` is the inverse reference/species temperature ratio and
must be positive. No threshold switches a sufficiently cold Maxwellian to fluid.
Multiple components of the same particle species may have separate temperatures,
densities and drifts; the input order is retained in every species-resolved output.

The kinetic susceptibility contains finite-Larmor-radius and resonant harmonic
response. It is not an arbitrary-distribution, relativistic, or general
perpendicular-drift solver. The limited all-kinetic finite-rate option is a
nonconservative Krook response, not a general Coulomb collision model.
See [numerical methods](numerics.md) and [limitations](limitations.md).

<a id="cold-fluid-populations"></a>

## Cold-fluid populations

Under d_ref, **exactly** `tauS=0` selects cold fluid. Negative tauS is invalid.
A cold population has zero thermal pressure; `alphS` is physically undefined and
ignored by its response. Supply a conventional placeholder (for example 1) for
clear provenance. A positive but tiny tauS remains kinetic.

Cold species retain finite inertia, signed cyclotron frequency, density,
charge and parallel drift. Their response includes oblique drift coupling, not
just the stationary Stix transverse terms. They support complex frequencies,
fully cold systems and a cold reference population. The kernel does not divide
by a temperature or call Z/Bessel functions. Exact Doppler/cyclotron poles are
reported, not regularized by artificial damping.

Cold populations have fluid density/velocity eigenfunctions when the Ex=1 field
solution and response are valid. **Kinetic** heating, kinetic harmonic heating
and FPC are unavailable for cold species. NaN means unavailable, not a measured
physical zero. Their reversible response still contributes to formal wave energy.

<a id="mixtures-and-species-ordering"></a>

## Mixtures and species ordering

PLUME stores ordinary susceptibility chi_s separately for each population.
The dielectric response is epsilon=I+sum_s chi_s; it adds the vacuum identity
once, then forms the electromagnetic wave matrix. The full plasma determines
the common field polarization; moments and kinetic absorption use that field
with the appropriate per-species response. See [Eigenfunctions](eigenfunctions.md)
and [Heating](heating.md).

Species_1 remains the physical reference even if cold. Cold-reference/warm
mixtures use the internal auxiliary scale described in [Normalization](normalization.md#cold-references-and-the-auxiliary-thermal-scale).
Do not reorder species merely to obtain a warm reference. Interleaved warm/cold
species and several warm populations are supported, subject to geometry and
workflow restrictions below.

<a id="capability-matrix"></a>

## Capability matrix

S = supported in the indicated scope; P = partially supported with conditions;
U = unsupported; NV = not generally validated. All cold columns require
`normalization='d_ref'` and `collision_type=0`. The kinetic column allows either
normalization. Support for dispersion is not support for every diagnostic.

| Workflow/geometry | Kinetic only | Fully cold | Mixed, warm reference | Mixed, cold reference |
|---|---|---|---|---|
| Oblique direct evaluation / supplied roots | S | S | S | S |
| Linear frequency-map root search | S | S | S | S |
| Serial / 2D wavevector scans | S | S | S | S |
| Repeated frequency maps | P: no supplementary diagnostics | P: same | P: same | P: same |
| Exact parallel kperp=0 | U: warning, historical singular path | S | S | S |
| Exact perpendicular kpar=0 single wavevector | U | P: dispersion supported; Ex=1 gauge may fail | U | U |
| Any wavevector scan containing kpar=0 | U | U | U | U |
| Historical plasma-parameter / radial scans | P: legacy restrictions | U | U | U |
| Finite-rate Krook response | P / NV off-diagonal accuracy | U | U | U |
| General collision operator | U | U | U | U |

For scans, even a hidden substep at kpar=0 is prohibited. Fully cold
perpendicular single-point O/X dispersion does not imply perpendicular scans.
Parallel longitudinal or perpendicular ordinary modes may have Ex=0 and thus
no output eigenfunction in PLUME's chosen gauge. All-kinetic kperp=0 remains
unsupported even though some mixed/cold parallel tests pass.

| Diagnostic | Kinetic only | Fully cold | Either mixed reference |
|---|---|---|---|
| E/B polarization | S for valid supported Ex=1 modes | S subject to gauge/poles | S subject to gauge/poles |
| Species U and density | S | S fluid moments | S kinetic/fluid per population |
| Total species kinetic absorption + seven channels | Corrected collisionless; legacy optional | U / NaN | Corrected warm species only; cold NaN |
| Complete formal wave energy | S with real-frequency weak-mode qualifications | S with same qualifications | S including all physical species |
| Corrected finite-rate heating | U | U | U |
| Gyro/Cartesian FPC | P: historical all-kinetic pathway | U | U, including warm members |
| General FPC energy/moment closure | NV | U | U |

<a id="equilibrium-checks-are-warnings"></a>

## Equilibrium checks are warnings

With the normalization symbols defined previously, charge and current imbalance
are proportional to \(\sum_s D_s/Q_s\) and \(\sum_s D_su_s/Q_s\), respectively.
Input validation warns at absolute magnitude at least 1e-6. A drifting reference
also produces a warning. None of these checks changes density, drift, species
order or reference frame, and neither neutrality nor zero current is mandatory.
Some nonstandard homogeneous equilibria require external charge, currents or
forces to maintain the assumed background. Users must make that equilibrium
interpretation explicit.

<a id="existing-examples"></a>

## Existing examples

Run from the repository root after `mkdir -p data/example`:

| Composition | Command | Inspect |
|---|---|---|
| Kinetic proton/electron | `./plume.e inputs/example/heating_corrected_kinetic.in` | Mode E/B/moments; kinetic Ps; .energy |
| Warm proton / cold electron | `./plume.e inputs/example/mixed_fields_warm_reference.in` | Three finite field/moment rows |
| Cold proton / warm electron | `./plume.e inputs/example/mixed_fields_cold_reference.in` | Same units, beta_R=0, no fictitious reference temperature |
| Fully cold | `./plume.e inputs/example/heating_corrected_cold.in` | Fluid moments, finite formal energy, NaN kinetic heat |
| Multiple interleaved populations | `./plume.e inputs/example/mixed_moments_multispecies.in` | All species in original order, signed drift currents |

[Examples and expected results](examples.md) links the actual files and identifies
which additional curriculum inputs have not yet been created.
