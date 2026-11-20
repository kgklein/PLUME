<a id="plasma-species-models-and-capabilities"></a>

# Plasma species, models and capabilities

[Manual home](index.md) · [Normalization](normalization.md) · [Modes and scans](modes.md)

<a id="kinetic-populations"></a>

PLUME supports bi-Maxwellian and cold populations in mixtures,
subject to the restrictions below.
If you use one or more cold populations, you **must** use the `d_ref` normalization.

## Kinetic populations

We define a warm population \(s\) as a drifting gyrotropic bi-Maxwellian
with number density \(n_{0s}\), 
parallel drift \(U_s\) and 
thermal speeds 
\(v_{t\parallel,s}=\sqrt{2T_{\parallel,s}/m_s}\),
\(v_{t\perp,s}=\sqrt{2T_{\perp,s}/m_s}\). 
Its equilibrium distribution is

\[
 f_{0s}(v_\parallel,v_\perp)=
 \frac{n_{0s}}{\pi^{3/2}v_{t\perp,s}^{2}v_{t\parallel,s}}
 \exp\!\left[-\frac{v_\perp^2}{v_{t\perp,s}^2}
              -\frac{(v_\parallel-U_s)^2}{v_{t\parallel,s}^2}\right].
\]

The magnetic field and drift are oriented along z. 
Kinetic anisotropy \(\alpha_s = T_{\perp,s}/T_{\parallel,s}>0\). 
In `d_ref`, \(theta_s\)=`tauS`>0 sets the population to be warm, 
while `tauS=0` sets the population to be cold.  
In `rho_ref` `tauS` is the reference-to-species temperature ratio and must be positive. 
*No threshold switches a sufficiently narrow Maxwellian to cold fluid.*
Multiple components of the same particle species may be separated into individual componets, with each population having separate temperatures,
densities and drifts.

<a id="cold-fluid-populations"></a>

## Cold-fluid populations

Under `d_ref`, setting `tauS=0` selects a cold fluid model. 
*Negative tauS is invalid.*
A cold population has zero thermal pressure, 
and thus `alphS` is physically undefined and ignored by its response. 
Supply it with 1 as a  placeholder.

Cold species retain finite inertia, signed cyclotron frequency, density,
charge and parallel drift. 
Their response includes oblique drift coupling, not just the stationary Stix transverse terms. 
The response supports complex frequencies, fully cold systems and a cold reference population. 

Cold populations have fluid density/velocity eigenfunctions when the Ex=1 field
solution and response are valid. 
**Kinetic** heating, kinetic harmonic heating and FPC are unavailable for cold species, and report NaN.
Their reversible response does, however, contribute to formal wave energy.

<a id="mixtures-and-species-ordering"></a>

## Mixtures and species ordering

PLUME stores ordinary susceptibility \(chi_s\) separately for each population.
The dielectric response is \(\epsilon=I+\sum_s \chi_s\), 
adding the vacuum identity and forming the electromagnetic wave matrix. 
The full plasma determines the common field polarization, 
with moments and kinetic absorption use that field with the appropriate per-species response. 
See [Eigenfunctions](eigenfunctions.md)
and [Heating](heating.md) for the full expressions.

`Species_1` remains the physical reference even if cold. 
Cold-reference/warm mixtures use the internal auxiliary scale described in [Normalization](normalization.md#cold-references-and-the-auxiliary-thermal-scale),
so there is no need to reorder species merely to obtain a warm reference. Interleaved warm/cold species and several warm populations are supported, subject to the capability below.

<a id="capability-matrix"></a>

## Capability matrix

A list of current and future planned capabilities.

- S = supported in the current release; 
- P = partially supported with conditions;
- U = unsupported; 
- NV = not fully validated. 

All cold columns require
`normalization='d_ref'` and `collision_type=0`. 
The fully kinetic solutions allow for either normalization. 

| Workflow/geometry | Kinetic only | Fully cold | Mixed, warm reference | Mixed, cold reference |
|---|---|---|---|---|
| Oblique direct evaluation / supplied roots | S | S | S | S |
| Frequency-map root search | S | S | S | S |
| Serial / 2D wavevector scans | S | S | S | S |
| Repeated frequency maps | P | P | P | P |
| Exact parallel kperp=0 | U  | S | S | U |
| Exact perpendicular kpar=0 single wavevector | U | P: Ex=1 gauge may fail | U | U |
| Any wavevector scan containing kpar=0 | U | U | U | U |
| Plasma-parameter scans | S: corrected collisionless `rho_ref`; historical restrictions otherwise | U | U | U |
| Radial scans | P: legacy restrictions | U | U | U |
| Finite-rate Krook response | P / NV off-diagonal accuracy | U | U | U |
| General collision operator | U | U | U | U |

For scans, even an intermediate step at kpar=0 is prohibited. 
Fully cold perpendicular single-point O/X dispersion does not imply perpendicular scans.
Parallel longitudinal or perpendicular ordinary modes may have Ex=0 and thus
no output eigenfunction in PLUME's chosen gauge. 

| Diagnostic | Kinetic only | Fully cold | Either mixed reference |
|---|---|---|---|
| E/B polarization | S for valid supported Ex=1 modes | S subject to gauge/poles | S subject to gauge/poles |
| Species U and density | S | S fluid moments | S kinetic/fluid per population |
| Total species kinetic absorption + seven channels | Corrected collisionless; legacy optional |  NaN | Corrected warm species only; cold NaN |
| Complete formal wave energy | S with real-frequency weak-mode qualifications | S with same qualifications | S including all physical species |
| Corrected finite-rate heating | U | U | U |
| Gyro/Cartesian FPC | P: historical all-kinetic pathway | U | U, including warm members |
| General FPC energy/moment closure | NV | U | U |

<a id="equilibrium-checks-are-warnings"></a>

## Equilibrium checks are warnings

With the symbols defined in [Normalization](normalization.md), charge and current
imbalance are proportional to \(\sum_s D_s/Q_s\) and \(\sum_s D_su_s/Q_s\), respectively.
Input validation warns at absolute magnitude at least 1e-6. 
A drifting reference also produces a warning. 
None of these checks changes density, drift, species order or reference frame, and neither charge neutrality nor zero current is mandatory.
Some nonstandard homogeneous equilibria require external charge, currents or
forces to maintain the assumed background. 
Users must make that equilibrium interpretation explicit.

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

[Examples and expected results](examples.md) lists additional input files for users to explore.
