[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15215514.svg)](https://doi.org/10.5281/zenodo.15215514)

# PLUME: Plasma in a Linear Uniform Magnetized Environment

This is the PLUME code: Plasma in a Linear Uniform Magnetized Environment.

<img src="./PLUME_logo.png" alt="drawing" width="200"/>
<img src="./qrcode_plume_github.png" alt="drawing" width="200"/>

## Authors

Kristopher Klein   (kgklein@arizona.edu)
Gregory Howes      (gregory-howes@uiowa.edu)
Collin Brown       (collin.crbrown@gmail.com)

## Contents

1. What is PLUME?
2. Acknowledgements
3. Installing the PLUME Code
4. Running the PLUME Code
5. JET-PLUME
6. License

## 1. What is PLUME?

Plasma in a Linear Uniform Magnetized Environment (PLUME) is a numerical code that solves the Vlasov-Maxwell dispersion relation in a hot magnetised plasma.
PLUME allows for any number of particle species or components, assuming each can be described by a bi-Maxwellian distribution with a defined density, velocity parallel to the mean magnetic fiedl, and parallel and perpendicular temperatures.
The solver is able to identify supported waves with any direction of propagation with respect to the background magnetic field.
This code is based upon a F90 adaptation by Greg Howes of a solver originally by Eliot Quataert.

The calculation follows Stix 'Waves in Plasmas' Chapter 10, Eqns. 66-73.
The dispersion relation for $\omega/\Omega_{ref}$ is dependent on four global dimensionless parameters:

- Reference plasma beta: $\beta_{\parallel,\mathrm{ref}} = \frac{8 \pi n_{\mathrm{ref}} T_{\parallel,\mathrm{ref}}}{B^2}$
- Perpendicular wavevector: $k_\perp \rho_{ref}$
- Parallel wavevector: $k_\parallel \rho_{ref}$
- Parallel reference thermal velocity: $\sqrt{2 T_{\parallel,ref}/m_{ref}}/c$

     and six dimensionless parameters for component $s$:

- Parallel Temperature Ratio: $T_{\parallel,ref}/T_{\parallel,s}$
- Mass Ratio: $m_{ref}/m_{s}$
- Temperature Anisotropy: $T_{\perp}/T_{\parallel}|s$
- Charge Ratio: $q_{ref}/q_{s}$
- Density Ratio: $n_{s}/n_{ref}$
- Relative Velocity: $v_{s,drift}/v_{A,ref}$

The code then varies defined parameters to construct dispersion relations as a function of wavevector $(k_\perp \rho_{ref},k_\parallel \rho_{ref})$ or plasma parameters for the identified solutions.

Supplementary calculation of the associated heating rates or eigenfunctions can also be calculated.

## 2. Acknowledgements

If you use the code for a science publication,
1. please provide the code website [github.com/kgklein/PLUME](https://github.com/kgklein/PLUME) in the acknowledgements,

2. cite the DOI of the code:
```
@software{PLUME_2025,
  author       = {{Klein}, K. G. and
                  {Howes}, G. G.},
  title        = {kgklein/PLUME: Zenodo release},
  month        = April,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v1.0.1},
  doi          = {10.5281/zenodo.15215514},
  url          = {https://doi.org/10.5281/zenodo.15215514}
}
```

3. and cite the code paper:

[Klein, K. G., Howes, G. G.,
and Brown, C. R.: PLUME: Plasma in a Linear Uniform Magnetized Environment, RNAAS, 2025](https://iopscience.iop.org/article/10.3847/2515-5172/add1c2)
   
##  3. Installing the PLUME code

For advice on the installation of the code, please check [`INSTALL.md`](./INSTALL.md).

##  4. Running the PLUME code

There are two ways to run PLUME

###  Command line:

The first way is to use the command line. PLUME works with input files that specify the plasma and numerical parameters for the calculation.
Details on the input parameters can be found in [`input.md`](./input.md).

The values for the plasma parameters are extracted from *.in file, appended after the executable program call, e.g.
```
./plume.e inputs/example/run_example.in
```

A tutorial on running PLUME can be found in [`tutorial.md`](./tutorial.md).
Details on the outputs of the code can be found in [`output.md`](./output.md).

## Reference-length normalization

As of Nov 2, 2026, PLUME accepts `normalization='rho_ref'` or `normalization='d_ref'` in `&params`.
*Existing input files without `normalization` specified default to `'rho_ref'`.*
The reference species is `species_1` regardless of normalization choice.
The two matched examples are
[`normalization_rho_ref.in`](inputs/example/normalization_rho_ref.in) and
[`normalization_d_ref.in`](inputs/example/normalization_d_ref.in).

| Input | `rho_ref` | `d_ref` |
| --- | --- | --- |
| `kperp`, `kpar` | $k\rho_R$, $\rho_R=v_{t\perp,R}/|\Omega_R|$ | $kd_R$, $d_R=c/\omega_{pR}$ |
| `vtp` | $v_{t\parallel,R}/c$ | $v_{A,R}/c$ |
| `betap` | $\beta_{\parallel,R}$ | $\beta_{\parallel,R}=\tau_R^{(d)}$ |
| `tauS` | $T_{\parallel,R}/T_{\parallel,s}$ | $2T_{\parallel,s}/(m_Rv_{A,R}^2)$ |
| `muS` | $m_R/m_s$ | same |
| `alphS` | $T_{\perp,s}/T_{\parallel,s}$ | same |
| `Qs` | $q_R/q_s$ (signed) | same |
| `Ds` | $n_s/n_R$ | same |
| `vvS` | $U_s/v_{A,R}$ | same |
| `Kn` | $(n_0\sigma\rho_R)^{-1}$ | $(n_0\sigma d_R)^{-1}$ |

For finite reference temperature,
$k d_R=k\rho_R/\sqrt{\beta_{\parallel,R}\alpha_R}$,
$v_{t\parallel,R}/c=\sqrt{\beta_{\parallel,R}}\,v_{A,R}/c$, and
$\tau_s^{(d)}=\beta_{\parallel,R}/\tau_s^{(\rho)}$.
The frequency remains $\omega/\Omega_R$.
`d_ref` permits temperature to be expressed without a reference thermal length.
Collisionless `d_ref` dispersion calculations support both warm-reference
kinetic/cold mixtures (Stage 2), fully cold plasmas (Stage 3A), and
cold-reference mixtures with warm species (Stage 3B.2). Exactly
`tauS=0` selects a cold species. For a fully cold plasma set **every** species'
`tauS=0` and `betap=0`; `vtp` remains the positive Alfvén speed divided by `c`.
Cold `alphS` placeholders are ignored. All-cold runs support exact parallel,
exact perpendicular, and oblique propagation, with signed parallel drifts.
Cold-reference mixtures require `betap=0`. Warm responses use an auxiliary
thermal scale selected from the first warm species, preserving the physical
reference and every species temperature. Exact parallel and oblique propagation
are supported; exact perpendicular propagation with warm species is deferred.

Use `option=1`, `use_map=.false.`, `nscan=0` for supplied root guesses, or
`option=-1` for determinant evaluation. Frequency maps are also supported.
Cold-containing electromagnetic fields (Stage 4B.1) and species density/velocity
moments (Stage 4B.2) are available. Corrected collisionless kinetic absorption
and complete wave energy are available (Stage 4B.3A); cold kinetic heating is
unavailable (NaN). FPC, cold collisions and plasma-parameter/radial scans remain deferred. Stage 4A enables collisionless **wavevector** scans,
including cold-reference mixtures, using the existing continuation framework. Charge and current balance are warning-only conditions.
See the [input restrictions](input.md), [cold-reference report](tests/COLD_REFERENCE_MIXED_REPORT.md),
and examples:

- [Cold protons/warm electrons](inputs/example/mixed_cold_protons_warm_electrons.in).
- [Warm protons/cold electrons](inputs/example/mixed_warm_protons_cold_electrons.in).
- [Fully cold parallel L mode](inputs/example/cold_proton_electron.in).
- [Fully cold perpendicular X mode](inputs/example/cold_proton_electron_perpendicular.in).

Roots, modes, and susceptibility-tensor output columns retain the input convention.
FPC headers and grids retain the derived thermal solver convention. `d_ref` runs write
`data/<dataName>/normalization_<outputName>.txt`
to identify the convention;
legacy numeric output files are unchanged.
See [input.md](input.md) and [output.md](output.md) for scan and column details.
The [source audit and implementation report](normalization.md) records the
architecture choice, validation results, and remaining cold-fluid work.

###  Python wrapper:

The second way is to use the python wrapper, which simply helps one write inputs, call the code using the command line, and load the data into python variables. One way to use the python wrapper is with jupyter notebook (to install, see [here](https://jupyter.org/install)).

For `d_ref`, pass `normalization='d_ref'` to `plume_input.set_params(...)`,
provide `vtp=v_A,ref/c`, use $kd_R$ for both wavevectors, and set the first
species `tauS` equal to `betap`. `make_species(...)` continues to take the
same argument names; `tauS` then means $2T_{\parallel,s}/(m_Rv_{A,R}^2)$.
The wrapper reads and writes the normalization string in input files.

To do so, call jupyter notebook,
```
jupyter notebook
```
and navigate to and open the [`examplelinfpc.ipynb`](./examplelinfpc.ipynb) notebook. This notebook shows how to use the main functions of both PLUME and JET-PLUME. 


### Mixed cold/kinetic wavevector scans (Stage 4A)

Use the existing `&scan_input_N` inputs with `option=1` for serial continuation,
`option=2` for two-dimensional grids, or `option=4` for independent frequency
maps. Cold-containing runs require collisionless `d_ref`. For heating use the default
`heating_convention='corrected'`.
`eigen=.true.` enables electromagnetic fields and species density/velocity
moments; unavailable or failed species moments are marked NaN. All-kinetic scans retain their existing normalization and
diagnostic support. Both signs of **nonzero** `kpar` are supported. Stage 4A
rejects perpendicular scan points; fully cold single-wavevector perpendicular
calculations remain supported.

Supported coordinates are `kpar`, `kperp`, vector endpoints, angle at fixed
magnitude, and magnitude at fixed angle. Grids use distinct `kpar/kperp` axes or
ordered `(theta,magnitude)` axes. Species temperatures, beta, Alfvén speed and
drifts stay fixed. Cold-reference mixtures keep the same first-warm-species
auxiliary scale. Failed continuation points produce NaN rows and retain the
last valid seed; inspect warnings and reduce the wavevector step near crossings
or poles. There is no new branch tracker.

Working examples (run from the repository root):

```bash
./plume.e inputs/example/mixed_scan_warm_reference.in
./plume.e inputs/example/mixed_scan_cold_reference.in
./plume.e inputs/example/fully_cold_scan.in
```

See the [scan input details](input.md#mixed-wavevector-scans-stage-4a),
[output conventions](output.md#mixed-wavevector-scan-output-stage-4a), and
[implementation/validation report](tests/MIXED_SCAN_REPORT.md).

### Electromagnetic eigenfunctions (Stage 4B.1)

The existing cofactor method now supplies E/B polarization for cold and mixed
plasmas, including cold references. It retains Cartesian x/y/z ordering and the
historical **Ex=1** amplitude/phase convention. Magnetic fields follow Faraday's
law with native physical wavevector normalization; auxiliary temperatures do
not set field amplitudes. Existing all-kinetic eigenfunction/moment/heating
calculations remain unchanged.

Use `eigen=.true., heating=.false.` in existing serial or 2D scan inputs.
Fields retain the existing columns. Stage 4B.2 fills density/velocity slots with
cold-fluid moments for cold species and existing kinetic current/continuity
moments for warm species. Corrected mixed heating and wave energy are described below; FPC remains deferred.
At a known single-wavevector eigenfrequency, `option=-1` prints E/B and species moments; the field-only
Fortran API is also available. Arbitrary off-mode frequencies do not have a
validated eigenfunction. Modes with Ex=0 cannot use the historical gauge and
report unavailable rather than adopting a different normalization.

```bash
./plume.e inputs/example/mixed_fields_warm_reference.in
./plume.e inputs/example/mixed_fields_cold_reference.in
./plume.e inputs/example/cold_fields_parallel.in
./plume.e inputs/example/cold_fields_single.in
```

See [conventions and API](input.md#electromagnetic-eigenfunctions-stage-4b1),
[field columns](output.md#electromagnetic-eigenfunction-conventions-stage-4b1),
and the [audit and independent validation report](tests/MIXED_EIGENFUNCTION_REPORT.md).

### Species density and velocity eigenfunctions (Stage 4B.2)

All supported cold/kinetic compositions now provide species velocity and density
perturbations, preserving species ordering and the existing Ex=1 field gauge.
Velocity is `deltaU/(c*deltaEx/B0)` and density is
`(delta n/n0)/(deltaEx/B0)`. Cold species solve the linearized momentum and
continuity equations, including the equilibrium-drift magnetic force. Warm
species retain the susceptibility-current/continuity method. Physical moments
never depend on the auxiliary thermal scale or cold anisotropy placeholders.

Existing `eigen=.true., heating=.false.` serial/2D scan files keep their columns.
A failure in one species produces NaN in that species' moment block, with a
species/status warning; valid E/B and other species remain available. An absent
zero-density species has no defined relative-density moment. Single-point
`option=-1` prints moments at a valid eigenfrequency. Stage 4B.3A adds corrected heating/energy below; FPC remains deferred; the all-kinetic historical path remains unchanged.

The four field examples above now also provide moments. A multiple cold/warm
population example is:

```bash
./plume.e inputs/example/mixed_moments_multispecies.in
```

See the [moment inputs/API](input.md#species-moments-stage-4b2),
[output conventions](output.md#species-moment-conventions-stage-4b2), and
[independent validation report](tests/COLD_FLUID_EIGENFUNCTION_REPORT.md).

### Plasma dispersion function backend

PLUME evaluates the Landau plasma dispersion function
$Z(\zeta)=i\sqrt{\pi}\,w(\zeta)$ using Steven G. Johnson's bundled
[Faddeeva implementation](third_party/faddeeva/README.md). The default
`z_backend='faddeeva'` fixes the historical evaluator's accuracy loss near
$|\zeta|=3$ and its unassigned-result failures. It requires only an ordinary
C compiler in addition to Fortran; no library installation or Python runtime
is needed. The upstream MIT license and attribution are bundled with the source.

These settings are independent:

```fortran
! Recommended defaults
z_backend = 'faddeeva'
bessel_mode = 'adaptive'

! Historical compatibility (select these instead)
! z_backend = 'legacy'
! bessel_mode = 'fixed128'
```

Fixed harmonics alone do not restore the old Z evaluator. The retained legacy
backend reproduces its historical approximations and failure behavior. Accurate
Z does not remove Bessel approximation errors or cancellation in kinetic sums.
Invalid or unrepresentable accurate evaluations return an invalid dispersion
status and cannot be accepted as roots. See the [input reference](input.md#plasma-dispersion-function-backend),
[working examples](inputs/example/faddeeva_accuracy.in), and
[accuracy, performance and compatibility report](tests/FADDEEVA_BACKEND_REPORT.md).

### Adaptive kinetic harmonics

Kinetic calculations now select a separate symmetric harmonic range for each
warm species. The defaults in `&params` are `bessel_mode='adaptive'`,
`bessel_tol=1.e-10`, and `bessel_nmax=512`. Cold species need no Bessel arrays.
Repeated frequency evaluations reuse coefficients at the same physical
perpendicular Bessel argument, including with auxiliary normalization.

Use `z_backend='legacy'` together with `bessel_mode='fixed128'` to
reproduce historical collisionless numerical results exactly with the tested
GNU toolchain. `bessel_mode='fixed', bessel_nfixed=128` is an equivalent spelling. Other fixed counts are available for convergence studies.
A Bessel-tail tolerance controls truncation, **not total susceptibility or root
error**. Tighten the tolerance and increase the maximum to check a result.
Reaching the adaptive maximum emits a warning and continues; unresolved points
are recorded in the `.bessel` sidecar. Set `bessel_diagnostics=.true.` for cache
statistics and species cutoffs.

The [large-argument example](inputs/example/adaptive_bessel.in) demonstrates a
mode requiring more than 128 harmonics. See the [input controls](input.md#bessel-harmonic-controls),
[metadata description](output.md#bessel-harmonic-metadata), and
[implementation and validation report](tests/ADAPTIVE_BESSEL_CACHE_REPORT.md).

## 5. JET-PLUME

<img src="./Jet-Plume_Logo.svg" alt="JetPlumeLogdrawing" width="200"/>

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique and linear theory. Please see the [README for JET-PLUME](./README-JETPLUME.md). JET-PLUME is automatically installed with PLUME. It's inclusion has no impact on the use of PLUME.

If you use JET-PLUME (i.e. not just PLUME) for a science publication, please follow the directions in the Acknowledgements section above AND in the Acknowledgements of your publication and cite the paper for JET-PLUME: [... Work in Progress]

With the creation of JET-PLUME, a wrapper to use PLUME and JET-PLUME in a jupyter notebook was created. It aids in the creation of input files and running of the code with said input files. Please see the [example notebook](./page/examplelinfpc.md) to see the key features of this wrapper. The use of the wrapper is entirely optional. NOTE: Please create and run notebooks from the main directory of the repository. This is necessary because the notebooks import key functions directly from the linfpclib module located there.

JET-PLUME is operated in a similar manner to PLUME, either using command line or the python wrapper. See Section 4 above on how to operate and [README for JET-PLUME](./README-JETPLUME.md) for the additional inputs.

## 6. License

BSD 2-Clause License

Copyright (c) 2025, Kristopher G. Klein and Gregory G. Howes and Collin R. Brown
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### Corrected heating and wave energy (Stage 4B.3A)

`heating_convention='corrected'` is now the default, independent of Z, Bessel,
and input normalization settings. Collisionless single-point (`option=-1`),
serial and 2D wavevector scans (`heating=.true.`) support kinetic, fully cold,
and mixed populations, including cold references. The full-plasma E/B fields
and species moments are preserved if heating is unavailable.

The species diagnostic is signed normalized absorption/emission,
`Ps = Im(E† chi_s(omega_r) E)/H`, approximately a contribution to
`-gamma/omega_r` for an isolated weak mode. It is neither dimensional power nor
irreversible thermodynamic heating. All seven selected historical channels
are retained; they are not an exhaustive partition. Cold species have NaN
kinetic heating, while contributing to wave energy through their susceptibility.

The complete formal wave energy includes vacuum electric, magnetic, and every
species' `d(omega_r*chi_H)/domega_r` contribution. A versioned `.energy`
companion reports `H_total`, `H_E`, `H_B`, individual species contributions,
signed fractions and the noncancelling electric contribution `f_E_abs`.
Signed fractions can be negative or exceed one. The latter is a contribution
measure, not a physical energy fraction. Check numerical and applicability
statuses before interpreting strongly growing/damped or reactive modes.
See [input controls](input.md#heating-and-wave-energy-stage-4b3a),
[output schema](output.md#corrected-energy-companion-plume-energy-v1), and the
[implementation report](tests/MIXED_HEATING_ENERGY_REPORT.md).

Exact historical heating reproduction requires all three independent settings:

```fortran
heating_convention='legacy'
z_backend='legacy'
bessel_mode='fixed128'
```

The legacy thermal heating implementation remains unchanged, including its
incomplete wave-energy normalizer and historical harmonic snapshot behavior.
Corrected heating requires `collision_type=0`; use explicit legacy mode for
previously supported all-kinetic collisional diagnostics. FPC retains its
historical implementation and availability restrictions. Corrected heating is
not added to maps, radial scans, or plasma-parameter scans.

Working examples (run `./plume.e <input>`):

- [All kinetic](inputs/example/heating_corrected_kinetic.in).
- [Warm-reference mixture](inputs/example/heating_corrected_warm_reference.in).
- [Cold-reference mixture](inputs/example/heating_corrected_cold_reference.in).
- [Fully cold wave energy](inputs/example/heating_corrected_cold.in).
- [Historical convention](inputs/example/heating_legacy.in).
