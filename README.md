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
Cold-containing eigenfunctions, heating, FPC, collisions, and parameter scans
are unavailable. Charge and current balance are warning-only conditions.
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
