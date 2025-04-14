[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15215514.svg)](https://doi.org/10.5281/zenodo.15215514)

# PLUME: Plasma in a Linear Uniform Magnetized Environment

This is the PLUME code: Plasma in a Linear Uniform Magnetized Environment.

<img src="./PLUME_logo.png" alt="drawing" width="200"/>
<img src="./qrcode_plume_github.png" alt="drawing" width="200"/>

## Authors

Kristopher Klein   (kgklein@arizona.edu)
Gregory Howes      (gregory-howes@uiowa.edu)

## Contents

1. What is PLUME?
2. Acknowledgements
3. Installing the PLUME Code
4. Running the PLUME Code
5. License

## 1. What is PLUME?

Plasma in a Linear Uniform Magnetized Environment (PLUME) is a numerical code that solves the Vlasov-Maxwell dispersion relation in a hot magnetised plasma.
PLUME allows for any number of particle species or components, assuming each can be described by a bi-Maxwellian distribution with a defined density, velocity parallel to the mean magnetic fiedl, and parallel and perpendicular temperatures.
The solver is able to identify supported waves with any direction of propagation with respect to the background magnetic field.
This code is based upon a F90 adaptation by Greg Howes of a solver originally by Eliot Quataert.

The calculation follows Stix 'Waves in Plasmas' Chapter 10, Eqns. 66-73.
The dispersion relation for $\omega/\Omega_{ref}$ is dependent on four global dimensionless parameters:

- Reference plasma beta: \(\beta_{\parallel,ref} = \frac{8 \pi n_{ref} T_{\parallel,ref}}{B^2}\)
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
and Brown, C. R.: PLUME: Plasma in a Linear Uniform Magnetized Environment, RNAAS, 2025 in prep](http://doi.org/10.1017/S0022377818000739)
   
##  3. Installing the PLUME code

For advice on the installation of the code, please check [`INSTALL.md`](./INSTALL.md).

##  4. Running the PLUME code

PLUME works with input files that specify the plasma and numerical parameters for the calculation.
Details on the input parameters can be found in [`input.md`](./input.md).

The values for the plasma parameters are extracted from *.in file, appended after the executable program call, e.g.
```
./plume.e inputs/example/run_example.in
```

A tutorial on running PLUME can be found in [`tutorial.md`](./tutorial.md).
Details on the outputs of the code can be found in [`output.md`](./output.md).

## 5. License

BSD 2-Clause License

Copyright (c) 2025, Kristopher G. Klein and Gregory G. Howes
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
