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

PLUME is a numerical code that solves the Vlasov-Maxwell dispersion
relation in hot magnetised plasma.
PLUME allows for any number of particle species, assuming each can be described by a bi-Maxwellian distribution with a defined density, velocity, and parallel and perpendicular temperature.
The solver is able to identify supported waves with any direction of propagation with respect to the background magnetic field.

This code uses an F90 adaptation (by Greg Howes) of the Hot Plasma 
       Dispersion Relation originally by Eliot Quataert.

 PLUME calculates the hot plasma dispersion relation for a plasma with 
       an arbitrary number of ion and electron species with relative drifts
       and bi-Maxwellian velocity distributions.
       The calculation follows Stix Chapter 10 eqn 66-73.
 The Dispersion relation for omega/Omega_ref
     is dependent on four global dimensionless parameters:

       betap: Plasma Reference Beta:               8 pi n_ref T_ref /B^2
       kperp: Perpendicular wavelength:         kperp rho_ref
       kpar : Parallel wavelength:              kparallel rho_ref
       vtp  : Parallel proton thermal velocity: sqrt(2 T_||ref/m_ref)/c

     and six dimensionless component parameters:

       tau_s : T_ref/T_s
       mu_s  : m_ref/m_s
       alph_s: T_perp/T_parallel|s
       Q_s   : q_ref/q_s
       D_s   : n_s/n_ref
       vv_s  : v_s,drift/v_Aref

The values for these parameters are extracted from *.in file, appended after
    the executable program call.

The code then varies defined parameters to construct dispersion relations as a function of wavevector or plasma parameter for the identified modes.
Supplementary calculation of the associated heating rates or eigenfunctions can also be calculated.

## 2. Acknowledgements

If you use the code for a science publication,
1. please provide the code website

lorem ipsum

in the acknowledgements,

2. cite the DOI of the code:
lorem ipsum

3. and cite the code paper:
   
lorem ipsum

##  3. Installing the PLUME code

For advice on the installation of the code, please check [`INSTALL.md`](./INSTALL.md)

##  4. Running the PLUME code

PLUME works with input files that specify the plasma and numerical parameters for
the calculation.

## 5. License

BSD 2-Clause License

Copyright (c) 2025, Kristopher G. Klein and Gregory Howes
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
