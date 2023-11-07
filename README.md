# PLUME (Plasma in a Linear Uniform Magnetized Environment)
Kristopher Klein<br />
kris.klein@gmail.com<br />
Lunar and Planetary Laboratory, University of Arizona<br />

PLUME calculates the hot plasma dispersion relation for a plasma with an arbitrary number of ion and electron species with relative drifts and bi-Maxwellian velocity distributions. The calculation follows Waves in Plasma by Stix, Chapter 10 equations 66-73.

This code uses an F90 adaptation (waves.f90 by Greg Howes) of the Hot Plasma Dispersion Relation originally by Eliot Quataert, extending thesolver from a proton-electron plasma with Maxwellian velocity distributions to the generalized case of n components with biMaxwellian velocity distributions and arbitrary parallel drift velocities.

# JET-PLUME (Judging Energy Transfer in a Plasma in a Linear Uniform Magnetized Environment) 
Collin Brown<br />
collin.crbrown@gmail.com<br />
University of Iowa<br />

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique. 

## Other Contributing Authors
Greg Howes<br />
Eliot Quataert<br />
Jason TenBarge<br />

## PLUME/JET-PLUME Setup
GFortran is needed to compile PLUME/JET-PLUME. To install, see [here](https://gcc.gnu.org/wiki/GFortranBinaries). (Alternatively, one may use IFort if appropriate changes are made to the Makefile).

To compile and run, first open a terminal, navigate to the desired folder and clone this repo:
```
git clone https://github.com/kgklein/PLUME/tree/linfpc
```
Open the newly downloaded folder:
```
cd PLUME
```
Clean any old builds:
```
make clean
```
Compile:
```
Make
```

From here, there are two ways to use PLUME or JET-PLUME. The first way is with the command line:

### Command Line

Make output folder (folder name is found in input file)
```
mkdir data/*datafoldername*
```
Run:
```
./plume.e *path*/*to*/inputputflnm.in
```
This will output data into the 'data/*datafoldername*' directory as text.

### Jupyter Notebook Wrapper

The second way is to use the jupyter notebook (to install, see [here](https://jupyter.org/install))

```
jupyter notebook
```

Then run 'examplelinfpc.ipynb'. The jupyter notebook wrapper will assist in making inputs, output folders, loading output, and plotting output. It is recommended for most use cases of PLUME or JET-PLUME. The notebook will provide examples of all key functionalities of the wrapper.

## PLUME/JET-PLUME Normalization

The Dispersion relation for omega/Omega_ref is dependent on four global parameters:
```
     betap: Plasma Reference Beta:                            8 pi n_ref T_parallel ref /B^2
     kperp: Perpendicular wavelength:                         k perp rho_ref
     kpar : Parallel wavelength:                              k parallel rho_ref
     vtp  : Normalized parallel proton thermal velocity:      sqrt(2 T_parallel,ref/m_ref c^2)
```

and six species dependent parameters:
```
     tau_s  : T_ref/T_s
     mu_s   : m_ref/m_s
     alph_s : T_perp/T_parallel|s
     Q_s    : q_ref/q_s
     D_s    : n_s/n_ref
     vv_s   : v_s,drift/v_Aref
```
The values for these parameters are extracted from *.in file, appended after the executable program call.

Time is normalized to the reference cyclotron velocity, $\Omega_{r} = q_{r} B/m_{r} c$, space is normalized to the reference gyroradius, $\rho_{r} = w_{\perp,ref}/\Omega_{r}$ and the relative drift speeds are normalized to the Alfven velocity calculated using the mass and density of the reference species, $v_{A,ref} = B/\sqrt{4 \pi n_{r} m_{r}}$.

## PLUME Routines

The main program (plume.f90) executes different subroutines from the disprels.f90 module depending on the input value of 'option'.

PLUME has two main procedures:

1.) Find roots of the dispersion relation within a given region of complex frequency space (omega_r,gamma)/Omega_ref
(Using the map_search routine with use_map=.true.)

2.) Calculate the dispersion relation as a function of varying parameters
(Using the om_scan routine with the map_search or an input guess providing the initial values for the frequencies)

The code will either find 'nroot_max' roots of the dispersion relation for the input parameters, or refine 'nroot_max' input guesses for such roots. The (4 + nspec * 6 ) parameters can be varied, with particular solutions being followed for the variations. 

## JET-PLUME Routines

JET-PLUME is a key subroutine of PLUME that is called by specifying the correct input value of 'option'. JET-PLUME has two procedures:

1.) compute the perturbed distribution function, $f_{s,1}(\mathbf{v},\mathbf{k},\omega)$ (real and imaginary parts), velocity-space energy transfer for all 3 components of $\mathbf{E}$ (i.e. $E_i$ $\in$  { $E_x$, $E_y$, $E_z$ }), $C_{E_i}(\mathbf{v})$, for specified $k_{||}$, and $k_{\perp}$ in 3D cartesian-space coordinates. If $\omega$ is specified, the quantites will be computed for each provided $\omega$. If $\omega$ is not specified, a PLUME subroutine will be called to attempt too compute up to the specified number of roots in the system, and the above quantities will be computed for each.

(Note without loss of generality, we choose coordiates such that $k_{\perp,1} = k_{\perp}$, thus the $\perp,1$ direction is in the plane of the wave and magnetic field.)

2.) Same as one, except in 2D 'gyro coordinates', i.e. $f_{s,1}(v_{||},v_{\perp},\mathbf{k},\omega)$, $C_{E_i}(\mathbf{v})(v_{||},v_{\perp},\mathbf{k},\omega)$, equal to integrating out the third coorindate, $\theta$, in cylindrical coordinates in velocity space.

## PLUME/JET-PLUME Input parameters
Input parameters are specified in *.in files and organized by namelist (see example_map_par.in). Here, we describe the inputs of each namelist. First, we provide the format of the namelist, and then we break down each input.

### Global Parameters:
Example global parameter namelist:
```
&params
betap=1.00
kperp=1.E-3
kpar=1.E-3
vtp=1.E-4
nspec=2
nscan=1
option=1
nroot_max=1
use_map=.false.
writeOut=.true.
dataName='test'
outputName='b1a1'
/
```

Global plasma parameters:
```
betap
kperp
kpar
vtp
```
used for calculation of dispersion relation. 'betap' is $\beta_{||}$, the total parallel plasma beta, 'kperp' is $k_{\perp}$, the total perpendicular wavenumber, 'kpar' is $k_{||}$, the parallel wavenumber, 'vtp' is $v_{th,r}/c$, the ratio of the reference species' (typically protons) thermal velocity to the speed of light.

Number of plasma species:
```
nspec
```
determines how many %species_n namelists to read in

Number of parameter scans:
```
nscan
```
determines how many %scan_n namelists to read in

Case selection
```
option
```
determines the type of calculation(s) to be performed with cases outlined in cases section (see below).

Number of dispersion roots to follow
```
nroot_max
```
determines max number of roots to follow (i.e. follow gradient and attempt to find root at)

Use mapping subroutine
```
use_map
```
if use_map==.false., read in nroot_max guess_n namelists<br />
if use_map==.true., map subroutine outputs nroot_max roots of the system code<br />

Output or suppress output to screen
```
writeOut
```

Location of data output
```
dataName=
```
data is output into folder 'data/_dataName_/' 
**Please make sure this folder exists before running PLUME or JET-PLUME**

Data name
```
outputName=
```
name of data that is put into above folder<br />

###Species parameters
Example namelist for species input (Note: 'species_1' is the reference species and is typically the proton species. There should be additional namelists for 'nspec' number of species with names of the form 'species_n'):
```
&species_1
tauS=1.0
muS=1.0
alphS=1.0
Qs=1.0
Ds=1.0
vvS=0.0
/
```

Species dependent plasma parameters
```
tauS
muS
alphS
Qs
Ds
vvS
```
Theses are species parameters for calculation of dispersion relation. 'tauS' is the ratio of the reference parallel temperature to the species parallel temperature, $\tau_s = T_{||,r}/T_{||,s}$. 'muS' is the mass ratio of the reference species to the spcies mass, $\mu_s = m_r/m_s$. 'alphS' is the ratio of the perpendicular species temperature to the parallel species temperature, $\aleph_s = T_{\perp,s}/T_{||,s}$. 'Qs' is the ratio of the reference species charge to the species charge, $Q_s = q_r/q_s$. 'Ds' is the ratio of the species number density to the refence species number density, $D_s = n_s/n_r$. 'vvS' is the species drift velocity normalized to the reference Alfven speed, $\bar{V}_s = V_s/V_{a,r}$. 

NOTE: v drift_reference should be 0.0 (Calculations are done in the proton rest frame)

**Make sure to insure that the plasma is:**
Quasineutral
```math
\sum_s n_s q_s = 0
```
and current free
```math
\sum_s n_s q_s V_s = 0
```
(NOTE: The code will only output warnings, but it will not halt calculation if either of these condition fail.)

###Map scan inputs
Example map scan inputs namelist (only used if use_map==T):
```
&maps
loggridw=.false.	
loggridg=.false.
omi     =-1.5E-3
omf     = 1.5E-3
gami    =-1.E-3
gamf    = 1.E-4
/
```
Log spacing
```
loggridw,loggridg
```
'loggridw', 'loggridg' determines log (T) or linear (F) spacing for real (w) or imaginary (g) scan

Grid bounds
```
omi,omf
gami,gamf
```
These variable determine the bounds of the grid that the code will use as start points to try and find roots from. <br />
'omi', 'omf' are the upper and lower bounds for real frequency <br />
'gami', 'gamf' are upper and lower bounds for imaginary frequency <br />

###Parameter scan inputs
Example namelist for a parameter scan (typically, only one scan_input is used, but it is possible for nscan namelists to be read in, each of which determine the nature of the plasma parameters scans):
```
&scan_input_1
scan_type=0
scan_style=0
swi=1.
swf=2.
swlog=.true.
ns=16
nres=16
heating=.true.
eigen=.true.
/
```

Type and style of scan
```
scan_type, scan_style;
```
Defines nature of parameter scans. Can scan over either global or specific species parameter.<br />
Style: -1- Global two component Scan:<br />
     Type: 0 k_0-> k_1<br />
           1 theta_0 -> theta_1<br />
           2 k_fixed angle<br />
Style: 0- Global Scan:<br />
     Type: 0 kperp<br />
           1 kpar<br />
           2 betap<br />
           3 vtp<br />
Style: 1-nspec (1 for reference, 'n' for nth species) -> Species Parameter Scan:<br />
    Type: 0 tau_s<br />
           1 mu_s<br />
           2 alph_s<br />
           3 Q_s<br />
           4 D_s<br />
           5 vv_s<br />

Parameter range
```
swi, swf 
```
'swi'/'swf' are the initial/final values of parameter range for selected paramter<br />
EXCEPTIONS:<br />
style_s=-1 are scans with multiple components<br />
type_s=0: scan from current value of (kperp,kpar) to (kperp,kpar)=(swi,swf)<br />
type_s=1: theta scan from current value of (k,theta) to (k,theta)=(k,swi), with swi in degrees<br />
type_s=2: k scan from current value of (kperp,kpar) to (k=swf) with constant theta=atan(kperp/kpar)<br />
<br />
When making a map of two parameters (om_double_scan, called with option=2) do the theta scan before the k_fixed angle scan.

Scan spacing
```
swlog
```
'swlog' is true for log spaced scan or false for linear spaced scan.

Scan number of steps
```
ns
```
'ns' is the number of output steps along scan.

Scan output steps
```
nres
```
'nres' is the number of steps between output steps.

Scan outputs
```
heating, eigen
```
'heating'/'eigen' determine if heating eigenfunctions/eigenvalues are calculated/output.

### JET-PLUME/FPC inputs
Example JET-PLUME/FPC namelist:
```
&fpc
vperpmin=0.01
vperpmax=3.0
vparmin=-3.0
vparmax=3.0
delv=0.15
vxmin=-3.0
vxmax=3.0
vymin=-3.0
vymax=3.0
vzmin=-3.0
vzmax=3.0
/
```

Gyro velocity space grid bounds
```
vperpmin,vperpmax,vparmin,vparmax
```
These correspond to the minimum and maximum values of $v_{\perp}$ and $v_{||}$ computed when computing the FPC in gyro coordinates. These values are ignored when using the cartesian routine.

Cartesian velocity space grid bounds
```
vxmin,vxmax,vymin,vymax,vzmin,vzmax
```
These correspond to the minimum and maximum values of $v_{x}$, $v_{y}$ and $v_{z}$ computed when computing the FPC in cartesian coordinates. These values are ignored when using the gyro routine.

Velocity grid spacing 
```
delv
```
'delv' is the spacing in between grid points in both routines. Spacing is equal in all directions.

## PLUME/JET-PLUME OPTIONS
The following values correspond to the value that should be passed to the 'option' input parameter to specify which routine is ran

### PLUME routines
0: Calculate Roots for input plasma parameters.

1: {Calculate Roots for input plasma parameters
    OR
   Read in guesses for frequency values, and refine}
   Scan over 'nscan' parameters, with range and type specified in *.in file
   Outputing 'nscan'x'nroot_max' files, each calculating the dispersion relation
   for mode n along for the variation of a given parameter

2: {Calculate Roots for input plasma parameters
    OR
   Read in guesses for frequency values, and refine}
   Scan over two parameters, with range and type specified in *.in file
   Produces maps of dispersion relations in (parameter 1, parameter 2) space
   of nroot_max modes

   !MAKE SURE TO SET NSCAN=2 FOR THIS OPTION

3:   !Replicating 
     !SAGA scan
     !from Gullveig (the precursor of this code)
     !A hardwired scan of 
     !  (k, theta) 
     !at a particular value of
     !  (betap, alph_p)

     !start at 
     !(beta,alph_p,kperp, kpar)=
     !(1.0, 1.0,   1.E-3,1.E-3)    

     !scan_1, scan_2 should be over theta, k_fixed, at desired resolution
     !scan_3, scan_4, scan_5  (beta_p, alph_p, (k1,k2))
     !-=-=-=-=-=
     !scan_1: theta (desired resolution) range_i:theta_final
              !Style: -1; Type: 1
     !scan_2: k_fixed angle (desired resolution) range_f:
              !Style: -1; Type: 2
     !scan_3: beta_p (1.0 -> desired beta_p)
              !Style:  0; Type: 2
     !scan_4: alph_p (1.0 -> desired alph_p)
              !Style:  1; Type: 2
     !scan_5: alph_p (1.E-3,1.E-3) -> 
              !(1.E-3 sin(theta_initial), 1.E-3 cos(theta_initial))
              !(range_i,range_f)
              !Style:  0; Type: 0     

!Initial Roots
&guess_1
g_om=9.9973E-04   
g_gam=-2.2572E-10 
/

&guess_2
g_om=2.0304E-03   
g_gam=-5.4273E-05
/

&guess_3
g_om=0.0000E+00   
g_gam=-7.2110E-04
/

&guess_4
g_om=1.1830E-03
g_gam=-7.3333E-04
/

4: Calculate map of complex frequency space for a scan of plasma parameters.

5: Find roots for parameters along a prescribed path.

### JET-PLUME routines

6: Compute FPC, fs1, in gyro coordinates ($C_{E_i} (v_{\perp},v_{par})$) for found roots if use_map is True and specified roots if False and roots are provided.

7: Compute FPC, fs1, in cartesian coordinates ($C_{E_i} (v_{\perp},v_{par})$) for found roots if use_map is True and specified roots if False and roots are provided.

## OUTPUT FORMAT 
(TODO: UPDATE TO INCLUDE NEW OUTPUT FORMAT FROM DR. HOWES' NEW POWER SPLIT!!!)

The power of each species by each mode can be computed using two different ways by changes the value of TODO in vars.f90 and recompiling the code.

The format for TODO is shown below:

The output format for parameter depends on the number of species in the plasma
FOR EACH MODE:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30
  n|species }|(if eigen==true) !31-32,33-34
{ P|species }|(if heat ==true) !35,36
(tau,mu,alph,q,D,vv)|species   !37-42,43-48
```

for three components:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36
  n|species }|(if eigen==true) !37-38,39-40,41-42
{ P|species }|(if heat ==true) !43,44,45
(tau,mu,alph,q,D,vv)|species   !46-51,52-57,58-63
```

for four components:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36,37-42
  n|species }|(if eigen==true) !43-44,45-46,47-48,49-50
{ P|species }|(if heat ==true) !51,52,53,54
(tau,mu,alph,q,D,vv)|species   !55-60,61-66,67-72,73-78
```

for five components:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36,37-42, 43-48
  n|species }|(if eigen==true) !49-50,51-52,53-54,55-56,57-58
{ P|species }|(if heat ==true) !59,60,61,62,63
(tau,mu,alph,q,D,vv)|species   !64-69,70-75,76-81,82-87,88-93
```

The value of the dispersion relation is written to the files:
```
          'data/',trim(dataName),'/',&
          trim(param),'_',int(1000.*scan(is)%range_i),&
          '_',int(1000.*scan(is)%range_f)
```
OR
```
          'data/',trim(dataName),'/',&
          trim(param),'_s',scan(is)%style_s,'_',int(1000.*scan(is)%range_i),&
          '_',int(1000.*scan(is)%range_f)
```
where 'param' is the name of the parameter being scanned. 

Note, the index of U_s onward will differ depending on the
number of plasma species and if eigen and heating are turned on.


Guess input
```
&guess_1
g_om=1.E-1
g_gam=-1.E-1
/
```

nroot_max namelist files read in (if map_true == F)
(g_om, g_gam) values for real and imaginary frequency guesses

## Eigenfunction Calculation:

The eigenfunctions $\mathbf{E}$,$\mathbf{B}$,$\mathbf{U_s}$,$\mathbf{n_s}$, and $\mathbf{P_s}$ are all calculated in the routine calc_eigen.

$E_x$, $E_y$, $E_z$ are found using linear algebra in manipulation of the equation. Note, it is common to write $E_x = E_{\perp,1}$, $E_y = E_{\perp,2}$, $E_z = E_{\perp,||}$ as done below.

```math
\lambda \cdot E = 0
```
where $\lambda$ is the wave vector tensor, defined as the matrix such that $\lambda \mathbf{E} \equiv \mathbf{n} \times {n} \times \mathbf{E} + \epsilon \cdot \mathbf{E} = 0$,  where $\mathbf{n} = \mathbf{k}c/\omega$ is the refractive index and the right expression is the wave vector equation, and $\epsilon = \mathbf{I} + \sum \chi_s$ is the dielectric tensor with species suseptibilities $\chi_s$, computed in Chapter 10 of Stix (1992).

We take $E_x = (1 + 0i)$, which yields
```math
E_z/E_x = -(\lambda_{21} \lambda_{32} - \lambda_{31} \lambda_{22})/(\lambda_{23} \lambda_{32} - \lambda_{33} \lambda_{22})
```
```math
E_y/E_x = -(Ez/Ex \lambda_{33} - \lambda_{31})/\lambda_{32}
```

We can express $B_x$,$B_y$,$B_z$ using Faraday's law,
```math
\nabla \times \mathbf{E} = -(1/c) \partial B/ \partial t
```

Rearrangement and expressing in our dimensionless units yields
```math
(1/E_x$)(B_x)    (1/(om vtp sqrt(alph_p) Ex)) (-kpar Ey
       By  =                                kpar Ex - kperp Ez
       Bz)                                  kperp Ey)
```

For the velocity fluctuations, we use the expression for the first order current fluctuations
```math
j_s = -i \omega \chi_s /(4 pi) . E = n_s q_s V_s
```

We choose a normalization by v_Ap, yielding (TODO: we might need to update this equation)
```math
V_s/v_Ap = -i om (Q_s/D_s)(v_{tp}^2 sqrt(\alpeh_p/\beta_p)) \chi_s . \mathbf{E}/E_x
```

Density is extracted from the linearized continuity equation (TODO: check below?)
```math
\partial n/\partial t + \nabla \cdot (U_s n_s)
```
```math
\omega (n_0 + n_s) - \mathbf{k} \cdot [(U_0+U_s)(n_0+n_s)] = 0
```
Taking the 1st order contribution (TODO: rewrite vv_s)
```math
n_s/n_0 = (\mathbf{k} \cdot U_s/v_Ap)/(\omega - k_{||} vv_s/sqrt(\beta_p \aleph_p)) * (1/\sqrt{\beta_p \aleph_p})
```

The heating calculation for species s is described in _Stix 1992, pg 289_ and _Quatart 1998_

The calculation of (PS) the energy per wave period injected into or extracted from species s per unit wave energy (W) is valid for om >> gam, and definately falls apart for om = 0 
(which happens with some frequency for the Alfven and Entropy modes) (TODO: update this eq below!! (include both forms!!)
```math
P_s = (E^* \cdot \chi_s^a |_{\gamma = 0} \cdot E)/4W
```

```math
\chi_s^a = (Chi_s - Chi_s^*)/(2i)
```
(TODO: add commentary  to the above two equations)
