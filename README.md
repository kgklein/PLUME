# PLUME (Plasma in a Linear Uniform Magnetized Environment)
Kristopher Klein
kris.klein@gmail.com
Lunar and Planetary Laboratory, University of Arizona

PLUME calculates the hot plasma dispersion relation for a plasma with an arbitrary number of ion and electron species with relative drifts and bi-Maxwellian velocity distributions. The calculation follows Waves in Plasma by Stix, Chapter 10 equatios 66-73.

This code uses an F90 adaptation (waves.f90 by Greg Howes) of the Hot Plasma Dispersion Relation originally by Eliot Quataert, extending thesolver from a proton-electron plasma with Maxwellian velocity distributions to the generalized case of n components with biMaxwellian velocity distributions and arbitrary parallel drift velocities.

# JET-PLUME (Judging Energy Transfer in a Plasma in a Linear Uniform Magnetized Environment) 
Collin Brown
collin.crbrown@gmail.com
University of Iowa

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique. 

## Other Contributing Authors
Greg Howes
Eliot Quataert
Jason TenBarge

## PLUME/JET-PLUME Setup
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

The second way is to use the jupyter notebook

```
jupyter notebook
```

Then run 'examplelinfpc.ipynb'. The jupyter notebook wrapper will assist in making inputs, output folders, loading output, and plotting output. It is recommended for most use cases of PLUME or JET-PLUME.

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

Time is normalized to the reference cyclotron velocity,
```
Omega_ref = q_ref B/m_ref c
```

space is normalized to the reference gyroradius,
```
rho_ref = w_perp,ref/Omega_ref
```

and the relative drift speeds are normalized to the Alfven velocity calculated using the mass and density of the reference species,
```
v_A,ref = B/sqrt{4 pi n_ref m_ref}
```


## PLUME Routines

The main program (plume.f90) executes different subroutines from the disprels.f90 module depending on the input value of 'option'.

PLUME has two main procedures:
+Find roots of the dispersion relation within a given region of complex frequency space (omega_r,gamma)/Omega_ref
(Using the map_search routine with use_map=.true.)

+Calculate the dispersion relation as a function of varying parameters
(Using the om_scan routine with the map_search or an input guess providing the initial values for the frequencies)


The code will either find 'nroot_max' roots of the dispersion relation for the input parameters, or refine 'nroot_max' input guesses for such roots. The (4 + nspec * 6 ) parameters can be varied, with particular solutions being followed for the variations. 

## JET-PLUME Routines

TODO!

## PLUME/JET-PLUME Input parameters
Input parameters are specified in *.in files and organized by namelist (see example_map_par.in). Here, we describe the inputs of each namelist.

Params:
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

GLOBAL PLASMA PARAMETERS
>betap, kperp, kpar, vtp

used for calculation of dispersion relation

NUMBER OF PLASMA SPECIES
>nspec

determines how many %species_n namelists to read in

NUMBER OF PARAMETER SCANS
>nscan

determines how many %scan_n namelists to read in

CASE SELECTION
> option

determines the type of calculation(s) to be performed with cases outlined in README.cases (TODO: move to below)

NUMBER OF DISPERSION ROOTS TO FOLLOW
> nroot_max

determines max number of roots to follow (i.e. follow gradient and attempt to find root at)

USE MAPPING SUBROUTINE
> use_map

if use_map==.false., read in nroot_max guess_n namelists
if use_map==.true., map subroutine outputs nroot_max roots of the system code

OUTPUT/SUPPRESS OUTPUT TO SCREEN
writeOut

LOCATION OF DATA OUTPUT
dataName=
data is output into folder 'data/_dataName_/' 
**Please make sure this folder exists before running PLUME or JET-PLUME**

outputName=
name of data that is put into above folder

header for multiparameter map (i.e. option 2) (TODO: check what this means- is it possible to have single parameter map? Thought plasma had to be neutral?)
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
nspec of these namelists will be read in

SPECIES DEPENDENT PLASMA PARAMETERS
>tauS,muS,alphaS,Qs,Ds,vvS
used for calculation of dispersion relation
NOTE: v drift_proton = 0.0 
(Calculations are done in the proton rest frame)

**Make sure to insure that the plasma is:**
QUASINEUTRAL
>Sum_s [n_s q_s] = 0
and 
CURRENT-FREE
>Sum_s [n_s q_s V_s] = 0

(Will output warnings, but will not halt calculation if either of these condition fail)

Map scan inputs
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
if use_map==T, determines bounds of map scan
_omi_, _omf_ upper and lower bounds for real frequency
_gami_, _gamf_ upper and lower bounds for imaginary frequency
_loggridw_,_loggridg_ determines log (T) or linear (F) spacing for real (w) or imaginary (g) scan

Parameter scan inputs
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
Typically, only one scan_input is used, but it is possible for nscan namelists to be read in, each of which determine the nature of the plasma parameters scans

scan_type, scan_style;
```
Defines nature of parameter scans:
Style: -1- Global two component Scan:
     Type: 0 k_0-> k_1           
           1 theta_0 -> theta_1
           2 k_fixed angle
Style: 0- Global Scan:
     Type: 0 kperp
           1 kpar
           2 betap
           3 vtp
Style: 1-nspec-> Species Scan:
    Type: 0 tau_s
           1 mu_s
           2 alph_s
           3 Q_s
           4 D_s
           5 vv_s
```

PARAMETER RANGE
_swi_, _swf_ are the initial, final values of parameter range
(EXCEPTIONS
style_s=-1 are scans with multiple components
type_s=0: scan from current value of (kperp,kpar) to (kperp,kpar)=(swi,swf)
type_s=1: theta scan from current value of (k,theta) to (k,theta)=(k,swi), with swi in degrees
type_s=2: k scan from current value of (kperp,kpar) to (k=swf) with constant theta=atan(kperp/kpar))

When making a map of two parameters (om_double_scan, called with option=2) do the theta scan before the k_fixed angle scan.

_swlog_ is log (T) or linear (F) scan of variables

_ns_ is the number of output steps along scan

_nres_ is the number of steps between output steps


heating, eigen determine if heating eigenfunctions/eigenvalues are calculated and output eigenfunctions and wavepower


OUTPUT FORMAT (TODO: UPDATE TO INCLUDE NEW OUTPUT FORMAT FROM DR. HOWES' NEW POWER SPLIT)
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

The value of the dispersion relation is written to the file:
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

Thus, the index of U_s onward will differ depending on the
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

**E**,**B**,**Us**,**ns**, and **Ps** are all calculated in the routine calc_eigen.

_Ex_,_y_,_z_ are found using linear algebra in manipulation of the equation (TODO: latex below if possible)

```
(wave vector tensor = lam) . E = 0
```
We take Ex = (1 + 0i), which yields
```
Ez/Ex = -(lam 21 lam 32 - lam 31 lam 22)/(lam 23 lam 32 - lam 33 lam 22)
Ey/Ex = -(Ez/Ex lam 33 - lam 31)/lam 32
```

We can express _Bx_,_y_,_z_ using Faraday's law
```
nabla x E = -(1/c)partial B/ partial t
```

Rearrangement and expressing in our dimensionless units yields
```
(1/Ex)(Bx    (1/(om vtp sqrt(alph_p) Ex)) (-kpar Ey
       By  =                                kpar Ex - kperp Ez
       Bz)                                  kperp Ey)
```
For the velocity fluctuations, we use the expression for the first order current fluctuations
```
j_s = -i omega Chi_s /(4 pi) . E = n_s q_s V_s
```

We choose a normalization by v_Ap, yielding
```
V_s/v_Ap = -i om (Qs/Ds)(vtp^2 sqrt(alph_p/beta_p)) Chi_s . E/Ex
```

Density is extracted from the linearized continuity equation
```
partial n/partial t + nabla . (Us ns)
```
```
omega (n0 + ns) - k . [(U0+Us)(n0+ns)] = 0
```
Taking the 1st order contribution
```
ns/n0 = (k . Us/v_Ap)/(om - kpar vv_s/sqrt(beta_p alph_p)) * (1/sqrt(beta_p alph_p))
```

The heating calculation for species s is described in _Stix 1992, pg 289_ and _Quatart 1998_

The calculation of (PS) the energy per wave period injected into or extracted from species s per unit wave energy (W) is valid for om >> gam, and definately falls apart for om = 0 
(which happens with some frequency for the Alfven and Entropy modes) (TODO: update this eq below?)
```
Ps = (E* . Chi_s^a |gam = 0 . E)/4W
```

```
Chi_s^a = (Chi_s - Chi_s*)/(2i)
```
-=--=-=-=-=- 
We have now added in the gyrokinetic dispersion relation!
