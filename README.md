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

<img src="https://github.com/kgklein/PLUME/blob/linfpc/Jet-Plume_Logo.svg" width=50% height=50%>

# 1.) PLUME/JET-PLUME Setup
GFortran is needed to compile PLUME/JET-PLUME. To install, see [here](https://gcc.gnu.org/wiki/GFortranBinaries) or run 'xcode-select --install' in terminal if using macOS. (Alternatively, one may use IFort if appropriate changes are made to the Makefile).

To compile, first open a terminal, navigate to the desired folder and clone this repo:
```
git clone https://github.com/kgklein/PLUME
```
Open the newly downloaded folder:
```
cd PLUME
```
Clean the old build:
```
make clean
```
Compile:
```
Make
```

From here, there are two ways to use PLUME or JET-PLUME. The first way is with the command line:

## Command Line

First, make an output folder (folder name is found in input file)
```
mkdir data/*datafoldername*
```
Run:
```
./plume.e *path*/*to*/inputputflnm.in
```
This will output data into the 'data/*datafoldername*' directory as text.

## Jupyter Notebook Wrapper

The second way is to use the jupyter notebook (to install, see [here](https://jupyter.org/install))

```
jupyter notebook
```

Then run 'examplelinfpc.ipynb'. The jupyter notebook wrapper will assist in making inputs, output folders, loading output, and plotting output. It is recommended for most use cases of PLUME or JET-PLUME. The notebook will provide examples of all key functionalities of the wrapper.

<br />
<br />
<br />
<br />
<br />

# 2.) PLUME/JET-PLUME Normalization

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

Time is normalized to the reference cyclotron velocity, $\Omega_{r} = q_{r} B/m_{r} c$, space is normalized to the reference gyroradius, $\rho_{r} = w_{\perp,r}/\Omega_{r}$ and the relative drift speeds are normalized to the Alfven velocity calculated using the mass and density of the reference species, $v_{A,r} = B/\sqrt{4 \pi n_{r} m_{r}}$.

<br />
<br />
<br />
<br />
<br />

# 3.) PLUME Routines

The main program (plume.f90) executes different subroutines from the disprels.f90 module depending on the input value of 'option'.

PLUME has two main procedures:

1.) Find roots of the dispersion relation within a given region of complex frequency space ($\omega_r$,$\gamma$)/$\Omega_{r}$, where $\Omega_{r} is the reference cyclotron frequency.
(Using the map_search routine with use_map=.true.).

2.) Calculate the dispersion relation as a function of varying parameters
(Using the om_scan routine with the map_search or an input guess providing the initial values for the frequencies).

The code will either find 'nroot_max' roots of the dispersion relation for the input parameters, or refine 'nroot_max' input guesses for such roots. The (4 + _nspec_ * 6 ) parameters can be varied, with particular solutions being followed for the variations. 

# JET-PLUME Routines

JET-PLUME is a key subroutine of PLUME that is called by specifying the correct input value of 'option'. JET-PLUME has two procedures:

1.) compute the perturbed distribution function, $f_{s,1}(\mathbf{v},\mathbf{k},\omega)$ (real and imaginary parts), velocity-space energy transfer for all 3 components of $\mathbf{E}$ (i.e. $E_i$ $\in$  { $E_x$, $E_y$, $E_z$ }), $C_{E_i}(\mathbf{v})$, for specified $k_{||}$, and $k_{\perp}$ in 3D cartesian-space coordinates. If $\omega$ is specified, the quantites will be computed for each provided $\omega$. If $\omega$ is not specified, a PLUME subroutine will be called to attempt too compute up to the specified number of roots in the system, and the above quantities will be computed for each.

(Note without loss of generality, we choose coordiates such that $k_{\perp,1} = k_{\perp}$, thus the $\perp,1$ direction is in the plane of the wave and magnetic field.)

2.) Same as one, except in 2D 'gyro coordinates', i.e. $f_{s,1}(v_{||},v_{\perp},\mathbf{k},\omega)$, $C_{E_i}(\mathbf{v})(v_{||},v_{\perp},\mathbf{k},\omega)$, equal to integrating out the third coorindate, $\theta$, in cylindrical coordinates in velocity space.

# 4.) PLUME/JET-PLUME OPTIONS
The following values correspond to the value that should be passed to the 'option' input parameter to specify which routine is ran.

## PLUME routines
0: Calculate Roots for input plasma parameters.

1: {Calculate Roots for input plasma parameters<br />
                   OR<br />
   Read in guesses for frequency values, and refine}<br />
   Scan over 'nscan' parameters, with range and type specified in *.in file, Outputing 'nscan'x'nroot_max' files, each calculating the dispersion relation for mode n along for the variation of a given parameter

2: {Calculate Roots for input plasma parameters<br />
                   OR<br />
   Read in guesses for frequency values, and refine}<br />
   Scan over two parameters, with range and type specified in *.in file. Produces maps of dispersion relations in (parameter 1, parameter 2) space of nroot_max modes. (**Make sure to set nscan = 2 for this option and to define 2 scan_input namelists**.)

3: Replicating SAGA scan from Gullveig (the precursor of this code). A hardwired scan of (k, theta) at a particular value of (betap, alph_p). This is used primarily for testing/ verification.

```
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
```

4: Calculate map of complex frequency space for a scan of plasma parameters.

5: Find roots for parameters along a prescribed path.

## JET-PLUME routines

6: Compute the field-particle correlation and fs1 in gyrotropic coordinates (e.g. $C_{E_i}(v_{\perp},v_{par})$ ) for found roots if use_map is True and specified roots if False and roots are provided.

7: Compute the field-particle correlation adn fs1 in cartesian coordinates (e.g. $C_{E_i}(v_{x},v_{y},v_{z})$ ) for found roots if use_map is True and specified roots if False and roots are provided.

<br />
<br />
<br />
<br />
<br />

# 5.) PLUME/JET-PLUME Input parameters
Input parameters are specified in *.in files and organized by namelist (see example_map_par.in). Here, we describe the inputs of each namelist. First, we provide the format of the namelist, and then we break down each input.

## Global Parameters:
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

**nspec**<br />
Number of plasma species<br />
Determines how many '%species_n' namelists to read in

**nscan**<br />
Number of parameter scans<br />
Determines how many %scan_n namelists to read in

**option**<br />
Case selection<br />
determines the type of calculation(s) to be performed with cases outlined in cases section (see below).

**nroot_max**<br />
Number of dispersion roots to follow<br />
determines max number of roots to follow (i.e. follow gradient and attempt to find root at)

**use_map**<br />
Use mapping subroutine
if use_map==.false., read in nroot_max guess_n namelists<br />
if use_map==.true., map subroutine outputs nroot_max roots of the system code

**writeOut**<br />
Output or suppress output to screen

**dataName**<br />
Location of data output<br />
data is output into folder 'data/_dataName_/'<br />
_Please make sure this folder exists before running PLUME or JET-PLUME_

**outputName**<br />
Data name<br />
name of data that is put into above folder

## Species parameters
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

**tauS,muS,alphS,Qs,Ds,vvS**<br />
Species dependent plasma parameters<br />
Theses are species parameters for calculation of dispersion relation. 'tauS' is the ratio of the reference parallel temperature to the species parallel temperature, $\tau_s = T_{||,r}/T_{||,s}$. 'muS' is the mass ratio of the reference species to the spcies mass, $\mu_s = m_r/m_s$. 'alphS' is the ratio of the perpendicular species temperature to the parallel species temperature, $\aleph_s = T_{\perp,s}/T_{||,s}$. 'Qs' is the ratio of the reference species charge to the species charge, $Q_s = q_r/q_s$. 'Ds' is the ratio of the species number density to the refence species number density, $D_s = n_s/n_r$. 'vvS' is the species drift velocity normalized to the reference Alfven speed, $\overline{V_s} = V_s / V_{a,r}$. 

_Make sure to insure that the plasma is:_<br />
Quasineutral
```math
\sum_s n_s q_s = 0
```
and current free
```math
\sum_s n_s q_s V_s = 0
```
(NOTE: The code will only output warnings, but it will not halt calculation if either of these condition fail.)

## Map scan inputs
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

**loggridw,loggridg**<br />
Log spacing<br />
'loggridw', 'loggridg' determines log (T) or linear (F) spacing for real (w) or imaginary (g) scan

**omi,omf,gami,gamf**
Grid bounds<br />
These variable determine the bounds of the grid that the code will use as start points to try and find roots from. <br />
'omi', 'omf' are the upper and lower bounds for real frequency <br />
'gami', 'gamf' are upper and lower bounds for imaginary frequency <br />

## Parameter scan inputs
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

**scan_type, scan_style**<br />
Type and style of scan<br />
Defines nature of parameter scans. Can scan over either global or specific species parameter.<br />
```
Style: -1- Global two component Scan:
     Type: 0 k_0-> k_1
           1 theta_0 -> theta_1
           2 k_fixed angle
Style: 0- Global Scan:
     Type: 0 kperp
           1 kpar
           2 betap
           3 vtp
Style: 1-nspec (1 for reference, 'n' for nth species) -> Species Parameter Scan:<br />
    Type: 0 tau_s
           1 mu_s
           2 alph_s
           3 Q_s
           4 D_s
           5 vv_s
```

**swi, swf**<br />
Parameter range<br />
'swi'/'swf' are the initial/final values of parameter range for selected paramter<br />
EXCEPTIONS:<br />
style_s=-1 are scans with multiple components and **changes the above definition of 'swi'/'swf'** to the below definition<br />
type_s=0: scan from current value of (kperp,kpar) to (kperp,kpar)=(swi,swf)<br />
type_s=1: theta scan from current value of (k,theta) to (k,theta)=(k,swi), with swi in degrees (note: theta = atan(kperp/kpar))<br />
type_s=2: k scan from current value of (kperp,kpar) to (k=swf) with constant theta=atan(kperp/kpar)<br />
(When making a map of two parameters (om_double_scan, called with option=2) do the theta scan before the k_fixed angle scan.)

**swlog**<br />
Scan spacing<br />
'swlog' is true for log spaced scan or false for linear spaced scan.

**ns**<br />
Scan number of steps<br />
'ns' is the number of output steps along scan.

**nres**<br />
Scan output steps<br />
'nres' is the number of steps between output steps.

**heating, eigen**<br />
Scan outputs<br />
'heating'/'eigen' determine if heating eigenfunctions/eigenvalues are calculated/output.

## Guess inputs
Examples guess input namelists (There should be additional namelists for 'nroot_max' number of guess with names of the form 'guess_n'):
```
&guess_1
g_om=1.E-1
g_gam=-1.E-1
/
```

**g_om,g_gam**<br />
Guess values<br />
'g_om', 'g_gam' values for real and imaginary frequency guesses. Note that g_gam<0 corresponds to damped wave in the sign convention used by PLUME.<br />
nroot_max namelist files read in (if map_true == F)

## JET-PLUME/FPC inputs
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

**vperpmin,vperpmax,vparmin,vparmax**<br />
Gyro velocity space grid bounds<br />
These correspond to the minimum and maximum values of $v_{\perp}$ and $v_{||}$ computed when computing the FPC in gyro coordinates. These values are ignored when using the cartesian routine.

**vxmin,vxmax,vymin,vymax,vzmin,vzmax**<br />
Cartesian velocity space grid bounds<br />
These correspond to the minimum and maximum values of $v_{x}$, $v_{y}$ and $v_{z}$ computed when computing the FPC in cartesian coordinates. These values are ignored when using the gyro routine.

**delv**<br />
Velocity grid spacing<br />
'delv' is the spacing in between grid points in both routines. Spacing is equal in all directions.

<br />
<br />
<br />
<br />
<br />

# 6.) OUTPUT FORMAT 

In this section, we describe the output format for each routine.

## Dispersion Relation/Eigenfunction Output
This output is created when options 1, 2, 3, or 5 is selected.

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
where 'param' is the name of the parameter being scanned. The file consists of many rows for each sample of the dispersion relation over some specified sweep parameter and columns for each parameter. The quantity each column corresponds to is specified below. Each row is self consistent, meaning that all values in that row correspond to a specific mode.

The power of each species by each mode can be computed using two different ways by changes the value of .new_low_n. in vars.f90 and _recompiling the code_. New users are recommended to keep new_low_n = .true.. This is the default format.

### Dispersion Relation/Eigenfunction Output
The format for .new_low_n. = .false. and is shown below (eigen = .false. omits the Bx,y,z Ex,y,z Ux,y,z output, heating = .false. omits the P output, and (heating = .false. or (new_low_n = .false. and .low_n. = .false.)) omits P_split output (note .new_low_n. takes priority if .low_n. and .new_low_n. are both true):

Note, when low_n is selected,
```
P_split|s = {psld,psttd,psn0,pscd}
```
where psld is the power due to landua damping, psttd is the power due to transit time damping, psn0 is the power due to n=0 modes, and pscd is the power due to cylcotron damping (all terms are for each species individually).

When new_low_n is selected,
```
P_split|s = {psttd1,psttd2,psld1,psld2,psn0,pscd}
```
where, p1ld1 is the off diagonal $\chi_zy$ term, psld2 is the on diagonal $\chi_zz$ term, pdttd1 is the off diagonal $\chi_yy$ term, psttd2 is the $\chi_yz$ term, psn0 is the power due to n=0 modes, and pscd is the power due to cylcotron damping (all terms are for each species individually).

The output format depends on the number of species in the plasma. (The output for species 1 is always first, species 2 is second, etc.)
For two species:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30
  n|species }|(if eigen==true) !31-32,33-34
{ P|species }|(if heat==true) !35,36
{Ps_split|species}(if heat==true and low_n||new_low_n=true) !37-40,41-44 or 37-42,43-48
(tau,mu,alph,q,D,vv)|species   !37-42,43-48 or 49-54,55-60
```

For three species:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36
  n|species }|(if eigen==true) !37-38,39-40,41-42
{ P|species }|(if heat ==true) !43,44,45
{Ps_split|species}(if heat==true and low_n||new_low_n=true) !46-49,50-53,54-57 or 46-51,52-57,58-63
(tau,mu,alph,q,D,vv)|species   !46-51,52-57,58-63 or 64-69,70-75,76-81
```

For four species:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36,37-42
  n|species }|(if eigen==true) !43-44,45-46,47-48,49-50
{ P|species }|(if heat ==true) !51,52,53,54
{Ps_split|species}(if heat==true and low_n||new_low_n=true)!55-58,59-62,63-66,67-70 or 55-60,61-66,67-72,73-78,79-84
(tau,mu,alph,q,D,vv)|species   !71-76,77-82,83-88,89-94 or 85-90,91-96,97-102,103-108
```

For five species:
```
kperp,kpar,betap,vtp, !1-4
omega,gamma,          !5-6
{ Bx,y,z              !7-12
  Ex,y,z              !13-18
  Ux,y,z|species      !19-24,25-30,31-36,37-42, 43-48
  n|species }|(if eigen==true) !49-50,51-52,53-54,55-56,57-58
{ P|species }|(if heat ==true) !59,60,61,62,63
{Ps_split|species}(if heat==true and low_n||new_low_n=true)!64-67,68-71,72-75,76-79,80-83 or 64-69,70-75,76-81,82-87,88-93
(tau,mu,alph,q,D,vv)|species   !84-89,90-95,96-101,102-107,108-113 or 94-99,100-105,106-111,112-117,118-123
```

The above pattern holds for _n_ species.

Note, the index of U_s onward will differ depending on the
number of plasma species and if eigen and heating are turned on.

## Root Map Output
This output is created when options 4 is selected. The value of the roots map is written to the files:
```
          'data/',trim(dataName),'/',&
          dispersion_*outputName*.map
```
AND
```
          'data/',trim(dataName),'/',&
          dispersion_*outputName*.roots
```

### dispersion_*outputName*.map
WIP...

```
index real, index imaginary, omega, log10(val(ir,ii)) !1-5
sign(1.,real(dal(ir,ii)))*log10(1.+abs(real(dal(ir,ii)))),&
sign(1.,aimag(dal(ir,ii)))*log10(1.+abs(aimag(dal(ir,ii))))
```
for all ir,ii in the grid where omega is the computed complex frequency, 'dal' is the value of the dispersion relation (i.e. the determinant of the matrix), and 'val' is the absolute value of the value fo the dispersion relation.

### dispersion_*outputName*.roots
The output format depends on the number of species in the plasma. (The output for species 1 is always first, species 2 is second, etc.)
For _n_ species:
```
kperp,kpar,betap,vtp, !1-4
{omega,gamma|roots} !{5+2*i_root,6+2*i_root}
(tau,mu,alph,q,D,vv)|species !{(7+2*n_roots+6*i_species,12+2*n_roots+6*i_species)}    
```

*outputName*.roots contains the list of the found roots. Note, that roots will sometimes be found outside of the original specified domain. Here, all current estimated roots are assumed to have at least approximately converged and will be output here. 

## JET-PLUME / FPC Output
These outputs are created when options 6 or 7 is selected.

Running JET-PLUME will create measurements of $C_{E_i}$, and $f_{s1}$ (the fourier coefficients of the perturbed portion of the distribution function), on the selected velocity space grid projection.

### Gyro Output
These outputs are created when option 6 is selected.

The value of  $C_{E_i}$, and $f_{s1}$ on a projected cartesian velocity grid is written to the files:
```
'data/',trim(dataName),'/',&
*outputName*.cpar.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.cperp.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.df1gyro.real.specie*num*.mode*num* 
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.df1gyro.imag.specie*num*.mode*num* 
```
for each species and each mode (i.e. selected $\mathbf{k}$ $\omega$ solution).

The files contain the correlation with respect to $E_{||}$, the correlation with respect to $E_{\perp}$, the real part of the perturbed distribution function fourier coefficients, and the imaginary part of the perturbed distribution function fourier coefficients respectively.

Each file contains rows and columns where the ith, jth element corresponds to the value on the ith, jth position on the velocity grid. $v_{||}$ evolves along the columns and $v_{\perp}$ evolves along the rows. The 0th,0th element (i.e. top left) corresponds to the minimum value of $v_{||}$ and $v_{\perp}$.

Each grid point is spaced by delv, including vmax and vmin. Warning!!!: If vmax and vmin bounds are not integer multiples of delv, a passive warning message will be output, but the code will continue running. This may cause off by 1 errors when loading.

### Cartesian Ouput
These outputs are created when option 7 is selected.

The value of  $C_{E_i}$, and $f_{s1}$ on a projected cartesian velocity grid is written to the files:
```
'data/',trim(dataName),'/',&
*outputName*.cparcart.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.cperp1.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.cperp2.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.df1.real.specie*num*.mode*num*
```
AND
```
'data/',trim(dataName),'/',&
*outputName*.df1.imag.specie*num*.mode*num*
```
for each species and each mode (i.e. selected $\mathbf{k}$ $\omega$ solution).

The files contain the correlation with respect to $E_{||}$, the correlation with respect to $E_{\perp,1}$, the correlation with respect to $E_{\perp,2}$, the real part of the perturbed distribution function fourier coefficients, and the imaginary part of the perturbed distribution function fourier coefficients respectively.

Each file contains all three projections of each quantity, starting with the $(v_{\perp,1},v_{\perp,2})$, $(v_{\perp,1},v_{||})$, $(v_{\perp,2},v_{||})$, separted by '---'.

For the $(v_{\perp,1},v_{\perp,2})$ projection, $v_{\perp,1}$ evolves with rows and $v_{\perp,2}$ evolves with columns. The 0th,0th element (i.e. top left) corresponds to the minimum value of $v_{\perp,1}$ and $v_{\perp,2}$.

For the $(v_{\perp,1},v_{||})$ projection, $v_{\perp,1}$ evolves with rows and $v_{||}$ evolves with columns. The 0th,0th element (i.e. top left) corresponds to the minimum value of $v_{\perp,1}$ and $v_{\perp,2}$.

For the $(v_{\perp,2},v_{||})$ projection, $v_{||}$ evolves with rows and $v_{\perp,2}$ evolves with columns. The 0th,0th element (i.e. top left) corresponds to the minimum value of $v_{\perp,1}$ and $v_{\perp,2}$.

Each grid point is spaced by delv, including vmax and vmin. Warning!!!: If vmax and vmin bounds are not integer multiples of delv, a passive warning message will be output, but the code will continue running. This may cause off by 1 errors when loading.
<br />
<br />
<br />
<br />
<br />

# 7.) Eigenfunction Calculation

The eigenfunctions $\mathbf{E}$, $\mathbf{B}$, $\mathbf{U_s}$, $\mathbf{n_s}$, and $\mathbf{P_s}$ are all calculated in the routine calc_eigen. Each eigenfunction, $A(\omega,\mathbf{k})$, is the complex Fourier coeffienct containing information about the amplitude and phase of the linear response of each quantity to linear perturbation of the incident mode.

$E_x$, $E_y$, $E_z$ are found using linear algebra in manipulation of the equation. (Note, it is common to write $A_x = A_{\perp,1}$, $A_y = A_{\perp,2}$, $A_z = A_{\perp,||}$ as done below. Here, $A_{||}$ is the component that is parallel to the external magnetic field, $A_{\perp,1}$ is the component in the plane of the incident wave and the external magentic field, $A_{\perp,2}$ is in the component normal to the plane of the incident wave and the external magnetic field.)

```math
\lambda \cdot E = 0
```
where $\lambda$ is the wave vector tensor, defined as the matrix such that $\lambda \mathbf{E} \equiv \mathbf{n} \times {n} \times \mathbf{E} + \epsilon \cdot \mathbf{E} = 0$,  where $\mathbf{n} = \mathbf{k}c/\omega$ is the refractive index and the right expression is the wave vector equation, and $\epsilon = \mathbf{I} + \sum \chi_s$ is the dielectric tensor with species suseptibilities $\chi_s$, computed in Chapter 10 of Stix (1992).

We take $E_x = (1 + 0i)$, which yields
```math
E_z/E_x = -(\lambda_{21} \lambda_{32} - \lambda_{31} \lambda_{22})/(\lambda_{23} \lambda_{32} - \lambda_{33} \lambda_{22})
```
and
```math
E_y/E_x = -(Ez/Ex \lambda_{33} - \lambda_{31})/\lambda_{32}
```

We can express $B_x$, $B_y$, $B_z$ using Faraday's law,
```math
\nabla \times \mathbf{E} = -(1/c) \partial B/ \partial t
```

Taking the eletric field fluctuations $\mathbf{E}$, we can next express the magnetic field fluctuations $\delta B$ using Faraday's Law, using our dimensionless scale
```math
\delta \mathbf{B} = \frac{\mathbf{k} \rho_r \times \delta \mathbf{E}}{\bar{\omega}\hat{w}_r \sqrt{\aleph_r}}
```

For the velocity fluctuations, we use the expression for the first order current fluctuations
```math
j_s = -i \omega \chi_s /(4 pi) /cdot E = n_s q_s V_s
```

We choose to normalize by v_Ap (the ion Alfven velocity), yielding
```math
V_s/v_Ap = -i \omega (Q_s/D_s)(v_{tp}^2 \sqrt{\aleph_p/\beta_p}) \chi_s \cdot \mathbf{E}/E_x
```

Density is extracted from the linearized continuity equation.
```math
\partial n/\partial t + \nabla \cdot (\mathbf{U}_s n_s) = 0 
```
```math
\omega (n_0 + n_s) - \mathbf{k} \cdot [(\mathbf{U}_0+\mathbf{U}_s)(n_0+n_s)] = 0
```
Taking the 1st order contribution
```math
n_s/n_0 = (\mathbf{k} \cdot \mathbf{U}_s/v_Ap)/(\omega - k_{||} \mathbf{V}_s/\sqrt{\beta_p \aleph_p}) * (1/\sqrt{\beta_p \aleph_p})
```

The heating calculation can be done in two ways, dependent on the value of '.new_low_n.' in vars.f90 (true by default).

if '.new_low_n.' is true, the heating calculation is as follows:

For landau damping,
```math
\lim_{\gamma/\omega \rightarrow 0} P_{LD,s} =  \frac{i \omega}{16 \pi} \big[ (\chi_{zz,s}^{(n=0)}-\chi_{zz,s}^{(n=0)*}E_{||}E_{||}^*+\chi_{zy,s}^{(n=0)E_{\perp,2}E_{||}^*}-\chi_{zy,s}^{(n=0)*}E_{\perp,2}^*E_{||})\big]
```
and for transit time damping, 
```math
\lim_{\gamma/\omega \rightarrow 0} P_{TTD,s} =  \frac{i \omega}{16 \pi} \big[ (\chi_{yy,s}^{(n=0)}-\chi_{yy,s}^{(n=0)*}E_{\perp,2}E_{\perp,2}^*+\chi_{yz,s}^{(n=0)E_{\perp,2}E_{||}^*}-\chi_{yz,s}^{(n=0)*}E_{\perp,2}^*E_{||})\big]
```
and for cyclotron damping,
```math
P_s^{\rm CD,n'}= \frac{\omega}{8 \pi}
\bigg[|E_{\perp,1}|^2\left(\chi_s^{xx} -\chi_s^{xx,*}\right)
+|E_{\perp,2}|^2\left(\chi_s^{yy} -\chi_s^{yy,*}\right)
+\left(E_{\perp,1}^*E_{\perp,2} - E_{\perp,2}^*E_{\perp,1}\right)\left(\chi_s^{xy} -\chi_s^{yx,*}\right)
\bigg]_{n=n',\omega=\omega_{\rm real}}
```

If '.new_low_n.' is false, the heating calculation for species s is described in _Stix 1992, pg 289_ and _Quatart 1998_. The calculation of (PS) the energy per wave period injected into or extracted from species s per unit wave energy (W) is valid for om >> gam, and definately falls apart for om = 0 (which happens with some frequency for the Alfven and Entropy modes) and takes the form
```math
P_s = (E^* \cdot \chi_s^a |_{\gamma = 0} \cdot E)/4W
```
where
```math
\chi_s^a = (\chi_s - \chi_s^*)/(2i).
```

<br />
<br />
<br />
<br />
<br />

# 8.) License
TODO!!!

# Papers

TODO: link relevant papers of PLUME/JET-PLUME. Be sure to note our upcoming paper whose appendix goes into the most depth about PLUME
