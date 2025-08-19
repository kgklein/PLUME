# JET-PLUME (Judging Energy Transfer in a Plasma in a Linear Uniform Magnetized Environment) 
[Collin Brown](https://www.collinrbrown.com)<br />
collin.crbrown@gmail.com<br />
University of Iowa; Naval Research Lab- Plasma Physics Division<br />

JET-PLUME is an extension to PLUME that predicts wave-particle energy transfer in velocity space using the field-particle correlation technique. 

## Other Contributing Authors
[Greg Howes](https://physics.uiowa.edu/people/gregory-g-howes)<br />
[Kris Klein](https://www.lpl.arizona.edu/faculty/kristopher-klein)<br />
[Jason TenBarge](https://plasmacenter.princeton.edu/people/jason-tenbarge)<br />

<img src="./Jet-Plume_Logo.svg" width=50% height=50%>

# 1.) JET-PLUME Setup and Running

The install for JET-PLUME is the same as PLUME. See the installation section in PLUME readme.

Running JET-PLUME is the same as PLUME except with additional inputs. See Running the PLUME code in PLUME readme.

# 2.) JET-PLUME / FPC New Inputs

One needs to *add* the following to a normal PLUME [input](./input.md).

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

*And be sure to select the correct routine.

# 3.) JET-PLUME Routines

JET-PLUME is a key subroutine of PLUME that is called by specifying the correct input value of 'option'. JET-PLUME has two procedures corresponding to routine 6 and 7. They work by doing the following for each root found by the root finder (Note, we typically assume only one root per input in the wrapper but this is not striclty necessary for command line usage):

6: Compute the field-particle correlation and fs1 in gyrotropic coordinates (e.g. $C_{E_i}(v_{\perp},v_{par})$ ) for found roots if use_map is True and specified roots if False and roots are provided. That is compute the perturbed distribution function, $f_{s,1}(\mathbf{v},\mathbf{k},\omega)$ (real and imaginary parts), velocity-space energy transfer for all 3 components of $\mathbf{E}$ (i.e. $E_i$ $\in$  { $E_x$, $E_y$, $E_z$ }), $C_{E_i}(\mathbf{v})$, for specified $k_{||}$, and $k_{\perp}$ in 3D cartesian-space coordinates. If $\omega$ is specified, the quantites will be computed for each provided $\omega$. If $\omega$ is not specified, a PLUME subroutine will be called to attempt too compute up to the specified number of roots in the system, and the above quantities will be computed for each.

(Note without loss of generality, we choose coordiates such that $k_{\perp,1} = k_{\perp}$, thus the $\perp,1$ direction is in the plane of the wave and magnetic field.)

7: Compute the field-particle correlation adn fs1 in cartesian coordinates (e.g. $C_{E_i}(v_{x},v_{y},v_{z})$ ) for found roots if use_map is True and specified roots if False and roots are provided. That is same as one, except in 2D 'gyro coordinates', i.e. $f_{s,1}(v_{||},v_{\perp},\mathbf{k},\omega)$, $C_{E_i}(\mathbf{v})(v_{||},v_{\perp},\mathbf{k},\omega)$, equal to integrating out the third coorindate, $\theta$, in cylindrical coordinates in velocity space.

# 4.) JET-PLUME OUTPUT FORMAT 

One set of these outputs are created per species per root when options 6 or 7 is selected.

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
