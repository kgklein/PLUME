# PLUME Input

This is a reference for the key input parameters used by PLUME (currently editing).

## Namelists in input files.

The following namelists and associated input parameters are read in by PLUME from the input file.

### *&params*  
General system parameters.

**`betap`**  
Initial reference parallel plasma beta $\beta_{ref,\parallel} = 8 \pi n_{ref} T_{ref,\parallel}/B^2$.

**`kperp`**  
Initial perpendicular wavevector $k_{\perp} \rho_{ref}$.

**`kpar`**  
Initial parallel wavevector $k_{\parallel} \rho_{ref}$.

**`vtp`**  
Initial reference parallel thermal velocity, normalized to c $v_{t,ref,\parallel}/c$.

**`nspec`**   
Number of plasma species or components.

**`nscan`**   
Number parameter or wavevector scans to execute.

**`option`**   
Determines set of scans to perform. Choice of:

- -1: Calculate the value of the dispersion relation
      at a single $(\omega/\Omega_{ref}, \gamma/\Omega_{ref})$
-  0: Calculate Roots for input plasma parameters.
-  1: Calculate Roots for input plasma parameters or Reads in root value
      and then scan over plasma parameters, with range and type specified in *.in file.
-  2: Calculate Roots for input plasma parameters or Reads in root value
      and then scan over two-dimensional plasma parameters space,
      with range and type specified in *.in file.
-  3: Deprecated...
-  4: Make multiple maps of complex frequency space.
-  5: Find roots for parameters along a prescribed path
      Path is set by solar wind models, with values calculated and
      output by helper function (in development, the radial scan function.).

**`nroot_max`**   
Number of dispersion solutions to find and follow.

**`use_map`**   
Choice of:  

- True: Searching for roots over a map in complex frequency space (see &maps namelist).  
- False: Input `nroot_max` guesses for solutions (see &guess_1 namelist).

**`low_n`**  
Logical to toggle on or off outputing the $n=0$ and $\pm 1$ resonances for
characterizing heating channels. *As made redundant by `new_low_n`, slated from removal.*

**`new_low_n`**  
*New* Logical to toggle on or off outputing the $n=0$ and $\pm 1$ resonances for
characterizing heating channels; correctly separates Landau, Transit, and Cyclotron damping.

**`writeOut`**  
Write or suppress output to screen.

**`dataName`**  
Subdirectory (below PLUME/data/) where outputs will be written.

**`outputName`**
String to be included in output files to identify specific run.

### *&guess_m*  
Initial guess of complex frequency for $m$th solution.  
Only used when `use_map`=.false.  
Need to have number of namelists equal to `nroot_max`.

**`g_om`**  
Guess for real solution $\omega_{r}/\Omega_{ref} $.

**`g_gam`**  
Guess for imaginary solution $\gamma/\Omega_{ref} $.


### *&maps*  
Range of complex frequencies for map_scan subroutine.  
Only used when `use_map`=.true.

**`loggridw`**  
Linear (F) or Log (T) spacing for $\omega_{r}/\Omega_{p}$ map search.
Spacing automatically calculated between `omi` and `omf`.  

**`loggridg`**  
Linear (F) or Log (T) spacing for $\gamma/\Omega_{p}$ map search.
Spacing automatically calculated between `gami` and `gamf`  

**`omi`**  
Smallest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`omf`**  
Largest $\omega_{r}/\Omega_{p}$ value for complex map search.

**`gami`**      
Smallest $\gamma/\Omega_{p}$ value for complex map search.

**`gamf`**  
Largest $\gamma/\Omega_{p}$ value for complex map search.

**`nr`**  
*TO BE ADDED* Number of $\omega_{r}/\Omega_{p}$ points in frequency grid.

**`ni`**  
*TO BE ADDED* Number of $\gamma/\Omega_{p}$ points in frequency grid.

### *&species_j*  
Species or component parameters list for distribution $f_{j}$.

**`tauS`**  
Relative parallel temperature $T_{ref,\parallel}/T_{j,\parallel}$.

**`muS`**  
Relative mass $m_{ref}/m_{j}$.

**`alphS`**  
Temperature anisotropy $T_{j,\perp}/T_{j,\parallel}$.

**`Qs`**  
Relative charge $q_{ref}/q_{j}$.

**`Ds`**  
Relative density $n_{j}/n_{ref}$.

**`vvS`**  
Drift speed parallel to mean magnetic field, normalized to the reference Alfven velocity  $v_{j,drift}/v_{A,ref}$, where $v_{A,ref}=B/\sqrt{4 \pi n_{ref} m_{ref}}$


### *&scan_input_l*
Inputs for scanning parameter space for $l$th scan.  

**`scan_type`**  
Type of parameter scan.

For `scan_style`=-1 (Global Two-Component Scan), options of:
- 0: Scan from $\textbf{k}_0 \rho_{ref}$ to $\textbf{k}_1 \rho_{ref}$. Scans from current value of $(k_\perp,k_\parallel) \rho_{ref}$ to $k_\perp \rho_p$=`swi` and $k_\parallel \rho_{ref}$=`swf`.
- 1: Scan from $\theta_0$ to $\theta_1$ with fixed $|k|\rho_{ref}$. Scans from current value of $(|k|\rho_{ref},\theta)$ to $(|k|\rho_{ref},$`swi`$)$, with `swi` in degrees.
- 2: Scan from $|k|_1\rho_{ref}$ to $|k|_1\rho_{ref}$ with a fixed value of $\theta$. Scan from current value of $(k_\perp,k_\parallel) \rho_{ref}$ to $|k|\rho_{ref}$=`swf` with constant $\theta = \atan (k_\perp/k_\parallel)$.

For `scan_style`=0 (Global Scan), options of:

- 0: $k_\perp \rho_{ref}$
- 1: $k_\parallel \rho_{ref}$
- 2: $\beta_{ref,\parallel}$
- 3: $v_{t,ref,\parallel}/c$

`swi` and `swf` represent the initial and final scan values for `scan_style`=0.

For `scan_style`=1-`nspec` (Parameter scan of component `scan_style`), options of:

- 0: $T_{ref,\parallel}/T_{j,\parallel}$.
- 1: $m_{ref}/m_{j}$.
- 2: $T_{j,\perp}/T_{j,\parallel}$.
- 3: $q_{ref}/q_{j}$.
- 4: $n_{j}/n_{ref}$.
- 5: $v_{j,drift}/v_{A,ref}$.

`swi` and `swf` represent the initial and final scan values for `scan_style`=1 through `nspec`.

**`scan_style`**  
Class of parameter scan. Options of:

- -1: Global Two-Component Scan
- 0: Global Scan
- 1-nspec: Single-Parameter Scan

**`swi`**  
Scan variable to define start of scan through parameter space.

**`swf`**  
Scan variable to define end of scan through parameter space.

**`swlog`**  
Use $\log_{10}$ (T) or linear (F) spacing.

**`ns`**  
Number of output scan values.

**`nres`**  
Resolution between output scan values.

**`heating`**  
Calculates heating rates if true.

**`eigen`**  
Calculates eigenfunctions if true.

**`tensor`**  
Outputs susceptibility tensor if true.     
