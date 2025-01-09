# PLUME Output

PLUME writes output solutions to the `/data` directory, to the subdirectory specified by `dataName`.
All output file names will have the `outputName` string included in the name to distinguish between distinct calculations. 

## dispersion_*outputName*.map

Value of the dispersion tensor $\mathcal{D}(\omega_{\textrm{r}},\gamma)$ on a defined complex frequency grid.  
Solutions to the dispersion relation satisfy $|\mathcal{D}|  =0$.
This file is generated from the *map_search* subroutine in *disprels.f90*, and invoked when `use_map` =.true. 

The data is ordered in columns as  
1. $\omega_r$  
2. $\gamma$   
3. $\log_{10} |\mathcal{D}|$  
4. Re $[|\mathcal{D}|]$  
5. Im $[|\mathcal{D}|]$  

The *&maps* namelist in *filename*.in determines the structure of *filename*.map.  
The range of $\omega_{\textrm{r}}/\Omega_p$ is from `omi` to `omi` with `nr` steps. Logorithmic or linear spacing is selected with `loggridw`. *`nr` and `ni` are currently hardcoded. Will add as user options.*
The range of $\gamma_{\textrm{r}}/\Omega_p$ is from `gami` to `gami` with `ni` steps. Logorithmic or linear spacing is selected with `loggridg`.

## dispersion_*outputName*.roots

Identified solutions to the dispersion relation $|\mathcal{D}|  =0$, calculated using *refine_guess* in *disprels.f90*.

The data is ordered in columns as:
1. $k_\perp \rho_{ref}$
2. $k_\parallel \rho_{ref}$
3. $\beta_{ref,\parallel}$
4. $v_{t,ref,\parallel}/c$
5. $\omega_{r}/\Omega_{ref}$
6. $\gamma/\Omega_{ref}$

This will be followed by 6`nspec` columns containing the parameter lists $\mathcal{P}_j$ for each species or component.
- $T_{ref,\parallel}/T_{j,\parallel}$.
- $m_{ref}/m_{j}$.
- $T_{j,\perp}/T_{j,\parallel}$.
- $q_{ref}/q_{j}$.
- $n_{j}/n_{ref}$.
- $v_{j,drift}/v_{A,ref}$.

Only the first `nroot_max` solutions will be identified and written to file.

## *filename*.scan_*scan_type_l*.root_m

The complex frequencies associated with solution *m* calculated from `om_scan`  
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_p$
2. $k_\parallel d_p$  
3. $\omega_{\textrm{r}}/\Omega_p$   
4. $\gamma/\Omega_p$   

See the *&scan_input* namelist description in the Quick Guide for details on determining the kind of wavevector scan.  
This same data structure is preserved for the output from `om_double_scan`.

## *filename*.eigen_*scan_type_l*.root_m

The eigenfunctions associated with solution *m* calculated from `om_scan` when `eigen` is set to .true.
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_p$
2. $k_\parallel d_p$  
3. $\omega_{\textrm{r}}/\Omega_p$   
4. $\gamma/\Omega_p$   
5. Re $[E_x]$ 
6. Im $[E_x]$ 
7. Re $[E_y]$ 
8. Im $[E_y]$ 
9. Re $[E_z]$ 
10. Im $[E_z]$ 
11. Re $[B_x]$ 
12. Im $[B_x]$ 
13. Re $[B_y]$ 
14. Im $[B_y]$ 
15. Re $[B_z]$ 
16. Im $[B_z]$  
17. [+6(is-1)] Re $[\delta U_{x,is}]$   
18. [+6(is-1)] Im $[\delta U_{x,is}]$   
19. [+6(is-1)] Re $[\delta U_{y,is}]$   
20. [+6(is-1)] Im $[\delta U_{y,is}]$   
21. [+6(is-1)] Re $[\delta U_{z,is}]$   
22. [+6(is-1)] Im $[\delta U_{z,is}]$   
17. [+6(`nspec`)+2(is-1)] Re $[\delta n_{is}]$   
18. [+6(`nspec`)+2(is-1)] Im $[\delta n_{is}]$ 

This same data structure is preserved for the output from `om_double_scan`.

## *filename*.heat_*scan_type_l*.root_m

The heating rates associated with solution *m* calculated from `om_scan` when `heating` is set to .true.
in the *&scan_input_l* namelist.

The data is ordered as  
1. $k_\perp d_p$
2. $k_\parallel d_p$  
3. $\omega_{\textrm{r}}/\Omega_p$   
4. $\gamma/\Omega_p$   
5. [+(is-1)] $\gamma_{is}/\omega$

This same data structure is preserved for the output from `om_double_scan`.
