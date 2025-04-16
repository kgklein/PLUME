# PLUME Output

PLUME writes output solutions to the `/data` directory, to the subdirectory specified by `dataName`.
All output file names will have the `outputName` string included in the name to distinguish between distinct calculations. 

## dispersion_*outputName*.map

Value of the dispersion tensor $\mathcal{\Lambda}(\omega_{\textrm{r}},\gamma)$ on a defined complex frequency grid.  
Solutions to the dispersion relation satisfy $|\mathcal{\Lambda}|  =0$.
This file is generated from the *map_search* subroutine in *disprels.f90*, and invoked when `use_map` =.true. 

The data is ordered in columns as  
1. $\omega_r$  
2. $\gamma$   
3. $\log_{10} |\mathcal{\Lambda}|$  
4. Re $[|\mathcal{\Lambda}|]$  
5. Im $[|\mathcal{\Lambda}|]$  

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

## *outputName*_*param*_*range*.modeN

The dispersion relation calculated from `om_scan` will have different formats based upon the number of species as well as the suplementary calculations invoked.

The *param* name will be assigned based on the kind of parameter scan performed, specifically:

*`scan_style`=-1:*

- `k` for scan from $\textbf{k}_0 \rho_{ref}$ to $\textbf{k}_1 \rho_{ref}$ (`scan_type`=0).
- `theta` for scan from $\theta_0$ to $\theta_1$ with fixed $|k|\rho_{ref}$ (`scan_type`=1).
- `k_fixed_theta` for scan from $|k|_1\rho_{ref}$ to $|k|_1\rho_{ref}$ with a fixed value of $\theta$ (`scan_type`=2).

*`scan_style`=0:*

- `kperp` for scan of $k_\perp \rho_{ref}$ (`scan_type`=0).
- `kpar` for scan of $k_\parallel \rho_{ref}$ (`scan_type`=1).
- `beta` for scan of $\beta_{ref,\parallel}$ (`scan_type`=2).
- `vtp` for scan of $v_{t,ref,\parallel}/c$ (`scan_type`=3).

*`scan_style`=1-`nspec`:*

- `tauS` for scan of $T_{ref,\parallel}/T_{j,\parallel}$ (`scan_type`=0).
- `muS` for scan of $m_{ref}/m_{j}$ (`scan_type`=1).
- `alphS` for scan of $T_{j,\perp}/T_{j,\parallel}$ (`scan_type`=2).
- `qs` for scan of $q_{ref}/q_{j}$ (`scan_type`=3).
- `Ds` for scan of $n_{j}/n_{ref}$ (`scan_type`=4).
- `Vs` for scan of $v_{j,drift}/v_{A,ref}$ (`scan_type`=5).

The initial and final values of the parameter scanned over will be reflected in the *range* portion of the file name.

One file for each mode followed will be produced.

The first six columns are the same regardless of supplemental calculations.
They are ordered as
1. $k_\perp \rho_{ref}$
2. $k_\parallel \rho_{ref}$
3. $\beta_{ref,\parallel}$
4. $v_{t,ref,\parallel}/c$
5. $\omega_{\textrm{r}}/\Omega_{ref}$   
6. $\gamma/\Omega_{ref}$

If `eigen` is set to true, the next set of columns will be the eigenfluctuations, with fields normalized to $E_x$, velocity normalized to $c E_x/B_0$ and, density normalized to $n_{0s} E_x/B_0$, 

7. Re $[B_x]$   
8. Im $[B_x]$   
9. Re $[B_y]$   
10. Im $[B_y]$   
11. Re $[B_z]$   
12. Im $[B_z]$   
13. Re $[E_x]$   
14. Im $[E_x]$   
15. Re $[E_y]$   
16. Im $[E_y]$   
17. Re $[E_z]$   
18. Im $[E_z]$   
19. [+6(j-1)] Re $[\delta U_{x,j}]$   
20. [+6(j-1)] Im $[\delta U_{x,j}]$   
21. [+6(j-1)] Re $[\delta U_{y,j}]$   
22. [+6(j-1)] Im $[\delta U_{y,j}]$   
23. [+6(j-1)] Re $[\delta U_{z,j}]$   
24. [+6(j-1)] Im $[\delta U_{z,j}]$   
19. [+6(`nspec`)+2(j-1)] Re $[\delta n_{j}]$   
20. [+6(`nspec`)+2(j-1)] Im $[\delta n_{j}]$ 

where `j` ranges from 1 to `nspec`.
The normalization follows Eqns. X-Y in Klein et al RNAAS 2025 [doi here]

If `heat` is set to true, the next set of columns will be the power absorption or emission from each component. If `new_low_n` is set to true, additional terms associated with Landau, Transit time, and Cyclotron heating will be output. If `eigen` is false, this data will start in the 7th column. If eigen is true, this data will start in the 18+8 `nspec`+1st column.

If `new_low_n` is true, we have

- 18+(8-6{!eigen})`nspec`+j. $P_j$
- 18+(9-6{!eigen})`nspec`+j. $P_j^{yy}$ (Transit Time Damping term 1).
- 18+(10-6{!eigen})`nspec`+j. $P_j^{yz}$ (Transit Time Damping term 2).
- 18+(11-6{!eigen})`nspec`+j. $P_j^{zy}$ (Landau Damping term 1).
- 18+(12-6{!eigen})`nspec`+j. $P_j^{zz}$ (Landau Damping term 2).
- 18+(13-6{!eigen})`nspec`+j. $P_j^{n=0}$ (sum of Landau and Transit Time Damping).
- 18+(14-6{!eigen})`nspec`+j. $P_j^{n=\pm 1}$ ($n=\pm 1$ Cyclotron Damping).
Here, {!eigen} (negation of eigen boolean) is equal to 1 if eigen is false and 0 if eigen is true.

This will be followed by 6`nspec` columns containing the parameter lists $\mathcal{P}_j$ for each species or component. 
- 19+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $T_{ref,\parallel}/T_{j,\parallel}$.
- 20+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $m_{ref}/m_{j}$.
- 21+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $T_{j,\perp}/T_{j,\parallel}$.
- 22+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $q_{ref}/q_{j}$.
- 23+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $n_{j}/n_{ref}$.
- 24+`noutperspec` `nspec`+6(j-1)-12({!eigen}). $v_{j,drift}/v_{A,ref}$.
Here, noutperspec = 0 is the number of additional outputs created by setting heating or eigen to true. If eigen and heating are false, then noutperspec = 7, if eigen is false and heating is true, then noutperspec = 8, and if eigen is true and if eigen and heating are true, then noutperspec=15. Note that {!eigen} (negation of eigen boolean) is equal to 1 if eigen is false and 0 if eigen is true.

This same data structure is preserved for the output from `om_double_scan`.
The file naming convention will include a value of *param* from both of the two parameters scanned, and the code will not output information about the range of parameters in the file name.
