file_map='../../data/example/dispersion_map_par.map'
file_roots='../../data/example/dispersion_map_par.map'

set contour
set pm3d
unset surface
set view map
set cntrparam levels 50

set xlabel 'omega_r/Omega_p'
set ylabel 'gamma/Omega_p'

set title 'Lambda[(k_{perp} rho_p=k_{par} rho_p = 0.001)]'

splot \
file_map u 1:2:3 w l title ''

pause -1