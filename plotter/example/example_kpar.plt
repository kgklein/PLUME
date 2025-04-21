file_in(nn)=sprintf(\
'../../data/example/map_par_kpar_1_10000.mode%d',nn)


set log
set grid

set xlabel 'k_{par} rho_p'
set ylabel '[omega_r/Omega_p,|gamma|/Omega_p]'

set key top left

set title 'Dispersion Relations for (k_{perp} rho_p = 0.001)'

plot \
file_in(1) u 2:5 w l lc rgb '#000000' title 'Alfven, omega_r' ,\
file_in(1) u 2:(-$6) w l lc rgb '#000000' dt 2 title 'Alfven, |gamma|' ,\
file_in(2) u 2:5 w l lc rgb '#aa0000' title 'Fast, omega_r' ,\
file_in(2) u 2:(-$6) w l lc rgb '#aa0000' dt 2 title 'Fast, |gamma|'

pause -1