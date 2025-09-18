load 'test.def'

#I/O Initialization

version=0

if (version==0){
data_dir='../../data/example'
name='guess_double_theta_k_fixed_theta'
file_in(mode)=sprintf('%s/%s.mode%d',data_dir,name,mode)

file_out=sprintf('example_double_v%2.2d.tex',version)
#Plasma Parameters
beta=1.0
anip=1.0
anie=1.0
tau=1.0

#wavevector frequency ranges
kperp_min=1E-5
kperp_max=1.E1

kpar_min=1E-5
kpar_max=1.E1

k_min=1.E-3
k_max=1.E1

theta_min=0
theta_max=90

}

pi=3.1415926

set contour
set pm3d
unset surface
unset clabel
set view map

set cntrparam levels 50
set log cb
set cbrange [1E-5:1E5]
set format cb '$10^{%L}$'
set log z
set zrange [1E-5:1E5]

set output file_out

#We will plot the data in both (kperp,kpar) and (|k|,theta)

set size 1.15,0.65

set multiplot

set size 0.55,0.65
set origin 0.,0.

#kperp,kpar

set log x
set xrange [kperp_min:kperp_max]
set format x '$10^{%L}$'

set log y
set yrange [kpar_min:kpar_max]
set format y '$10^{%L}$'

set label sprintf('\texttt{PLUME} Alfven Dispersion Relation; scan over $|k|\rho_p$ and $\theta$.') at graph 0.0,5.2


set xlabel '$k_\perp \rho_p$' offset 0,0.75
set ylabel '$k_\parallel \rho_p$'

splot \
file_in(1) u ($1):($2):($5) w l lc rgb '#000000' title ''

unset label
set origin 0.5,0.

#|k|,theta

set log x
set xrange [k_min:k_max]
set format x '$10^{%L}$'

unset log y
set yrange [theta_min:theta_max]
set format y '${%g}$'

set xlabel '$|k| \rho_p$'
set ylabel '$\theta$'


splot \
file_in(1) u (sqrt($1**2+$2**2)):(atan($1/$2)*180/pi):5 w l lc rgb '#000000' title ''


unset multiplot