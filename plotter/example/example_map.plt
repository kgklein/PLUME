load 'test.def'

#I/O Initialization

version=0

if (version==0){
data_dir='../../data/example'
name='dispersion_map_par'
file_in(type)=sprintf('%s/%s.%s',data_dir,name,type)
file_out=sprintf('example_map_v%2.2d.tex',version)
#Plasma Parameters
kperp=0.001
kpar=0.001
beta=1.0
anip=1.0
anie=1.0
tau=1.0

#Complex frequency ranges
om_min=-3E-3
om_max=-om_min
gam_min=-3E-3
gam_max=1.E-3
mult=1E3
}

om_min=om_min*mult
om_max=om_max*mult
gam_min=gam_min*mult
gam_max=gam_max*mult

set output file_out

#We will overplot the identified roots on
#the dispersion contours.

set multiplot

#Set frequency ranges
set xrange [om_min:om_max]
set yrange [gam_min:gam_max]

set size 0.65,0.75
set origin 0.,0.

#Contour Plot Parameters
set contour
set pm3d
unset surface
set view 0.,0.
set cntrparam levels 50
unset clabel
unset colorbox

set xlabel sprintf('$\omega_{\textrm{r}}/\Omega_p \times 10^{%2.1f}$',log10(mult)) \
offset 0,1.
set label sprintf('$\gamma_{\textrm{r}}/\Omega_p \times 10^{%2.1f}$',log10(mult)) \
at graph -0.1,0.5 rotate center front

set label sprintf('\texttt{PLUME} Dispersion Relation') at graph 0.0,1.6

set label sprintf('$k_\parallel \rho_p=%2.3f$',kpar) at graph -0.85,1.4
set label sprintf('$k_\perp \rho_p=%2.3f$',kperp) at graph -0.85,1.15

set label sprintf('$\beta_{\parallel,p}=%2.4f$',beta) at graph -0.1,1.4
set label sprintf('$\frac{T_{\perp}}{T_{\parallel}}|_p=%2.3f$',anip) at graph -0.1,1.15

set label sprintf('$\frac{T_{\parallel,p}}{T_{\parallel,e}}=%2.3f$',tau) at graph 0.65,1.4
set label sprintf('$\frac{T_{\perp,e}}{T_{\parallel,e}}=%2.3f$',anie) at graph 0.65,1.15


splot \
file_in('map') u ($1*mult):($2*mult):3 w l title ''

unset xlabel
unset ylabel
unset contour
unset pm3d
set surface

set pointsize 0.75
unset label
unset xlabel
unset ylabel
set format x ' '
set format y ' '

set size 0.4455,0.511
set origin 0.098,0.1275

unset label

plot \
file_in('roots') u ($5*mult):($6*mult) lc rgb '#dd0000' pt 7 title '' ,\
0. lc rgb '#00aa00' dt 0 title ''

unset label

unset multiplot