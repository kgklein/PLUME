load 'test.def'

#I/O Initialization

version=0

if (version==0){
data_dir='../../data/example'
name='map_par_kpar_1_100000'
file_in(nn)=sprintf('%s/%s.mode%d',data_dir,name,nn)
file_out=sprintf('example_kpar_v%2.2d.tex',version)
#Plasma Parameters
kperp=0.001
kpar=0.001
beta=1.0
anip=1.0
anie=1.0
tau=1.0
k_min=1.E-3
k_max=1.E2

#Limits
om_min=2.E-3
om_max=2.E1

g_min=2.E-6
g_max=1.E0

gam_min=2.E-8
gam_max=1.E0

E_min=-1.1E0
E_max=1.1E0
}

#Size and placement of panels
xx=0.5
dx=0.325

yy=0.4
dy=0.225

base_color='#000000'
n0_color='#009e73'
n1_color='#e68f00'

set log x
set xrange [k_min:k_max]

#Helper Functions
amp(x,y)=sqrt(x**2.+y**2.)
amp_six(x,y,a,b,c,d)=sqrt(x**2.+y**2.+a**2.+b**2.+c**2.+d**2.)

Eright(Exr,Exi,Eyr,Eyi)=\
((Exr+Eyi)**2.+(Exi-Eyr)**2.)**(0.5)/(2.)**0.5

Eleft(Exr,Exi,Eyr,Eyi)=\
((Exr-Eyi)**2.+(Exi+Eyr)**2.)**(0.5)/(2.)**0.5

Ep(Exr,Exi,Eyr,Eyi)=\
(Eright(Exr,Exi,Eyr,Eyi)-Eleft(Exr,Exi,Eyr,Eyi))/\
(Eright(Exr,Exi,Eyr,Eyi)+Eleft(Exr,Exi,Eyr,Eyi))

#-=-=-=-=-
#-=-=-=-=-
set origin 0,0
set output file_out

set size 1.5,1.45
set grid
#-=-=-=-=-

set multiplot

set size xx,yy
set key top left sample 0.5

#Mode One
mm=1
kx=0

#Real Frequency
unset label
ky=4
set origin kx*dx,ky*dy

set label sprintf('\texttt{PLUME} Dispersion Relation') at graph 0.05,1.75
set label sprintf('$k_\perp \rho_p=%4.3f$',kperp) at graph 1.25,1.75
set label sprintf('$T_{\perp,p}/T_{\parallel,p}=%2.2f$',anip) at graph 2.,1.65
set label sprintf('$T_{\perp,e}/T_{\parallel,e}=%2.2f$',anie) at graph 2.,1.9
set label sprintf('$T_{\parallel,p}/T_{\parallel,e}=%2.2f$',tau) at graph 2.0,1.4
set label sprintf('$\beta_p=%2.2f$',beta) at graph 1.25,1.4

set label sprintf('Mode %d',mm) at graph 0.05,1.1

set format x ''
set log y
set format y '$10^{%L}$'
set yrange [om_min:om_max]
set label '$|\omega_{\textrm{r}}|/\Omega_p$' at graph -0.3,0.15 rotate by 90

plot \
file_in(mm) u 2:($5) w l lc rgb base_color dt 1 title '$\omega_{\textrm{r}}>0$' ,\
file_in(mm) u 2:(-$5) w l lc rgb base_color dt 2 title '$\omega_{\textrm{r}}<0$'

#Total Damping Rate
unset label
ky=3
set origin kx*dx,ky*dy

set format x ''
set log y
set format y '$10^{%L}$'
set yrange [g_min:g_max]
set label '$|\gamma|/\Omega_p$' at graph -0.3,0.15 rotate by 90

plot \
file_in(mm) u 2:($6) w l lc rgb base_color dt 1 title '$\gamma>0$' ,\
file_in(mm) u 2:(-$6) w l lc rgb base_color dt 2 title '$\gamma<0$'

#Proton Damping Rates
unset label
ky=2
set origin kx*dx,ky*dy

set format x ''
set log y
set format y '$10^{%L}$'
set yrange [gam_min:gam_max]
set label '$|\gamma_p|/|\omega_{\textrm{r}}|$' at graph -0.3,0.15 rotate by 90

plot \
file_in(mm) u 2:(-$35/sgn($5)) w l lc rgb base_color dt 1 title '$\gamma_p>0$' ,\
file_in(mm) u 2:($35/sgn($5)) w l lc rgb base_color dt 2 title '$\gamma_p<0$' ,\
file_in(mm) u 2:(-$41/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($41/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$42/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($42/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Electron Damping Rates
unset label
ky=1
set origin kx*dx,ky*dy

set format x ''
set log y
set format y '$10^{%L}$'
set yrange [gam_min:gam_max]
set label '$|\gamma_e|/|\omega_{\textrm{r}}|$' at graph -0.3,0.15 rotate by 90

plot \
file_in(mm) u 2:(-$36/sgn($5)) w l lc rgb base_color dt 1 title '$\gamma_e>0$' ,\
file_in(mm) u 2:($36/sgn($5)) w l lc rgb base_color dt 2 title '$\gamma_e<0$' ,\
file_in(mm) u 2:(-$47/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($47/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$48/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($48/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Polarization
unset label
ky=0
set origin kx*dx,ky*dy

set format x '$10^{%L}$'

unset log y
set format y '${%g}$'
set yrange [E_min:E_max]

set label '$\mathcal{P}_{E}^{xy}$' at graph -0.3,0.35 rotate by 90
set label '$k_\parallel d_p$' at graph 0.5,-0.25

plot \
file_in(mm) u ($2):(Ep($13,$14,$15,$16)) \
w l lc rgb base_color dt 1 title ''

#Mode Two
mm=2
kx=1

#Real Frequency
unset label
ky=4
set origin kx*dx,ky*dy

set label sprintf('Mode %d',mm) at graph 0.05,1.1

set format x ''
set log y
set format y ''
set yrange [om_min:om_max]

plot \
file_in(mm) u 2:($5) w l lc rgb base_color dt 1 title '$\omega_{\textrm{r}}>0$' ,\
file_in(mm) u 2:(-$5) w l lc rgb base_color dt 2 title '$\omega_{\textrm{r}}<0$'

#Total Damping Rate
unset label
ky=3
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [g_min:g_max]

plot \
file_in(mm) u 2:($6) w l lc rgb base_color dt 1 title '$\gamma>0$' ,\
file_in(mm) u 2:(-$6) w l lc rgb base_color dt 2 title '$\gamma<0$'

#Proton Damping Rates
unset label
ky=2
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$35/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($35/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$41/sgn($5)) w l lc rgb n0_color dt 1 title '$\gamma_p^{n=0}>0$' ,\
file_in(mm) u 2:($41/sgn($5)) w l lc rgb n0_color dt 2 title '$\gamma_p^{n=0}<0$' ,\
file_in(mm) u 2:(-$42/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($42/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Electron Damping Rates
unset label
ky=1
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$36/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($36/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$47/sgn($5)) w l lc rgb n0_color dt 1 title '$\gamma_e^{n=0}>0$' ,\
file_in(mm) u 2:($47/sgn($5)) w l lc rgb n0_color dt 2 title '$\gamma_e^{n=0}<0$' ,\
file_in(mm) u 2:(-$48/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($48/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Polarization
unset label
ky=0
set origin kx*dx,ky*dy

set format x '$10^{%L}$'

unset log y
set format y ''
set yrange [E_min:E_max]

set label '$k_\parallel d_p$' at graph 0.5,-0.25

plot \
file_in(mm) u ($2):(Ep($13,$14,$15,$16)) \
w l lc rgb base_color dt 1 title ''

#Mode Three
mm=3
kx=2

#Real Frequency
unset label
ky=4
set origin kx*dx,ky*dy

set label sprintf('Mode %d',mm) at graph 0.05,1.1

set format x ''
set log y
set format y ''
set yrange [om_min:om_max]

plot \
file_in(mm) u 2:($5) w l lc rgb base_color dt 1 title '$\omega_{\textrm{r}}>0$' ,\
file_in(mm) u 2:(-$5) w l lc rgb base_color dt 2 title '$\omega_{\textrm{r}}<0$'

#Total Damping Rate
unset label
ky=3
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [g_min:g_max]

plot \
file_in(mm) u 2:($6) w l lc rgb base_color dt 1 title '$\gamma>0$' ,\
file_in(mm) u 2:(-$6) w l lc rgb base_color dt 2 title '$\gamma<0$'

#Proton Damping Rates
unset label
ky=2
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$35/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($35/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$41/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($41/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$42/sgn($5)) w l lc rgb n1_color dt 1 title '$\gamma_p^{n=\pm 1}>0$' ,\
file_in(mm) u 2:($42/sgn($5)) w l lc rgb n1_color dt 2 title '$\gamma_p^{n=\pm 1}<0$'

#Electron Damping Rates
unset label
ky=1
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$36/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($36/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$47/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($47/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$48/sgn($5)) w l lc rgb n1_color dt 1 title '$\gamma_e^{n=\pm 1}>0$' ,\
file_in(mm) u 2:($48/sgn($5)) w l lc rgb n1_color dt 2 title '$\gamma_e^{n=\pm 1}<0$'

#Polarization
unset label
ky=0
set origin kx*dx,ky*dy

set format x '$10^{%L}$'

unset log y
set format y ''
set yrange [E_min:E_max]

set label '$k_\parallel d_p$' at graph 0.5,-0.25

plot \
file_in(mm) u ($2):(Ep($13,$14,$15,$16)) \
w l lc rgb base_color dt 1 title ''

#Mode Four
mm=4
kx=3

#Real Frequency
unset label
ky=4
set origin kx*dx,ky*dy

set label sprintf('Mode %d',mm) at graph 0.05,1.1

set format x ''
set log y
set format y ''
set yrange [om_min:om_max]

plot \
file_in(mm) u 2:($5) w l lc rgb base_color dt 1 title '$\omega_{\textrm{r}}>0$' ,\
file_in(mm) u 2:(-$5) w l lc rgb base_color dt 2 title '$\omega_{\textrm{r}}<0$'

#Total Damping Rate
unset label
ky=3
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [g_min:g_max]

plot \
file_in(mm) u 2:($6) w l lc rgb base_color dt 1 title '$\gamma>0$' ,\
file_in(mm) u 2:(-$6) w l lc rgb base_color dt 2 title '$\gamma<0$'

#Proton Damping Rates
unset label
ky=2
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$35/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($35/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$41/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($41/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$42/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($42/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Electron Damping Rates
unset label
ky=1
set origin kx*dx,ky*dy

set format x ''
set log y
set format y ''
set yrange [gam_min:gam_max]

plot \
file_in(mm) u 2:(-$36/sgn($5)) w l lc rgb base_color dt 1 title '' ,\
file_in(mm) u 2:($36/sgn($5)) w l lc rgb base_color dt 2 title '' ,\
file_in(mm) u 2:(-$47/sgn($5)) w l lc rgb n0_color dt 1 title '' ,\
file_in(mm) u 2:($47/sgn($5)) w l lc rgb n0_color dt 2 title '' ,\
file_in(mm) u 2:(-$48/sgn($5)) w l lc rgb n1_color dt 1 title '' ,\
file_in(mm) u 2:($48/sgn($5)) w l lc rgb n1_color dt 2 title ''

#Polarization
unset label
ky=0
set origin kx*dx,ky*dy

set format x '$10^{%L}$'

unset log y
set format y ''
set yrange [E_min:E_max]

set label '$k_\parallel d_p$' at graph 0.5,-0.25

plot \
file_in(mm) u ($2):(Ep($13,$14,$15,$16)) \
w l lc rgb base_color dt 1 title ''


unset multiplot