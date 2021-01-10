set encoding iso_8859_1
set key out
set key box opaque height 0.6 #width 0.50 
set key horiz samplen 3. spacing 1.01 font ",20"
set xrange [-6:-4]
set yrange [29.5:31.0]
set xtics 1; set mxtics 10
set ytics 0.5; set mytics 5
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L_{/Symbol e} [erg s^{-1}])'
set xlabel offset 0.0,0.3
set ylabel offset 1.,0.0

#-------------------------------------------------------------
# CONSTANTS:

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
#dist = 1.45*3.0856e21
dist = 1.72*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
a = 4.0*pi*dist**2
fac = 1.0e9*h/eV

#-------------------------------------------------------------
# FITTING :

x1=-5.0
g(x)=b0+b1*(x-x1)+b2*(x-x1)**2+b3*(x-x1)**3+b4*(x-x1)**4+b5*(x-x1)**5+b6*(x-x1)**6+b7*(x-x1)**7
fit[-6.1:-3.9] g(x) "../radsyn_sw_tot.dat" u (log10($1)):(log10($2)) via b0,b1,b2,b3,b4,b5,b6,b7
g1(x)=x>-3.?0 : x<-6.?0 : g(x)

#-----------------------------------------
#
#     PLOT
#
#----------------------------------------

# OJO QUE Ben2010 tomó una distancia de 1.4 kpc y yo ahora de 1.72 kpc!

plot log10(10**g1(x)) t 'Total' lw 3,\
     '../radsyn_sw_tot.dat' u (log10($1)):(log10($2)) t 'sync' w l lw 2.4 dt 2,\
     './data_radio2010.dat' u 1:2:3:4:5:6 w xyerrorbars pt 7 lc 7 t '2010',\
     './data_radio2020_2.3.dat' u 1:2:3:4:5:6 w xyerrorbars pt 6 lc 6 t '2020, 2.3',\
     './data_radio2020_1.5.dat' u 1:2:3:4:5:6 w xyerrorbars pt 4 lc 5 t '2020, 1.5'


set terminal postscript colour solid enhanced dl 1.25 lw 2.2 "Helvetica" 22
set output 'SED_radio.eps'
replot
!./fixbb SED_radio.eps
