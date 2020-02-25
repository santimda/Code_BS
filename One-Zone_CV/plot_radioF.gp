set encoding iso_8859_1
set terminal postscript colour solid enhanced dl 1. lw 2.8 font "Helvetica,19"
set output 'radio_F.ps'

set multiplot
#####################################################
# LARGE PLOT
set key at -2.,-12.2
set grid
set font "Helvetica, 20"
set xrange [-6:5]
set yrange [-14.5:-12]
set xtics 3; set mxtics 3
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e}  F({/Symbol e}) [erg s^{-1} cm^{-2}])'
set pointsize 1.5


set object rectangle center -5.4,-13.8 size .9, .25
set arrow from -4.9,-13.8 to -1,-13.7 front lt -1 

dist = 2.9*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
a = 4.0*pi*dist**2.

plot  './data_GMRT_GWB.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*($2)*10**(-26))):(log10(1+($3/1000.0)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 lc 7 t 'GMRT(GWB)',\
	'./radsyn.dat' u (log10($1)):(log10($2/a)) t 'Model' w l lw 2 lc 1

#	'./radsyn_0.75.dat' u (log10($1)):(log10($2/a)) t '{/Symbol g}=0.75' w l lw 2 lc 2,\
#	'./radsyn_0.01.dat' u (log10($1)):(log10($2/a)) t '{/Symbol g}=0.01' w l lw 2 lc 3

#####################################################
# INSET PLOT
set origin .5, .2
set size .35,.35
clear
unset key
unset grid
unset object
unset arrow
unset xlabel
unset ylabel

set xtics .5; set mxtics 5; set xtics font "Helvetica,16"
set ytics .1; set mytics 1; set ytics font "Helvetica,16"
set bmargin 1
set tmargin 1
set lmargin 3
set rmargin 1
set xrange [-5.8:-5]
set yrange [-13.9:-13.7]

dist = 2.9*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
a = 4.0*pi*dist**2.

plot  './data_GMRT_GWB.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*($2)*10**(-26))):(log10(1+($3/1000.0)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 lc 7 t 'GMRT(GWB)',\
	'./radsyn.dat' u (log10($1)):(log10($2/a)) t 'Model' w l lw 2 lc 1

unset multiplot


!./fixbb radio_F.ps
!mv radio_F.ps radio_F.eps


