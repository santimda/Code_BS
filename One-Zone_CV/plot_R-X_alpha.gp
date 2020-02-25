set encoding iso_8859_1
set size square
set key top left spacing 1.1
set xrange [-6:5]
set yrange [-15:-10.5]
#set xtics 0.5; set mxtics 5
set xtics 2; set mxtics 2
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} F({/Symbol e}) [erg s^{-1} cm^{-2}])'; set ylabel offset 0.5,0
set pointsize 1.5

dist = 2.88*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
a = 4.0*pi*dist**2.

set object 1 rect from 2.7,-15 to 4,-10.5 fc rgb "grey" behind 
set obj 1 fillstyle noborder

plot  './data_GMRT_inner.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*($2)*10**(-26))):(log10(1+($3/1000.)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 lc 7 t 'GMRT',\
	'./radsyn_0.75_2.1.dat' u (log10($1)):(log10($2/a)) t '{/Symbol a}=2.1' w l lw 2 lc 1,\
	'./radsyn_0.75_2.0.dat' u (log10($1)):(log10($2/a)) t '{/Symbol a}=2.0' w l lw 2 lc 2,\
	'./radsyn_0.75_1.9.dat' u (log10($1)):(log10($2/a)) t '{/Symbol a}=1.9' w l lw 2 lc 3


set terminal postscript colour solid enhanced dl 1. lw 2.8 "Helvetica" 22
set output 'R-X_alpha.ps'
replot
!./fixbb R-X_alpha.ps
#!epstopdf R-X_alpha.ps
