set encoding iso_8859_1
set key bottom left
set xrange [2:5]
set yrange [-15:-12]
set xtics 0.5; set mxtics 5
set ytics 0.5; set mytics 5
set xlabel 'log_{10} (E [eV])'
set ylabel 'log_{10} (E F(E) [erg s^{-1} cm^{-2}])'
set pointsize 1.5

dist = 2.9*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
a = 4.0*pi*dist**2.

plot  './data_GMRT.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*($2)*10**(-26))):(log10(1+($3)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 lc 7 t 'GMRT',\
	'./sens_nustar.dat' u (log10($1)):(log10($2)) t 'NuSTAR' w l lt 1 lc rgb 'grey70' lw 5,\
	'./sens_xmm.dat' u (log10($1)):(log10($2)) t 'XMM' w l lt 1 lc rgb 'grey40' lw 5,\
	'./radsyn_lowB.dat' u (log10($1)):(log10($2/a)) t 'low B' w l lw 2 lc 1,\
	'./radic_lowB.dat' u (log10($1)):(log10($2/a)) t '' w l lw 2 lc 1 smooth bezier,\
	'./radsyn_highB.dat' u (log10($1)):(log10($2/a)) t 'high B' w l lw 2 lc 2,\
	'./radic_highB.dat' u (log10($1)):(log10($2/a)) t '' w l lw 2 lc 2 smooth bezier,\

#	'./radic_cmb.dat' u (log10($1)):(log10($2/a)) t 'IC-CMB' w l lw 2 smooth bezier,\
#	'./radbr.dat' u (log10($1)):(log10($2/a)) t 'rel. br.' w l lw 2 smooth bezier,\
#	'./radpp.dat' u (log10($1)):(log10($2/a)) t 'p-p' w l lw 2

set terminal postscript colour solid enhanced dl 1. lw 2.8 "Helvetica" 20
set output 'F_compare.ps'
replot
!./fixbb F_compare.ps
