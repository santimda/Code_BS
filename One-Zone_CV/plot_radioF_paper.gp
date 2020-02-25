set encoding iso_8859_1
set size square # ratio 0.8
set key top right 
set xrange [-6:4.1]
set yrange [-14.1:-12.5]
set xtics 3; set mxtics 3
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e}  F({/Symbol e}) [erg s^{-1} cm^{-2}])'
set pointsize 1.5

dist = 2.88*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
a = 4.0*pi*dist**2

plot  './data_GMRT_GWB.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*($2)*10**(-26))):(log10(1+($3/1000.0)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 lc -1 t 'uGMRT',\
	'./radsyn.dat' u (log10($1)):(log10($2/a)) t 'Model' w l lw 3 lc rgb "grey50"

set terminal postscript colour enhanced dl 1.25 lw 2. "Helvetica" 22
set output 'model_SED.ps'
replot
!./fixbb model_SED.ps
!mv model_SED.ps model_SED.eps


