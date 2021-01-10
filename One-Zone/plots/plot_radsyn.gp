set encoding iso_8859_1
set key top left
set xrange [-2:4]
set yrange [:]
set xtics 1
set mxtics 4
set xlabel 'log({/Symbol n} [GHz])'
set ylabel 'log(E L(E) [erg/s])'
h = 6.6260755e-27
eV = 1.602176e-12
f = eV/h

plot '../radsyn_FS.dat' u (log10(($1)*f/1.e9)):(log10($2)) t 'Sync + R-T' w l lw 3,\
	'../radsyn_FS.dat' u (log10(($1)*f/1.e9)):(log10($3)) t 'Sync' w l lw 2 

set terminal postscript colour solid enhanced dl 1. lw 2.8 "Helvetica" 20
set output 'radsyn.ps'
replot
!./fixbb radsyn.ps
