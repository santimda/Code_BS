set xrange [0:2]
set yrange [1:]
set ytics 0.2
set mytics 2
set xtics 0.2
set mxtics 2
set key box opaque width 0.95 height 1.0 top left
set key samplen 3. spacing 1.1 font ",18"

plot '../CD.dat' u 1:2 w l lw 3 lc 1 t 'alpha=0.01',\
	'../CD_aprox.dat' u 1:2 w l lw 3 lc 2 dt 2 t 'alpha=0'

set terminal postscript color solid enhanced dl 2. lw 2.0 "Helvetica" 20
set output 'r_vs_theta.eps'
replot

!./fixbb r_vs_theta.eps