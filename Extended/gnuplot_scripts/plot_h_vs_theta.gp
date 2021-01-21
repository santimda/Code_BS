set xrange [0:2]
set yrange [0:1]
set ytics 0.2
set mytics 2
set xtics 0.2
set mxtics 2
set key box opaque width 0.95 height 1.0 top left
set key samplen 3. spacing 1.1 font ",18"

plot '../Hnorm_sw.dat' u 1:2 w l lw 3 lc 1 t 'H/R' 

set terminal postscript color solid enhanced dl 2. lw 2.0 "Helvetica" 20
set output 'h_vs_theta.eps'
replot

!./fixbb h_vs_theta.eps