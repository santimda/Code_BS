set encoding utf8
set title '' #'Hydrodinamics according to V. code'
set key box outside opaque horizontal center top
set key samplen 0.7 spacing 1. font ",18"
set xrange [0.:3.14]
set yrange [0.:1.05]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
set ytics 0.2; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''

plot 'hydro.dat' u 1:2 t '{/Symbol r}/{/Symbol r}_{0}' w l lw 2.4 dt 3,\
     'hydro.dat' u 1:3 t 'p/p_0' w l lw 2.4 dt 4,\
     'hydro.dat' u 1:4 t 'vs/v_*' w l lw 2.4 lc rgb "coral" dt 5,\
     'hydro.dat' u 1:5 t 'c_{s}/v_*' w l lw 2.4 lc rgb "blue" dt 5

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'hydro.eps'
replot
!./fixbb hydro.eps