set encoding utf8
set title '' #'Compression factor and L_{inj} in RS'
set xrange [-0.1:]
set yrange [:6.]
set y2range [27.:]
set xtics 0.5; set mxtics 2; set ytics font ",16"
set ytics 0.5; set mytics 2; set ytics font ",16"
set y2tics 0.5; set y2tics font ",16"
set key bottom right
set key font ",18"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel 'Compression factor' font ",18"
set y2label 'log(L_{inj} [erg s^{-1}])'; set y2label offset 2.2,0 font",18"

plot '../termo_RS.dat' u 8:10 t 'r' w l lw 2.4 lc 8,\
     '../Linj_RS.dat' u 1:(log10($3)) axis x1y2 t 'Log(L_{inj})' wl lw 2.4 lc 7 

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica,22"
set output 'Linj_compres_RS.eps'
replot
!./fixbb Linj_compres_RS.eps