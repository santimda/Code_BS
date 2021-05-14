set encoding utf8
set title '' #'Total luminosity along the shocks in terms of L'
set key box top left opaque height 0.6
set key samplen 2.8 spacing 1. font ",18"
set xrange [0:3.14]
set yrange [0:1.8]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
set ytics 0.3; set mytics 3; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel 'L_{tot}/L_{sh}'; set ylabel offset 0.0,-0.2 font ",18"

plot '../Linj_RS.dat' u 1:5 t 'RS' w l lw 2.4 lc rgb '#a2142f',\
     '../Linj_FS.dat' u 1:5 t 'FS' w l lw 2.4 lc rgb '#aadc32'

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'L_inj_tot.eps'

replot
!./fixbb L_inj_tot.eps