set encoding utf8
set title '' #'Total luminosity along the shocks in terms of L'
set key box top left opaque height 0.6
set key samplen 2.8 spacing 1. font ",20"
set xrange [0:3.14]
set yrange [0:]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",18"
#set ytics 0.3; set mytics 3; set ytics font ",16"
set ytics 1; set mytics 4; set ytics font ",18"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",20"
set ylabel 'L_{tot}/L_{sh}'; set ylabel offset 0.0,-0.2 font ",20"


plot 'L_theta.dat' u 1:3 t 'sin^2' w l lw 2,\
	'L_theta.dat' u 1:4 t 'sin^3' w l lw 2


#plot 'L_theta.dat' u 1:4 t 'sin^3' w l lw 2

#'L_theta.dat' u 1:2 t 'sin' w l lw 2,\

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica,21"
set output 'L_theta.eps' 
replot