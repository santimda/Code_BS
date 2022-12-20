set encoding iso_8859_1
set size ratio 0.85
set title 'Shock width'
#set xtics 100
set xrange [0.:2.4]
set yrange [-0.1:0.5]
set ytics 1; set mytics 10
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
set ylabel 'H({/Symbol q})/R({/Symbol q})'
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
#set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'H.eps'

plot '../H_FS.dat' u 1:3 t 'Isothermal region' w l lw 2.4 dt 1,\
     '../H_FS.dat' u 1:4 t 'Hot layer' w l lw 2.4 dt 1,\
     '../H_RS.dat' u 1:3 t 'Reverse shock' w l lw 2.4 dt 1

!./fixbb H.eps
!epstopdf H.eps