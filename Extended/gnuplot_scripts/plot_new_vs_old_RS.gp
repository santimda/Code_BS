set encoding utf8
set size ratio 0.8
set title '' #'Comparison btw P, rho and v_t new vs old RS'
set key box opaque width 0.4 height 0.5 top left
set key samplen 0.9 spacing 1.2 font ",20" 
set xrange [0.:3.]
set yrange [0.5:2.5]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
set ytics 1.; set mytics 10; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2 font ",18"
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3, '10' 2.63, '20' 2.84) font ",16"
#set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'new_vs_old_RS.eps'

f(x) = 1

plot '../new_vs_old_RS.dat' u 1:3 t 'P' w l lw 2.4 lc rgb 'forest-green',\
     '../new_vs_old_RS.dat' u 1:4 t '{/Symbol r}' w l lw 2.4 lc rgb 'dark-blue',\
     '../new_vs_old_RS.dat' u 1:2 t 'v_t' w l lw 2.4 lc rgb 'orange-red',\
     f(x) w l t '' lc 'black' lw 2. dt 2


!./fixbb new_vs_old_RS.eps
!epstopdf new_vs_old_RS.eps