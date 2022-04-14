set encoding utf8
set title '' #'Comparison btw P, rho and v_t new vs old FS'
set key box opaque width 0.4 height 0.5 top right
set key samplen 0.7 spacing 1. font ",20" 
set xrange [0:3.]
set yrange [0:]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
set ytics 2.; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2 font ",18"
#set x2tics ('0' 0, '' 0.49, '1' 0.92, '' 1.27, '2' 1.54, '' 1.75, '3' 1.91, '' 2.04, '4' 2.15, '' 2.23, '5' 2.3)
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3) font ",16"
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'new_vs_old_FS.eps'

f(x) = 1

plot '../new_vs_old_RS.dat' u 1:3 t 'P' w l lw 2.4 lc 7,\
     '../new_vs_old_RS.dat' u 1:4 t '{/Symbol r}' w l lw 2.4 lc 3,\
     '../new_vs_old_RS.dat' u 1:2 t 'v_t' w l lw 2.4 lc 1,\
     f(x) w l t '' lc 'black' lw 2. dt 2


!./fixbb new_vs_old_FS.eps
!epstopdf new_vs_old_FS.eps