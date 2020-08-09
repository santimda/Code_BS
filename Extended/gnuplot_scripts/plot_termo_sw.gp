set encoding iso_8859_1
set title '' #'Thermodynamical quantities along the shocked wind'
set key box opaque width 0.95 height 1.1 center left
set key samplen 4. spacing 1.2 font ",20"
set xrange [0:2.3]
set yrange [0:1.05]
set xtics 0.5; set mxtics 5; set xtics nomirror 
set ytics 0.2; set mytics 2
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2
#set x2tics ('0' 0, '' 0.49, '1' 0.92, '' 1.27, '2' 1.54, '' 1.75, '3' 1.91, '' 2.04, '4' 2.15, '' 2.23, '5' 2.3)
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3)

plot '../termonorm_sw.dat' u 1:5 t 'T/T_0' w l lw 2.4 dt 1,\
     '../termonorm_sw.dat' u 1:6 t 'B/B_{eq_0}' w l lw 2.4 lc 7 dt 2,\
     '../termonorm_sw.dat' u 1:2 t '{/Symbol r}/{/Symbol r}_{0}' w l lw 2.4 dt 3,\
     '../termonorm_sw.dat' u 1:3 t 'P/P_0' w l lw 2.4 dt 4,\
     '../termonorm_sw.dat' u 1:4 t 'v_t /v_w' w l lw 2.4 lc rgb "coral" dt 5
     

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'termo_sw.eps'
replot
!./fixbb termo_sw.eps
