set encoding utf8
set title '' #'Thermodynamical quantities along reverse shock'
set key box opaque width 0.4 height 0.5 center left
set key samplen 0.7 spacing 1. font ",18" 
set xrange [0:3.14]
set yrange [0:1.05]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
set ytics 0.2; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2 font ",18"
#set x2tics ('0' 0, '' 0.49, '1' 0.92, '' 1.27, '2' 1.54, '' 1.75, '3' 1.91, '' 2.04, '4' 2.15, '' 2.23, '5' 2.3)
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3) font ",16"
set label "NAME" at graph 0.4,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'termo_RS.eps'

plot '../termonorm_RS.dat' u 1:5 t 'T/T_0' w l lw 2.4 dt 1,\
     '../termonorm_RS.dat' u 1:6 t 'B/B_0' w l lw 2.4 lc 7 dt 2,\
     '../termonorm_RS.dat' u 1:2 t '{/Symbol r}/{/Symbol r}_{0}' w l lw 2.4 dt 3,\
     '../termonorm_RS.dat' u 1:3 t 'P/P_0' w l lw 2.4 dt 4,\
     '../termonorm_RS.dat' u 1:4 t 'v_t /v_w' w l lw 2.4 lc rgb "coral" dt 5,\
     '../termonorm_RS.dat' u 1:7 t 'v_{p} /v_w' w l lw 2.4 lc rgb "green" dt 5

!./fixbb termo_RS.eps
!epstopdf termo_RS.eps
