set encoding utf8
set size ratio 0.8
set title
set key box opaque width 0.4 height 0.5 top left
set key samplen 0.7 spacing 1.2 font ",18"
set xrange [0.:]
set yrange [:]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
#set ytics 5.;
set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2 font ",18"
#set x2tics ('0' 0, '' 0.49, '1' 0.92, '' 1.27, '2' 1.54, '' 1.75, '3' 1.91, '' 2.04, '4' 2.15, '' 2.23, '5' 2.3)
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3, '10' 2.63, '20' 2.84) font ",16"
set label "NAME  FS" at graph 0.7,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'radadi_FS.eps'

f(x) = 0.

plot '../rad_adi_FS.dat' u 1:(log($2)) t 'FS' w l lc 1 lw 2.4 dt 1,\
     f(x) w l t '' lc 'black' lw 2. dt 2

!./fixbb radadi_FS.eps
!epstopdf radadi_FS.eps