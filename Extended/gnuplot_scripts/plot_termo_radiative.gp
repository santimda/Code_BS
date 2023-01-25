set encoding utf8
set size ratio 0.85
set title '' #'Thermodynamical quantities along forward shock'
#set key box opaque height 0.5
set key samplen 0.7 spacing 1.2 font ",18"
set xrange [0.:2.2]
#set yrange [0.:1.05]
set logscale y
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
#set ytics 0.2; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel ''
set x2label '{/Symbol D}s/R_0'; set x2label offset 0.0,-0.2 font ",18"
#set x2tics ('0' 0, '' 0.49, '1' 0.92, '' 1.27, '2' 1.54, '' 1.75, '3' 1.91, '' 2.04, '4' 2.15, '' 2.23, '5' 2.3)
set x2tics ('0' 0, '1' 0.92, '2' 1.54, '3' 1.91, '4' 2.15, '5' 2.3) font ",16"

set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'termo_radiative.eps'

plot '../termo_radiative_norm.dat' u 1:2 t '{/Symbol r}_{HL}/{/Symbol r}_{ISM}' w l lw 2.4 dt 3,\
     '../termo_radiative_norm.dat' u 1:3 t '{/Symbol r}_{IL}/{/Symbol r}_{ISM}' w l lw 2.4 dt 1,\
     '../termo_radiative_norm.dat' u 1:4 t 'T/T_{ISM}' w l lw 2.4 dt 1,\

!./fixbb termo_radiative.eps
!epstopdf termo_radiative.eps
