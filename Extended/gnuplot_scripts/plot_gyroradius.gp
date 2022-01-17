set encoding iso_8859_1
#set title 'Gyroradius divided by shock width'
#set xtics 100
set xrange [0.:3.14]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
set ylabel 'log(H({/Symbol q})/r_{g}({/Symbol q}))'
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'rg_vs_H.eps'

plot '../Gyroradius_e_FS.dat' u 2:(log10($4)) t 'electrons FS' w l lw 2.4 dt 1,\
     '../Gyroradius_e_RS.dat' u 2:(log10($4)) t 'electrons RS' w l lw 2.4 dt 1,\
     '../Gyroradius_p_FS.dat' u 2:(log10($4)) t 'protons FS' w l lw 2.4 dt 1,\
     '../Gyroradius_p_RS.dat' u 2:(log10($4)) t 'protons RS' w l lw 2.4 dt 1

!./fixbb rg_vs_H.eps.eps
!epstopdf rg_vs_H.eps