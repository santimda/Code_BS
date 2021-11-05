set encoding utf8
set title '' #'Total luminosity along the shocks in terms of L'
set key top left opaque height 0.6
set key samplen 2.8 spacing 1. font ",18"
set xrange [0.0:3.141]
set yrange [0.0:]
set xtics 0.5; set mxtics 5; set xtics nomirror font ",16"
set ytics 5.; set mytics 2; set ytics font ",16"
set xlabel '{/Symbol q} [rad]'; set xlabel offset 0,0.25 font ",18"
set ylabel 'L_{tot}/L_{sh}'; set ylabel offset 0.0,-0.2 font ",18"


set terminal postscript color solid enhanced dl 1.25 lw 2.4 "Helvetica" 21
set output 'L_inj_tot.eps'

plot '../Linj_RS.dat' u 1:5 t 'RS' w l lw 3. lc "magenta" ,\
     '../Linj_FS.dat' u 1:5 t 'FS' w l lw 3. lc "blue"

!./fixbb L_inj_tot.eps
!epstopdf L_inj_tot.eps